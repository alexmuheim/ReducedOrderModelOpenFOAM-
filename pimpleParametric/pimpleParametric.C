/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include <chrono>
#include "Utilities.H"
#include "IOhelper.H"
#include "MathROM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "parametric Transient solver using code pimpleParametric.C inside folder pimpleParametric"
    );
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    /***********************************/
    /**********CLASS INITIALIZ**********/
    /***********************************/
    utilities utils; //helper class


    #include "createParameters.H" //helper that reads dictionaries and saves datas 
    Info.level = infoLevel; //change in debugdict to have info of solver or not
    utils.debugMsg("starting parametric simulation with " + parameterSpace + " as parameter", debug); //function that prints out only if debug is on

    DynamicList<scalar>      uBC;
    DynamicList<scalar>      fBC;

    /****************************************/
    /*********** PARAMETER LOOP  ************/
    /****************************************/
    for(label par = 0; par < parSize; ++par)
    {  
        //resetting time each simulation
        runTime.setTime(0,0); 

        //changing viscosity if the parameter chosen is viscosity
        if(parameterSpace == "viscosity")
        {
            transportProperties.set("nu", nuList[par]);
            transportProperties.regIOobject::write();
        }
        
        #include "postProcess.H"
        #include "addCheckCaseOptions.H"
        #include "initContinuityErrs.H"
        #include "createDyMControls.H"
        #include "createFields.H"
        #include "createUfIfPresent.H"
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"

        //changing patch value if parameter space is dirichlet
        label parPatchID = mesh.boundary().findPatchID(parPatch);
        if(parameterSpace == "dirichlet")
        {
            vector newVelocityBC = dirichletList[par];
            // Loop through all boundary patches
            if(parPatchID == -1)
            {
                FatalError << "Patch '" << parPatch << "' not found in mesh!" 
                        << "\nAvailable patches are: " << mesh.boundaryMesh().names()
                        << exit(FatalError);
            }
            if(parPatchID != -1)
            {
                auto& pf = U.boundaryFieldRef()[parPatchID];
                if(isType<fixedValueFvPatchVectorField>(pf))
                {
                    refCast<fixedValueFvPatchVectorField>(pf) == newVelocityBC;
                }
                else
                {
                    FatalError << "Patch '" << parPatch << "' is not a fixedValue boundary condition!"
                            << "\nPatch type is: " << U.boundaryField()[parPatchID].type()
                            << exit(FatalError);
                }
            }
        }
        

        turbulence->validate();

        if (!LTS)
        {
            #include "CourantNo.H"
            #include "setInitialDeltaT.H"
        }

        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        
        Info<< "\nStarting time loop\n" << endl;
        
        while (runTime.run())
        {
            #include "readDyMControls.H"

            if (LTS)
            {
                #include "setRDeltaT.H"
            }
            else
            {
                #include "CourantNo.H"
                #include "setDeltaT.H"
            }

            ++runTime;

            Info<< "Time = " << runTime.timeName() << nl << endl;

            // --- Pressure-velocity PIMPLE corrector loop
            while (pimple.loop())
            {
                if (pimple.firstIter() || moveMeshOuterCorrectors)
                {
                    // Do any mesh changes
                    mesh.controlledUpdate();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            // Calculate absolute flux
                            // from the mapped surface velocity
                            phi = mesh.Sf() & Uf();

                            #include "correctPhi.H"

                            // Make the flux relative to the mesh motion
                            fvc::makeRelative(phi, U);
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                #include "UEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    laminarTransport.correct();
                    turbulence->correct();
                }
            }
            
            if(runTime.outputTime())
                {   
                    uBC.append(U.boundaryField()[parPatchID][0].x());
                    fBC.append(phi.boundaryField()[parPatchID][0]);
                    volVectorField USave
                    (
                        IOobject("U" + name(par+1), runTime.timeName(), mesh,
                                IOobject::NO_READ, IOobject::AUTO_WRITE),
                        U
                    );
                    volScalarField pSave
                    (
                        IOobject("p" +name(par+1), runTime.timeName(), mesh,
                                IOobject::NO_READ, IOobject::AUTO_WRITE),
                        p
                    );
                    USave.write();
                    pSave.write();
                }

            runTime.write(); 
            runTime.printExecutionTime(Info);
        }
    }

    //reconstruction
    return 0;
}


// ************************************************************************* //
