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
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    /***********************************/
    /**********CLASS INITIALIZ**********/
    /***********************************/
    utilities utils; //helper class
    IO IO; 
    math math(mesh);
    Info << "Using POD.C inside PODpost folder" << endl;

    #include "createParameters.H" //helper that reads dictionaries and saves datas 
    Info.level = infoLevel; //change in debugdict to have info of solver or not
    utils.debugMsg("starting parametric simulation with " + parameterSpace + " as parameter", debug); //function that prints out only if debug is on
    IO.createFolderEnvironment("POD", mesh, runTime, debug);

    /*****************************************************/
    /*** initialize the snapshots lists for U, p, phi ***/
    DynamicList<vectorField> uSnap;      //velocity
    DynamicList<scalarField> pSnap;     //pressure
    DynamicList<scalarField> fSnap;
    DynamicList<scalar>      uBC;
    DynamicList<scalar>      fBC;
    
    if(parSize > 1)
    {
        for(label i = 1; i <= parSize; ++i)
        {

            utils.debugMsg(debug, "getting snapshots for parameter ", i);
            runTime.setTime(startTime,0);
            Info << "time is: " << runTime.timeName() << endl;
            volScalarField pSnapField
            (
                IOobject
                (
                    "p" + name(i),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            volVectorField USnapField
            (
                IOobject
                (
                    "U" + name(i),
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            surfaceScalarField fSnapField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(USnapField) & mesh.Sf()
            );
            pSnap.append(pSnapField.internalField());
            uSnap.append(USnapField.internalField());
            fSnap.append(fSnapField.internalField());

            while(runTime.loop())
            {
                if(runTime.outputTime())
                {
                    Info << "time is: " << runTime.timeName() << endl;
                    volScalarField pSnapField
                    (
                        IOobject
                        (
                            "p" + name(i),
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    volVectorField USnapField
                    (
                        IOobject
                        (
                            "U" + name(i),
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    surfaceScalarField fSnapField
                    (
                        IOobject
                        (
                            "phi",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        fvc::interpolate(USnapField) & mesh.Sf()
                    );
                    pSnap.append(pSnapField.internalField());
                    uSnap.append(USnapField.internalField());
                    fSnap.append(fSnapField.internalField());
                }
            }
        }
    }
    else
    {
            runTime.setTime(startTime,0);
            Info << "Not parametric...time is: " << runTime.timeName() << endl;
            volScalarField pSnapField
            (
                IOobject
                (
                    "p",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            volVectorField USnapField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            surfaceScalarField fSnapField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::interpolate(USnapField) & mesh.Sf()
            );
            pSnap.append(pSnapField.internalField());
            uSnap.append(USnapField.internalField());
            fSnap.append(fSnapField.internalField());

            while(runTime.loop())
            {
                if(runTime.outputTime())
                {
                    Info << "time is: " << runTime.timeName() << endl;
                    volScalarField pSnapField
                    (
                        IOobject
                        (
                            "p",
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    volVectorField USnapField
                    (
                        IOobject
                        (
                            "U" ,
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    surfaceScalarField fSnapField
                    (
                        IOobject
                        (
                            "phi",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::NO_WRITE
                        ),
                        fvc::interpolate(USnapField) & mesh.Sf()
                    );
                    pSnap.append(pSnapField.internalField());
                    uSnap.append(USnapField.internalField());
                    fSnap.append(fSnapField.internalField());
                }
            }
    }
        /*********************************************************************************/
    /*********************** STARTING POD COMPUTATION *******************************/
    /*******************************************************************************/



    //convert snpashots lists to eigen in order to do POD math
    utils.debugMsg("converting snap to eigen", debug);
    Eigen::MatrixXd uSnapEigen = IO.foamListToEigenMatrix(uSnap, debug, "USnap");
    Eigen::MatrixXd pSnapEigen = IO.foamListToEigenMatrix(pSnap, debug, "pSnap");
    Eigen::MatrixXd fSnapEigen = IO.foamListToEigenMatrix(fSnap, debug, "fSnap");
    Eigen::VectorXd uBCEigen   = IO.foamListToEigenVector(uBC, debug, "uBC");
    Eigen::VectorXd fBCEigen   = IO.foamListToEigenVector(fBC, debug, "fBC");  

    //taking mean field and if necessary scaling it with boundary value to create the lifitng function
    utils.debugMsg("taking mean field", debug);
    Eigen::VectorXd uMeanEigen      = math.arithmeticMeanRowWise(uSnapEigen, debug, "U");
    Eigen::VectorXd pMeanEigen      = math.arithmeticMeanRowWise(pSnapEigen, debug, "p");
    /* Eigen::VectorXd uMeanEigen = uSnapEigen.col(0);
    Eigen::VectorXd pMeanEigen = pSnapEigen.col(0); */
    /* label Np = pSnap[0].size();
    Eigen::VectorXd pMeanEigen =  Eigen::VectorXd::Zero(Np); 
 */
    //saving mean in foam format and store it in folder 
    IO.EigenVecWriteToDisk<volVectorField>(uMeanEigen, mesh, runTime, "U","uMean", "POD/Mean",debug, "velocity");
    IO.EigenVecWriteToDisk<volScalarField>(pMeanEigen, mesh, runTime, "p","pMean", "POD/Mean",debug, "pressure");
    Info << "Restoring runTime after saving mean to: " << runTime.timeName() << nl; 

    //to get the mean of the flux, the mean of u is used
    volVectorField uMeanField
    (
        IOobject("uMean", "POD/Mean", mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh
    );
 
    surfaceScalarField fMeanField
    (
        IOobject("fMean", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), 
        fvc::interpolate(uMeanField) & mesh.Sf()
    );
    const label oldIndex = runTime.timeIndex();
    const instant oldInst(runTime.time().value(), runTime.timeName());
    //trick runtime to go in relpath
    runTime.setTime(instant(runTime.time().value(),"POD/Mean"), oldIndex);
    fMeanField.write();
    // Restore runtime I had
    runTime.setTime(oldInst, oldIndex);
    
    DynamicList<scalarField> fMeanList;
    fMeanList.append(fMeanField.internalField());
    Eigen::VectorXd fMeanEigen = IO.foamListFieldToEigenVector(fMeanList, debug, "fMean");
    

    //creating the fluctuation field with homogen BCs 
    utils.debugMsg("compute fluctuations", debug);
    uSnapEigen = math.computeFluct(uSnapEigen, uMeanEigen, debug, "U");
    fSnapEigen = math.computeFluct(fSnapEigen, fMeanEigen, debug, "f");
    pSnapEigen = math.computeFluct(pSnapEigen, pMeanEigen, debug, "p");

    //create correlation matrix to do the snapshots method SVD (only for p and u) 1/Ns*<u_i, u_j>_L2
    utils.debugMsg("compute correlation matrix", debug);
    Eigen::MatrixXd Cu = math.correlationMatrix2(uSnapEigen,debug, "velocity");
    Eigen::MatrixXd Cp = math.correlationMatrix2(pSnapEigen, debug, "pressure");

    //compute SVD U = USiV^T and extrapolate eigenVectors (columns of matrix V) and eigenvalues 
    //already truncated here 
    utils.debugMsg("computing SVD", debug);
    auto [EWu, EVu] =  math.SVDTruncation(Cu, nModesU, debug, "velocity"); // --> a little bit less memory usage (you will get exact same modes as complete SVD)
    auto [EWp, EVp] =  math.SVDTruncation(Cp, nModesp, debug, "pressure");

    // === DIAGNOSTICS: singular/eigen values, energy, conditioning ===
    auto printSV = [&](const Eigen::VectorXd& EW, const std::string& tag)
    {
        const int m = static_cast<int>(EW.size());
        double maxSV = EW(0);
        double minSV = EW(m-1);
        double cond  = (minSV > 0) ? maxSV/minSV : std::numeric_limits<double>::infinity();

        Info<< "\n[" << tag << "] singular/eigen values (first 10 or all):" << nl;
        for (int i=0; i<std::min(m,10); ++i)
            Info<< "  s["<<i<<"] = " << EW(i) << nl;

        // cumulative energy
        double tot = EW.sum();
        double cum = 0.0;
        int modes_99 = m;  // modes needed for 99% energy
        bool found_99 = false;
        
        Info<< "["<<tag<<"] cumulative energy:" << nl;
        for (int i=0; i<m; ++i)
        {
            cum += EW(i);
            double energy_fraction = (tot > 0) ? cum/tot : 0;
            
            if (i < 30 || i == m-1)
                Info<< "  k="<<i+1<<"  E="<< energy_fraction << nl;
            
            // find first mode that reaches 99% energy
            if (!found_99 && energy_fraction >= thresholdEnergy)
            {
                modes_99 = i + 1;
                found_99 = true;
            }
        }

        Info<< "["<<tag<<"] s_max="<<maxSV<<" s_min="<<minSV<<"  cond="<<cond<<nl;

        // choose a safe numerical rank (optional suggestion)
        const double eps = std::numeric_limits<double>::epsilon();
        const double tol = maxSV * eps * m;
        int r = 0; while (r<m && EW(r) > tol) ++r;
        //Info<< "["<<tag<<"] suggested numerical rank r="<<r<<" with tol="<<tol<<nl;
        
        // report 99% energy threshold
        Info<< "["<<tag<<"] modes for "<< thresholdEnergy*100 << "% energy: "<<modes_99<<nl;
    };

    printSV(EWu, "U");
    printSV(EWp, "p");

    // extrapolates and saves modes 
    utils.debugMsg("extrapolating modes", debug);
    Eigen::MatrixXd uModesEigen = math.getModes2(EVu, EWu, uSnapEigen, debug, "velocity");
    Eigen::MatrixXd fModesEigen = math.getModes2(EVu, EWu, fSnapEigen, debug, "flux"); // this is not orthornmal by construction since using POD of velocity (but it doesnt have to be)
    Eigen::MatrixXd pModesEigen = math.getModes2(EVu, EWu, pSnapEigen, debug, "pressure");
    Eigen::MatrixXd pModesEigen2 = math.getModes2(EVp, EWp, pSnapEigen, debug, "pressure");

    Info << "modes computed... saving to disk..." << nl;
    
    //save POD modes for p, U, phi in foam format to file and sets BCs to be homogenous
    IO.EigenMatrixWriteToDiskBCchange<volVectorField>(uModesEigen, vector(0,0,0), mesh, runTime, "U", "UMode", "POD", 1, 1, debug, "velocity");
    IO.EigenMatrixWriteToDiskBCchange<surfaceScalarField>(fModesEigen, 0.0, mesh, runTime, "phi", "fMode", "POD", 1, 1, debug, "flux");
    IO.EigenMatrixWriteToDiskBCchange<volScalarField>(pModesEigen,0.0, mesh, runTime, "p", "pMode", "POD", 1, 1, debug, "pressure");
    IO.EigenMatrixWriteToDiskBCchange<volScalarField>(pModesEigen2,0.0, mesh, runTime, "p", "pMode2", "POD", 1, 1, debug, "pressure2");

    //reconstructing solution
    Eigen::MatrixXd UFluctRecon = math.reconstructSolution(uModesEigen, EWu.cwiseSqrt(), EVu);
    Eigen::MatrixXd pFluctRecon = math.reconstructSolution(pModesEigen2, EWp.cwiseSqrt(), EVp); 
    //std::cout << "fluct field recon constructed" << std::endl;
    Eigen::MatrixXd Urecon = UFluctRecon.colwise() + uMeanEigen;
    Eigen::MatrixXd precon = pFluctRecon.colwise() + pMeanEigen;
    //std::cout << "sol recon constructed" << std::endl;
    IO.EigenMatrixWriteToDiskRunTime<volVectorField>(Urecon, mesh, runTime,startTime, "U", "UReconPost", debug, "Ureconstruction Saved to file");
    IO.EigenMatrixWriteToDiskRunTime<volScalarField>(precon, mesh, runTime,startTime, "p", "pReconPost", debug, "preconstruction Saved to file");
    return 0;
}


// ************************************************************************* //
