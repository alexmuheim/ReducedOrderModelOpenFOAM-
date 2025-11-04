/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    utilitiesROM

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "RectangularMatrix.H"
#include <Eigen/Dense>
#include <utility> 
#include <chrono>
#include "IOhelper.H"
#include "MathOffline2.H"
#include "MathOffline3.H"


// Return true if this patch type should contribute to boundary integrals
inline bool isPhysicalPatchType(const word& t)
{
    // extend as needed (eg add "wedge", "symmetryPlane", "cyclic", "processor", …)
    return !(t == "empty" || t == "processor" || t == "processorCyclic"
          || t == "cyclic" || t == "cyclicAMI"
          || t == "symmetry" || t == "symmetryPlane" || t == "wedge");
}

// Build a patch-id set (auto-skips non-physical like 'empty')
inline labelHashSet makePhysicalPatchSet
(
    const fvMesh& mesh,
    const UList<word>& patchNames = UList<word>()   // optional explicit list
)
{
    const polyBoundaryMesh& bnd = mesh.boundaryMesh();
    labelHashSet ids;

    if (patchNames.empty())
    {
        // include all physical patches
        forAll(bnd, p)
            if (isPhysicalPatchType(bnd[p].type()))
                ids.insert(p);
    }
    else
    {
        // include only requested names that are also physical
        forAll(patchNames, k)
        {
            const label pid = bnd.findPatchID(patchNames[k]);
            if (pid >= 0 && isPhysicalPatchType(bnd[pid].type()))
                ids.insert(pid);
        }
    }
    return ids;
}

// Sum a surfaceScalarField over selected boundary patches
inline scalar sumOnPatches(const surfaceScalarField& s, const labelHashSet& patchIDs)
{
    scalar acc = 0.0;
    const auto& sbf = s.boundaryField();
    forAllConstIter(labelHashSet, patchIDs, it)
    {
        const label pid = it.key();
        acc += gSum(sbf[pid]);   // MPI-safe reduction
    }
    return acc;
}

// N(i,j) = ∫_Γ (n × ∇χ_i) · (curl φ_j) dS   (auto-skips 'empty' etc.)
inline RectangularMatrix<scalar>
assembleNij(const PtrList<volScalarField>& pModes,        // χ_i
            const PtrList<volVectorField>& uModes,        // φ_j
            const UList<word>& patchNames = UList<word>())// optional subset Γ
{
    const fvMesh& mesh = uModes[0].mesh();
    const label Np = pModes.size();
    const label Nu = uModes.size();

    RectangularMatrix<scalar> N(Np, Nu, 0.0);

    // face unit normals
    surfaceVectorField n(mesh.Sf()/mesh.magSf());

    // physical patches only (empty/cyclic/processor/symmetry/wedge are excluded)
    const labelHashSet patchIDs = makePhysicalPatchSet(mesh, patchNames);

    for (label i=0; i<Np; ++i)
    {
        tmp<volVectorField> tGradPi = fvc::grad(pModes[i]);
        surfaceVectorField gradPi_f = fvc::interpolate(tGradPi());
        surfaceVectorField nXgradPi = n ^ gradPi_f;

        for (label j=0; j<Nu; ++j)
        {
            tmp<volVectorField> tCurlUj = fvc::curl(uModes[j]);
            surfaceVectorField curlUj_f = fvc::interpolate(tCurlUj());

            surfaceScalarField integrand = (nXgradPi & curlUj_f) * mesh.magSf();
            N(i,j) = sumOnPatches(integrand, patchIDs);
        }
    }
    return N;
}



int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    IO IO;
    utilities utils;

    #include "createParameters2.H"
    #include "createFields.H"

    mathOffline mathOff(mesh, debug);
    mathOffline2 mathOff2(mesh, debug);
    fileName offlinePath = runTime.path()/"offline";

    Info << "Using offlineROM2.C code" << endl;

    clockTime computeTimeOffline;
    if (isDir(offlinePath))
    {
        if (!fileHandler().rmDir(offlinePath, /*recursive*/ true))
        {
            FatalErrorInFunction << "Failed to remove " << offlinePath << exit(FatalError);
        }
    }
    mkDir(offlinePath);

    //residual creation
    if(residual == "galerkin")
    {
        Info << "using galerkin without any stabilization and without pressure term" << endl;
        #include "residualGalerkin/createResidualGalerkin.H"
    }
    else if(residual == "galerkinPressure")
    {
        Info << "using galerkin without any stabilization but with pressure term" << endl;
        #include "residualGalerkin/createResidualGalerkinPressure.H"
    }
    else if(residual == "PPE")
    {
        Info << "using galerkin with PPE stabilization" << endl;
        //#include "residualGalerkin/createResidualGalerkinPPE.H"
        //#include "residualGalerkin/createResidualGalerkinPPE2.H"
        #include "residualGalerkin/createResidualGalerkinFirstOrder.H"
    }
    else
    {
        FatalErrorInFunction << "stabilization must be: galerkin, galerkinPressure " << exit(FatalError); 
    }


    return 0;
}



   // ************************************************************************* //


