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
#include "saveSolution.H"
//#include "ResidualConstruction.H"
//#include "ResidualConstruction2eqs.H"
//#include "NewtonRaphsonSolver2eqs.H"
//#include "ResidualConstruction2eqsTime.H"
#include "NewtonRaphsonSolver2eqs.H"
//#include "NewtonRaphsonSolver2eqstimeResidual.H"
//#include "NewtonRaphsonSolver.H"
//#include "NewtonRaphsonSolverWithTimeLoop.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Put this in a header/namespace you include once (e.g., namespace rom)
inline Foam::RectangularMatrix<Foam::scalar>
solveSmallLinearSPD(const Foam::RectangularMatrix<Foam::scalar>& A,
                    const Foam::RectangularMatrix<Foam::scalar>& b)
{
    using namespace Eigen;
    const Foam::label n = A.m();
    MatrixXd Ae(n,n); VectorXd be(n);
    for (Foam::label i=0; i<n; ++i){
        be[i] = b[i][0];
        for (Foam::label j=0; j<n; ++j) Ae(i,j) = A[i][j];
    }

    // Try LLT first (SPD)
    LLT<MatrixXd> llt(Ae);
    VectorXd xe;
    if (llt.info() == Success) {
        xe = llt.solve(be);
        if (llt.info() != Success) {
            Foam::Info<< "LLT failed at solve, falling back to LU\n";
            xe = Ae.fullPivLu().solve(be);
        }
    } else {
        Foam::Info<< "A not SPD (gauge issue?). Falling back to LU\n";
        xe = Ae.fullPivLu().solve(be);
    }

    // (Optional) warn if rank-deficient (Neumann/gauge not fixed)
    JacobiSVD<MatrixXd> svd(Ae);
    double rcond = svd.singularValues().tail(1)(0) / svd.singularValues()(0);
    if (rcond < 1e-12)
        Foam::Info<< "Warning: PPE matrix near-singular (rcond ~ " << rcond << ")\n";

    Foam::RectangularMatrix<Foam::scalar> x(n,1,0.0);
    for (Foam::label i=0; i<n; ++i) x[i][0] = xe[i];
    return x;
}


int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createParameters.H"

    //mathResidual mathRes(mesh, debug);
    IO io;
    utilities utils;
    saveSolution save(debug);
    Info << "Using galerkin.C inside galerkin projection folder" << endl;

    fileName offlinePath = runTime.path()/"offline";

    #include "createFields.H"

    //need to declare due to scope of if and due to the fact that ex takes ptr
    rom::Mat B, C, extra1, extra2, extra3, extra4, extra5, A, D, F, extra6, extra7, extra8, extra9, N, E1fo, F1fo, F2fo;
    rom::Expr ex, ex2;
    if(residual == "galerkin")
    {
        Info << "Using galerkin without stabilization" << endl;
        #include "createResidualGalerkin.H"
    }
    else if(residual == "galerkinPressure")
    {
        Info << "Using galerkin with pressure term but no stabilization " << endl;
        #include "createResidualGalerkinPressure.H"
    }
    else if(residual == "PPE")
    {
        Info << "using PPE to stabilize ROM " << endl;
        #include "createResidualPPE.H"
        //#include "createResidualPPEintByParts.H"
    }
    else
    {
        FatalErrorInFunction << "stabilization must be: galerkin, galerkinPressure or PPE " << exit(FatalError); 
    }

    #include "initialCondition.H"
    #include "solverOptions.H"
    
    /* rom::NewtontotalTime solver2(ex, a0, dt, tEnd, opt2, trj);
    solver2.solve();
    io.saveFoamMatrixToDisk<scalar>(timeCoeff,offlinePath/"timeCoeffNewton"); */

    if(residual == "PPE")
    {
        clockTime totalTime;
        clockTime solverTime;
        clockTime ioTime;
        scalar totalSolverTime = 0.0;
        scalar totalIOTime = 0.0;

        rom::NewtonSolver solver(ex, ex2, a0, b0, dt, tEnd, opt);
        utils.debugMsg(debug, "solver constructed");
        //Newton solver + saving solution
        runTime.setTime(startTime,0); 
        //saving 0th timestep
        /* if(startTime == 0)
        {
        save.saveSolutionToDisk<volVectorField>(uRom,uMean,uModes,a0);
        save.saveSolutionToDisk<volScalarField>(pRom,pMean,pModes,b0);
        uRom.write();
        pRom.write();
        utils.debugMsg(debug, "1st timestep saved");
        } */
        ioTime.timeIncrement();
        save.saveSolutionToDisk<volVectorField>(uRom,uMean,uModes,a0);
        save.saveSolutionToDisk<volScalarField>(pRom,pMean,pModes,b0);
        uRom.write();
        pRom.write();
        totalIOTime += ioTime.timeIncrement();
        utils.debugMsg(debug, "1st timestep saved");
        while(runTime.loop())
        {
            /* if(debug)
            {Info <<"timeloop at time " << runTime.value() << endl;} */
            //Info << runTime.value() << endl;
            solverTime.timeIncrement();
            auto [timeCoeff_a, timeCoeff_b] = solver.solve();
            totalSolverTime += solverTime.timeIncrement();

            static label iterCount = 0;
            iterCount++;

            if (debug && (iterCount <= 5 || iterCount % 10 == 0)) {
                Info << "\n=== Time step " << iterCount << " ===" << nl;
                
                // Check Jacobian
                rom::checkJacobian2(ex2, timeCoeff_a, timeCoeff_b);
                
                // Evaluate actual residuals
                auto [R2, scale2] = rom::residual(ex2, timeCoeff_a, timeCoeff_b);
                Info << "ex2 residual norm: " << rom::l2(R2) << nl;
                
                // If you have access to Newton iteration count
                // Info << "Newton iterations: " << solver.getIterCount() << nl;
            }

            //RectangularMatrix<scalar> timeCoeff_a = solver.solveInverse();
            //if(runTime.write())
            if(runTime.outputTime())
            {
                ioTime.timeIncrement();
                save.saveSolutionToDisk<volVectorField>(uRom, uMean, uModes, timeCoeff_a);
                save.saveSolutionToDisk<volScalarField>(pRom, pMean, pModes, timeCoeff_b);
                uRom.write();
                pRom.write();
                totalIOTime += ioTime.timeIncrement();
                //Info << "saving time" << runTime.value() << endl;
                // After reconstruction, check solution quality
                if (debug) {
                    // Check if reconstructed field satisfies conservation
                    scalar divU = fvc::div(uRom)().weightedAverage(mesh.V()).value();
                    Info << "Divergence of reconstructed velocity: " << divU << nl;
                    
                    // Check solution magnitude
                    Info << "Max |U|: " << max(mag(uRom)).value() << nl;
                    Info << "Max |p|: " << max(mag(pRom)).value() << nl;
                    
                    // Compare with reference if available
                    // if (refSolutionExists) { ... }
                    // After reconstruction, check pressure BCs
                    Info << "\n=== Pressure Boundary Conditions ===" << nl;
                    forAll(pRom.boundaryField(), patchI) {
                        Info << "Patch " << pRom.boundaryField()[patchI].patch().name() 
                            << ": type = " << pRom.boundaryField()[patchI].type() 
                            << ", mean = " << average(pRom.boundaryField()[patchI]) << nl;
                    }

                    // Compare gradients at boundaries
                    volVectorField gradPRom = fvc::grad(pRom);
                    Info << "Max |grad(p)| at boundaries:" << nl;
                    forAll(gradPRom.boundaryField(), patchI) {
                        Info << "  " << gradPRom.boundaryField()[patchI].patch().name() 
                            << ": " << max(mag(gradPRom.boundaryField()[patchI])) << nl;
                    }

                }
                //Info << "***NewtonLoop**** this is a at time t = " << runTime.value()  << " --> a = " << timeCoeff_a[0][0] << endl;
            }
            
        }
        Info << "\n=== Timing Summary ===" << nl;
        Info << "total solver + saving time " << totalTime.elapsedTime() << "s" <<endl;
        Info << "Total solver time: " << totalSolverTime << " s" << nl;
        Info << "Total I/O time: " << totalIOTime << " s" << nl;
        Info << "Total elapsed time: " << solverTime.elapsedTime() << " s" << nl;
        Info << "Solver efficiency: " << (totalSolverTime / solverTime.elapsedTime() * 100) << " %" << nl;
    }
    else
    {
        clockTime totalTime;
        clockTime solverTime;
        clockTime ioTime;
        scalar totalSolverTime = 0.0;
        scalar totalIOTime = 0.0;

        rom::NewtonSolver solver(ex, a0, dt, tEnd, opt);
        
        //Newton solver + saving solution
        runTime.setTime(startTime,0); 
        //saving 0th timestep
        /* if(startTime == 0)
        {
        save.saveSolutionToDisk<volVectorField>(uRom,uMean,uModes,a0);
        save.saveSolutionToDisk<volScalarField>(pRom,pMean,pModes,b0);
        uRom.write();
        pRom.write();
        utils.debugMsg(debug, "1st timestep saved");
        } */

        ioTime.timeIncrement();
        save.saveSolutionToDisk<volVectorField>(uRom,uMean,uModes,a0);
        save.saveSolutionToDisk<volScalarField>(pRom,pMean,pModes,b0);
        uRom.write();
        pRom.write();
        utils.debugMsg(debug, "1st timestep saved");
        totalIOTime += ioTime.timeIncrement();
        while(runTime.loop())
        {
            //Info <<"timeloop at time " << runTime.value() << endl;
            solverTime.timeIncrement();
            RectangularMatrix<scalar> timeCoeff_a = solver.solve().first;
            totalSolverTime += solverTime.timeIncrement();
            /* if(debug)
            {Info <<"timeloop at time " << runTime.value() << endl;} */
            //RectangularMatrix<scalar> timeCoeff_a = solver.solveInverse();
            //if(runTime.write())
            if(runTime.outputTime())
            {
                ioTime.timeIncrement();
                save.saveSolutionToDisk<volVectorField>(uRom, uMean, uModes, timeCoeff_a);
                save.saveSolutionToDisk<volScalarField>(pRom, pMean, pModes, timeCoeff_a);
                uRom.write();
                pRom.write();
                totalIOTime += ioTime.timeIncrement();
                //Info << "***NewtonLoop**** this is a at time t = " << runTime.value()  << " --> a = " << timeCoeff_a[0][0] << endl;
            }
            
        }
        Info << "\n=== Timing Summary ===" << nl;
        Info << "total solver + saving time " << totalTime.elapsedTime() << "s" <<endl;
        Info << "Total solver time: " << totalSolverTime << " s" << nl;
        Info << "Total I/O time: " << totalIOTime << " s" << nl;
        Info << "Total elapsed time: " << solverTime.elapsedTime() << " s" << nl;
        Info << "Solver efficiency: " << (totalSolverTime / solverTime.elapsedTime() * 100) << " %" << nl;
        Info << "Mean of pressure mode 0: " << pMeanMode0 << nl;
    }

    
    return 0;
}



   // ************************************************************************* //


