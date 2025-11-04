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
#include "mathError.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createParameters.H"

    //mathResidual mathRes(mesh, debug);
    IO io;
    utilities utils;
    mathError math;

    const instantList& times = utils.getTimeDirs(runTime,true,debug);
    label Nt = times.size();
    label Ns = ROMfieldU.size();
    Info << times << endl;

    Eigen::MatrixXd errVelocity = Eigen::MatrixXd::Zero(Nt+1, Ns);
    Eigen::MatrixXd errPressure = Eigen::MatrixXd::Zero(Nt+1, Ns);
    Eigen::MatrixXd errGradPressure = Eigen::MatrixXd::Zero(Nt+1, Ns);
    
    runTime.setTime(startTime,0);
    label ti = 1;
    for(label i = 0; i < Ns; ++i)
    {
    #include "createInitFields.H"
    scalar errorVelocity = math.errorVector(u0Fom,u0Diff);
    scalar errorGradPressure = math.errorVector(p0Grad,p0DiffGrad);
    scalar errorPressure = math.errorScalar(p0Fom,p0Diff);
    if (!std::isfinite(errorVelocity)) errorVelocity = 0;
    if (!std::isfinite(errorPressure)) errorPressure = 0;
    if (!std::isfinite(errorGradPressure)) errorGradPressure = 0;
    errVelocity(0,i) = errorVelocity * 100;
    errPressure(0,i) = errorPressure * 100;
    errGradPressure(0,i) = errorGradPressure * 100;
    }
    while(runTime.loop())
    {
        if(runTime.outputTime())
        {
            for(label i = 0; i < Ns; ++i)
            {
            #include "createFields.H"
            scalar errorVelocity = math.errorVector(uFom,uDiff);
            scalar errorGradPressure = math.errorVector(pGrad,pDiffGrad);
            scalar errorPressure = math.errorScalar(pFom,pDiff);
            if (!std::isfinite(errorVelocity)) errorVelocity = 0;
            if (!std::isfinite(errorPressure)) errorPressure = 0;
            if (!std::isfinite(errorGradPressure)) errorGradPressure = 0;

            errVelocity(ti,i) = errorVelocity * 100;
            errPressure(ti,i) = errorPressure * 100;
            errGradPressure(ti,i) = errorGradPressure * 100;
            }
            //Info << runTime.timeName() << endl;
            ++ti;
        }
    }
    errVelocity.conservativeResize(ti, Ns);
    errPressure.conservativeResize(ti, Ns);
    errGradPressure.conservativeResize(ti, Ns);
    scalar saveTime = dt * writeInterval;
    std::cout<< errVelocity << std::endl;
    std::cout << errPressure << std::endl;
    utils.plotVsTime2(errVelocity,startTime,saveTime,namesU,"time [s]","Error U [%]","error U vs time");
    utils.plotVsTime2(errPressure,startTime,saveTime,namesP,"time","Error p[%]","error p vs time");
    utils.plotVsTime2(errGradPressure,startTime,saveTime,namesP,"time[s]","Error grad p[%]","error grad p vs time");

    return 0;
}



   // ************************************************************************* //


