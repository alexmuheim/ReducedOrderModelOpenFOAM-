/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

Class
    POD

Description
    Compute the SVD 

SourceFiles
    POD.C

\*---------------------------------------------------------------------------*/

#ifndef SAVESOLUTION_H
#define SAVESOLUTION_H
#pragma once
#include "fvCFD.H"
#include <Eigen/Dense>
#include "Utilities.H"

namespace Foam
{
    class saveSolution: public utilities
    {
        private: 
            bool debug;
        public:
            explicit saveSolution(bool debug) 
            : debug(debug)
            {}                         // Constructor
            saveSolution(const saveSolution&) = delete;          // no copy constructor
            void operator=(const saveSolution&) = delete;     // no assignment
            virtual ~saveSolution() = default;                // default destructor

            /****************************************/
            /*********** MEMBER FUNCTIONS ***********/
            /****************************************/
            template<typename T>
            void saveSolutionToDisk(T& field, const T& meanField, const PtrList<T>& modes, const RectangularMatrix<scalar>& timeCoeff, word key = "")
            {
                label nModes = modes.size();
                auto& internalField = field.internalFieldRef();
                internalField = meanField.internalField();
                for(label m = 0; m < nModes; ++m)
                {
                    const scalar coeff = timeCoeff[m][0];
                    internalField += coeff * modes[m].internalField();
                }
                field.correctBoundaryConditions();
            }


            template<typename T>
            void saveSolutionToDiskPressure(T& field, const T& meanField, const PtrList<T>& modes, const RectangularMatrix<scalar>& timeCoeff, word key = "")
            {
                label nModes = modes.size();
                auto& internalField = field.internalFieldRef();
                internalField = meanField.internalField();
                for(label m = 0; m < nModes; ++m)
                {
                    const scalar coeff = timeCoeff[m][0];
                    internalField += coeff * modes[m].internalField();
                }
                // Instead of using cell[0], compute the average or use a fixed reference
                scalar pRef = gAverage(internalField);  // or gMin, or a fixed value
                Info << "Adjusting pressure by reference value: " << pRef << endl;

                forAll(internalField, i)
                {
                    internalField[i] -= pRef;
                }
                field.correctBoundaryConditions();
            }
            

    }; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
