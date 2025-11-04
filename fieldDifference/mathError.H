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

#ifndef MATHERROR_H
#define MATHERROR_H

#include "fvCFD.H"
#include <Eigen/Dense>

namespace Foam
{
class mathError
{
    public:
        mathError() = default;                         // Constructor
        mathError(const mathError&) = delete;          // no copy constructor
        void operator=(const mathError&) = delete;     // no assignment
        virtual ~mathError() = default;                // default destructor


        /****************************************/
        /*********** MEMBER FUNCTIONS ***********/
        /****************************************/
       scalar L2prodVector(const volVectorField& field)
       {
            const auto& contraction = field & field;
            return fvc::domainIntegrate(contraction).value();
       }

       scalar L2prodScalar(const volScalarField& field)
       {
            const auto& contraction = field * field;
            return fvc::domainIntegrate(contraction).value();
       }

       scalar errorVector(const volVectorField& fom, const volVectorField& diff)
       {
            scalar num = L2prodVector(diff);
            scalar den = L2prodVector(fom);
            scalar error  = std::sqrt(num/den);
            return error;
       }

       scalar errorScalar(const volScalarField& fom, const volScalarField& diff)
       {
            scalar num = L2prodScalar(diff);
            scalar den = L2prodScalar(fom);
            scalar error  = std::sqrt(num/den);
            return error;
       }

};
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
