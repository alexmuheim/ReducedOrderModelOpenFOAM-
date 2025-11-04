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

#ifndef MATHROM_H
#define MATHROM_H

#include "fvCFD.H"
#include <Eigen/Dense>
#include "Utilities.H"

namespace Foam
{
class mathOffline: public utilities
{
    private: const fvMesh& mesh_; //reference to mesh
             bool debug;
    public:
        explicit mathOffline(const fvMesh& mesh, bool debug) 
        : mesh_(mesh),
        debug(debug)
        {}                         // Constructor
        mathOffline(const mathOffline&) = delete;          // no copy constructor
        void operator=(const mathOffline&) = delete;     // no assignment
        virtual ~mathOffline() = default;                // default destructor

        /****************************************/
        /*********** MEMBER FUNCTIONS ***********/
        /****************************************/

        PtrList<volVectorField> fieldGradient(const PtrList<volScalarField>& fields, const word& key=word::null)
        {
            PtrList<volVectorField> grads(fields.size());
            forAll(fields, i)
            {
                tmp<volVectorField> tg = fvc::grad(fields[i]);
                grads.set(i, new volVectorField(tg));  // materialize tmp into heap object
                if (debug)
                    Info<< "grad(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return grads;
        }

        /* template<typename T>
        PtrList<T> fieldLaplacian(const PtrList<T>& fields, const word& key=word::null)
        {
            PtrList<T> lapl(fields.size());

            forAll(fields, i)
            {
                tmp<T> tl = fvc::laplacian(fields[i]);
                lapl.set(i, new T(tl));  // materialize tmp into heap object
                if (debug)
                    Info<< "laplacian(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return lapl;
        } */

        //prova nuova con div(grad) al posto di laplacian --> consistent
        template<typename T>
        PtrList<T> fieldLaplacian(const PtrList<T>& fields, const word& key=word::null)
        {
            PtrList<T> lapl(fields.size());

            forAll(fields, i)
            {
                auto tg = fvc::grad(fields[i]);
                tmp<T> tl = fvc::div(tg());
                lapl.set(i, new T(tl));  // materialize tmp into heap object
                if (debug)
                    Info<< "laplacian(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return lapl;
        }

        template<typename T1, typename T2>
        PtrList<T1> fieldDiv(const PtrList<T2>& fields, const word& key=word::null)
        {
            PtrList<T1> divs(fields.size());
            forAll(fields, i)
            {
                tmp<T1> td = fvc::div(fields[i]);
                divs.set(i, new T1(td));  // materialize tmp into heap object
                if (debug)
                    Info<< "div(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return divs;
        }
        
        //iterates over flux field
        PtrList<volVectorField> convectiveTermFlux(const PtrList<volVectorField>& vecFields, const PtrList<surfaceScalarField>& surfFields, const word& key=word::null)
        {
            PtrList<volVectorField> conv(surfFields.size());
            forAll(surfFields, i)
            {
                tmp<volVectorField> tc = fvc::div(surfFields[i],vecFields[0]); // uMean fixed, vary fMode_i
                conv.set(i, new volVectorField(tc));  // materialize tmp into heap object
                if (debug)
                    Info<< "div(" << surfFields[i].name() << ", " << vecFields[0].name() << ") stored [" << key << "]\n";
            }
            return conv;
        }

        //iterates over vectorfield
        PtrList<volVectorField> convectiveTermVec(const PtrList<volVectorField>& vecFields, const PtrList<surfaceScalarField>& surfFields, const word& key=word::null)
        {
            PtrList<volVectorField> conv(vecFields.size());
            forAll(vecFields, i)
            {
                tmp<volVectorField> tc = fvc::div(surfFields[0],vecFields[i]);
                conv.set(i, new volVectorField(tc));  // materialize tmp into heap object
                if (debug)
                    Info<< "div(" << surfFields[0].name() << ", " << vecFields[i].name() << ") stored [" << key << "]\n";
            }
            return conv;
        }

         //third order tensor with divergence operator 
        template<typename T1, typename TFlux>
        RectangularMatrix<scalar> thirdOrder(const PtrList<volVectorField>& proj,  const T1& field1, const TFlux& fluxField, word key = "")
        {
            label size = proj.size();
            RectangularMatrix<scalar> A(size, size*size, 0.0);
            //lambda function to get the index of the flattened matrix
            auto idx = [size](label n, label m) -> label { return n*size + m;};
            if constexpr (std::is_same<TFlux, PtrList<volVectorField>>::value)
            {
                for(label m = 0; m < size; ++m)
                { 
                    surfaceScalarField flux = linearInterpolate(fluxField[m]) & mesh_.Sf();
                    for(label n = 0; n < size; ++n)
                    {
                        const auto& div = fvc::div(flux, field1[n]);
                        for(label l = 0; l < size; ++l) 
                        {
                            A(l, idx(n,m)) = dotProductL2(proj[l], div);
                        }
                    }
                }
            }
            else if constexpr (std::is_same<TFlux, PtrList<surfaceScalarField>>::value)
            {
                for(label m = 0; m < size; ++m)
                { 
                    const auto& flux = fluxField[m];
                    for(label n = 0; n < size; ++n)
                    {
                        const auto& div = fvc::div(flux, field1[n]);
                        for(label l = 0; l < size; ++l) 
                        {
                            A(l, idx(n,m)) = dotProductL2(proj[l], div);
                        }
                    }
                }
            }
            debugMsg(debug, "non linear term computed");
            return A;
        }

        template<typename T1, typename TFlux>
        RectangularMatrix<scalar> thirdOrderPressProj(const PtrList<volScalarField>& proj,  const T1& field1, const TFlux& fluxField, word key = "")
        {
            label size = proj.size();
            RectangularMatrix<scalar> A(size, size*size, 0.0);
            //lambda function to get the index of the flattened matrix
            auto idx = [size](label n, label m) -> label { return n*size + m;};
            for(label m = 0; m < size; ++m)
            { 
                const auto& flux = fluxField[m];
                for(label n = 0; n < size; ++n)
                {
                    const auto& div = fvc::div(flux, field1[n]);
                    const auto& divDiv = fvc::div(div);
                    for(label l = 0; l < size; ++l) 
                    {
                        A(l, idx(n,m)) = dotProductL2Scalar(proj[l], divDiv);
                    }
                }
            }
            debugMsg(debug, "non linear term computed");
            return A;
        }

        template<typename T1, typename TFlux>
        RectangularMatrix<scalar> thirdOrderPressProj2(const PtrList<volScalarField>& proj,  
                                                    const T1& field1, 
                                                    const TFlux& fluxField, 
                                                    word key = "")
        {
            label size = proj.size();
            RectangularMatrix<scalar> A(size, size*size, 0.0);
            auto idx = [size](label n, label m) -> label { return n*size + m;};
            
            // Compute gradient of pressure modes ONCE, outside the loops
            PtrList<volVectorField> gradPModes = fieldGradient(proj, "gradient of pressure modes");
            
            for(label m = 0; m < size; ++m)
            { 
                const auto& flux = fluxField[m];
                for(label n = 0; n < size; ++n)
                {
                    const auto& convTerm = fvc::div(flux, field1[n]);  // This is ∇·(φ_n⊗φ_m)
                    
                    for(label l = 0; l < size; ++l) 
                    {
                        // Inner product: <∇χ_l, ∇·(φ_n⊗φ_m)>
                        A(l, idx(n,m)) = dotProductL2(gradPModes[l], convTerm);
                    }
                }
            }
            debugMsg(debug, "PPE non linear term computed");
            return A;
        }

        template<typename T>
        RectangularMatrix<scalar> firstOrder(const PtrList<T>& proj, const PtrList<T>& field, word key = "")
        {
            label size = proj.size();
            RectangularMatrix<scalar> A(size, 1, 0.0);
            if constexpr (std::is_same<T, volVectorField>::value)
            {    
                for(label l = 0; l < size; ++l)
                {
                    A(l,0) = dotProductL2(proj[l],field[0]);
                }
            }
            else if constexpr (std::is_same<T, volScalarField>::value)
            {    
                for(label l = 0; l < size; ++l)
                {
                    A(l,0) = dotProductL2Scalar(proj[l], field[0]);
                }
            }
            debugMsg(debug, key, " computed");
            return A;
        }

        template<typename T>
        RectangularMatrix<scalar> secondOrder(const PtrList<T>& proj, const PtrList<T>& field, word key = "")
        {
            label size = proj.size();
            label size2 = field.size();
            RectangularMatrix<scalar> A(size, size2, 0.0);
            if constexpr (std::is_same<T, volVectorField>::value)
            {    
                for(label m = 0; m < size2; ++m)
                {
                    const auto& rhs = field[m];
                    for(label l = 0; l < size; ++l)
                    {
                        A(l,m) = dotProductL2(proj[l], rhs);
                    }
                }
            }
            else if constexpr (std::is_same<T, volScalarField>::value)
            {    
                for(label m = 0; m < size2; ++m)
                {
                    const auto& rhs = field[m];
                    for(label l = 0; l < size; ++l)
                    {
                        A(l,m) = dotProductL2Scalar(proj[l], rhs);
                    }
                }
            }
            debugMsg(debug, key, " computed");
            return A;
        }

        scalar dotProductL2(const volVectorField& field1, const volVectorField& field2, bool debug = false, word key = "")
        {
            auto contraction = field1 & field2;
            //debugMsg(debug, "dotProd of ", key, " computed successfully");
            return fvc::domainIntegrate(contraction).value();
        };

        scalar dotProductL2Scalar(const volScalarField& field1, const volScalarField& field2, bool debug = false, word key = "")
        {
            auto contraction = field1 * field2;
            //debugMsg(debug, "dotProd of ", key, " computed successfully");
            return fvc::domainIntegrate(contraction).value();
        };
}; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
