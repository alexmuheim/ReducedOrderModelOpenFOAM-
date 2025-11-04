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

#ifndef MATHROM2_H
#define MATHROM2_H

#include "fvCFD.H"
#include <Eigen/Dense>
#include "Utilities.H"

namespace Foam
{
class mathOffline2: public utilities
{
    private: 
        const fvMesh& mesh_;
        bool debug;
    
    public:
        explicit mathOffline2(const fvMesh& mesh, bool debug) 
        : mesh_(mesh), debug(debug)
        {}
        
        mathOffline2(const mathOffline2&) = delete;
        void operator=(const mathOffline2&) = delete;
        virtual ~mathOffline2() = default;

        /****************************************/
        /*********** MEMBER FUNCTIONS ***********/
        /****************************************/

        // Keep these helper functions for when you need the fields explicitly
        PtrList<volVectorField> fieldGradient(const PtrList<volScalarField>& fields, const word& key=word::null)
        {
            PtrList<volVectorField> grads(fields.size());
            forAll(fields, i)
            {
                tmp<volVectorField> tg = fvc::grad(fields[i]);
                grads.set(i, new volVectorField(tg));
                if (debug)
                    Info<< "grad(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return grads;
        }

       /*  template<typename T>
        PtrList<T> fieldLaplacian(const PtrList<T>& fields, const word& key=word::null)
        {
            PtrList<T> lapl(fields.size());
            forAll(fields, i)
            {
                tmp<T> tl = fvc::laplacian(fields[i]);
                lapl.set(i, new T(tl));
                if (debug)
                    Info<< "laplacian(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return lapl;
        } */

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
                divs.set(i, new T1(td));
                if (debug)
                    Info<< "div(" << fields[i].name() << ") stored [" << key << "]\n";
            }
            return divs;
        }
        
        PtrList<volVectorField> convectiveTermFlux(const PtrList<volVectorField>& vecFields, 
                                                    const PtrList<surfaceScalarField>& surfFields, 
                                                    const word& key=word::null)
        {
            PtrList<volVectorField> conv(surfFields.size());
            forAll(surfFields, i)
            {
                tmp<volVectorField> tc = fvc::div(surfFields[i], vecFields[0]);
                conv.set(i, new volVectorField(tc));
                if (debug)
                    Info<< "div(" << surfFields[i].name() << ", " << vecFields[0].name() 
                        << ") stored [" << key << "]\n";
            }
            return conv;
        }

        PtrList<volVectorField> convectiveTermVec(const PtrList<volVectorField>& vecFields, 
                                                   const PtrList<surfaceScalarField>& surfFields, 
                                                   const word& key=word::null)
        {
            PtrList<volVectorField> conv(vecFields.size());
            forAll(vecFields, i)
            {
                tmp<volVectorField> tc = fvc::div(surfFields[0], vecFields[i]);
                conv.set(i, new volVectorField(tc));
                if (debug)
                    Info<< "div(" << surfFields[0].name() << ", " << vecFields[i].name() 
                        << ") stored [" << key << "]\n";
            }
            return conv;
        }

        /****************************************/
        /******** ITHACA-FV STYLE METHODS *******/
        /****************************************/

        // Mass matrix: M[i][j] = ∫ φ_i · φ_j dV
        RectangularMatrix<scalar> massMatrix(const PtrList<volVectorField>& modes, 
                                              const word& key = "")
        {
            label size = modes.size();
            RectangularMatrix<scalar> M(size, size, 0.0);
            
            for(label i = 0; i < size; i++)
            {
                for(label j = 0; j < size; j++)
                {
                    M(i, j) = fvc::domainIntegrate(modes[i] & modes[j]).value();
                }
            }
            
            debugMsg(debug, "Mass matrix ", key, " computed");
            return M;
        }

        // Stiffness/Diffusion matrix: B[i][j] = ∫ φ_i · ∇²φ_j dV
        RectangularMatrix<scalar> diffusiveMatrix(const PtrList<volVectorField>& modes, 
                                                   const word& key = "")
        {
            label size = modes.size();
            RectangularMatrix<scalar> B(size, size, 0.0);
            
            for(label i = 0; i < size; i++)
            {
                for(label j = 0; j < size; j++)
                {
                    B(i, j) = fvc::domainIntegrate(modes[i] & fvc::laplacian(modes[j])).value();
                }
            }
            
            debugMsg(debug, "Diffusive matrix ", key, " computed");
            return B;
        }

        // Pressure gradient matrix: K[i][j] = ∫ φ_u[i] · ∇φ_p[j] dV
        RectangularMatrix<scalar> pressureGradientMatrix(const PtrList<volVectorField>& uModes,
                                                          const PtrList<volScalarField>& pModes,
                                                          const word& key = "")
        {
            label Nu = uModes.size();
            label Np = pModes.size();
            RectangularMatrix<scalar> K(Nu, Np, 0.0);
            
            for(label i = 0; i < Nu; i++)
            {
                for(label j = 0; j < Np; j++)
                {
                    K(i, j) = fvc::domainIntegrate(uModes[i] & fvc::grad(pModes[j])).value();
                }
            }
            
            debugMsg(debug, "Pressure gradient matrix ", key, " computed");
            return K;
        }

        // Divergence matrix: D[j][i] = ∫ φ_p[j] * (∇·φ_u[i]) dV
        RectangularMatrix<scalar> divergenceMatrix(const PtrList<volScalarField>& pModes,
                                                    const PtrList<volVectorField>& uModes,
                                                    const word& key = "")
        {
            label Np = pModes.size();
            label Nu = uModes.size();
            RectangularMatrix<scalar> D(Np, Nu, 0.0);
            
            for(label j = 0; j < Np; j++)
            {
                for(label i = 0; i < Nu; i++)
                {
                    D(j, i) = fvc::domainIntegrate(pModes[j] * fvc::div(uModes[i])).value();
                }
            }
            
            debugMsg(debug, "Divergence matrix ", key, " computed");
            return D;
        }

        // Pressure Laplacian: P[i][j] = ∫ ∇φ_p[i] · ∇φ_p[j] dV
        RectangularMatrix<scalar> pressureLaplacian(const PtrList<volScalarField>& pModes,
                                                     const word& key = "")
        {
            label Np = pModes.size();
            RectangularMatrix<scalar> P(Np, Np, 0.0);
            
            for(label i = 0; i < Np; i++)
            {
                for(label j = 0; j < Np; j++)
                {
                    P(i, j) = fvc::domainIntegrate(fvc::grad(pModes[i]) & 
                                                    fvc::grad(pModes[j])).value();
                }
            }
            
            debugMsg(debug, "Pressure Laplacian ", key, " computed");
            return P;
        }

        // Convection tensor: C[i][j][k] = ∫ φ_i · ∇·(φ_j ⊗ φ_k) dV
        template<typename TFlux>
        RectangularMatrix<scalar> convectionTensor(const PtrList<volVectorField>& proj,  
                                                    const PtrList<volVectorField>& field, 
                                                    const TFlux& fluxField, 
                                                    const word& key = "")
        {
            label size = proj.size();
            RectangularMatrix<scalar> C(size, size*size, 0.0);
            auto idx = [size](label n, label m) -> label { return n*size + m; };
            
            if constexpr (std::is_same<TFlux, PtrList<volVectorField>>::value)
            {
                for(label m = 0; m < size; m++)
                { 
                    surfaceScalarField flux = fvc::interpolate(fluxField[m]) & mesh_.Sf();
                    for(label n = 0; n < size; n++)
                    {
                        for(label l = 0; l < size; l++) 
                        {
                            C(l, idx(n,m)) = fvc::domainIntegrate(proj[l] & 
                                            fvc::div(flux, field[n])).value();
                        }
                    }
                }
            }
            else if constexpr (std::is_same<TFlux, PtrList<surfaceScalarField>>::value)
            {
                for(label m = 0; m < size; m++)
                { 
                    for(label n = 0; n < size; n++)
                    {
                        for(label l = 0; l < size; l++) 
                        {
                            C(l, idx(n,m)) = fvc::domainIntegrate(proj[l] & fvc::div(fluxField[m], field[n])).value();
                        }
                    }
                }
            }
            
            debugMsg(debug, "Convection tensor ", key, " computed");
            return C;
        }

        // PPE convection tensor: F[l][n,m] = ∫ ∇χ_l · ∇·(φ_n ⊗ φ_m) dV
        template<typename TFlux>
        RectangularMatrix<scalar> ppeConvectionTensor(const PtrList<volScalarField>& pModes,  
                                                       const PtrList<volVectorField>& uModes, 
                                                       const TFlux& fluxField, 
                                                       const word& key = "")
        {
            label size = pModes.size();
            RectangularMatrix<scalar> F(size, size*size, 0.0);
            auto idx = [size](label n, label m) -> label { return n*size + m; };
            
            for(label m = 0; m < size; m++)
            { 
                for(label n = 0; n < size; n++)
                {
                    for(label l = 0; l < size; l++) 
                    {
                        // ∫ ∇χ_l · ∇·(φ_n ⊗ φ_m) dV
                        F(l, idx(n,m)) = fvc::domainIntegrate(fvc::grad(pModes[l]) & 
                                        fvc::div(fluxField[m], uModes[n])).value();
                    }
                }
            }
            
            debugMsg(debug, "PPE convection tensor ", key, " computed");
            return F;
        }

        // Generic second order with operator: A[i][j] = ∫ φ_i · Op(φ_j) dV
        template<typename T>
        RectangularMatrix<scalar> secondOrder(const PtrList<T>& proj, 
                                               const PtrList<T>& field, 
                                               const word& key = "")
        {
            label size1 = proj.size();
            label size2 = field.size();
            RectangularMatrix<scalar> A(size1, size2, 0.0);
            
            if constexpr (std::is_same<T, volVectorField>::value)
            {    
                for(label i = 0; i < size1; i++)
                {
                    for(label j = 0; j < size2; j++)
                    {
                        A(i, j) = fvc::domainIntegrate(proj[i] & field[j]).value();
                    }
                }
            }
            else if constexpr (std::is_same<T, volScalarField>::value)
            {    
                for(label i = 0; i < size1; i++)
                {
                    for(label j = 0; j < size2; j++)
                    {
                        A(i, j) = fvc::domainIntegrate(proj[i] * field[j]).value();
                    }
                }
            }
            
            debugMsg(debug, key, " computed");
            return A;
        }

        // First order: A[i] = ∫ φ_i · f dV
        template<typename T>
        RectangularMatrix<scalar> firstOrder(const PtrList<T>& proj, 
                                              const PtrList<T>& field, 
                                              const word& key = "")
        {
            label size = proj.size();
            RectangularMatrix<scalar> A(size, 1, 0.0);
            
            if constexpr (std::is_same<T, volVectorField>::value)
            {    
                for(label i = 0; i < size; i++)
                {
                    A(i, 0) = fvc::domainIntegrate(proj[i] & field[0]).value();
                }
            }
            else if constexpr (std::is_same<T, volScalarField>::value)
            {    
                for(label i = 0; i < size; i++)
                {
                    A(i, 0) = fvc::domainIntegrate(proj[i] * field[0]).value();
                }
            }
            
            debugMsg(debug, key, " computed");
            return A;
        }

        // Keep old functions for backward compatibility if needed
        scalar dotProductL2(const volVectorField& field1, 
                           const volVectorField& field2, 
                           bool debug = false, 
                           const word& key = "")
        {
            return fvc::domainIntegrate(field1 & field2).value();
        }

        scalar dotProductL2Scalar(const volScalarField& field1, 
                                  const volScalarField& field2, 
                                  bool debug = false, 
                                  const word& key = "")
        {
            return fvc::domainIntegrate(field1 * field2).value();
        }
}; 
}

#endif
