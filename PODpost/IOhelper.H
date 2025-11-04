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

#ifndef IOHELPER_H
#define IOHELPER_H

#include "fvCFD.H"
#include <Eigen/Dense>
#include "Utilities.H"
#include "calculatedFvsPatchFields.H"
#include "fixedValueFvsPatchFields.H"

namespace Foam
{
class IO: public utilities
{
    public:
        IO() = default;                         // Constructor
        IO(const IO&) = delete;          // no copy constructor
        void operator=(const IO&) = delete;     // no assignment
        virtual ~IO() = default;                // default destructor

        /****************************************/
        /*********** MEMBER FUNCTIONS ***********/
        /****************************************/
        
        //openfoam DynamicList to eigen matrix
        template<typename T>
        Eigen::MatrixXd foamListToEigenMatrix(const DynamicList<T>& foamList, bool debug = false, word key = "")
        {
            //list stores with rows being the snapshots and columns being the dofs --> in eigen we want the opposite
            label Ns = foamList.size();
            label nDof = foamList[0].size();
            label nRowsEigen = nDof;
            label nComponents = 1;
            if constexpr (std::is_same<T,vectorField>::value)
            {
                nComponents = foamList[0][0].size();
                nRowsEigen *= nComponents;
            } 

            Eigen::MatrixXd eigenMatrix(nRowsEigen, Ns);
            for(label i = 0; i < Ns; ++i) //rows of foam matrix
            {
                for(label j = 0; j < nDof; ++j) //cols of foam matrix
                {
                    if constexpr (std::is_same<T, vectorField>::value)
                    {
                        for(label k = 0; k < nComponents; ++k)
                        {
                            eigenMatrix(j * nComponents + k, i) = foamList[i][j][k];
                        } 
                    }
                    else if constexpr (std::is_same<T, scalarField>::value || std::is_same<T,vector>::value)
                    {
                        eigenMatrix(j,i) = foamList[i][j];
                    }
                }
            }
            debugMsg("eigen Matrix *" + key + "* has dimensions: (" + std::to_string(eigenMatrix.rows()) + ", " + std::to_string(eigenMatrix.cols()) + ")", debug);
            return eigenMatrix;
        } 

        //openfoam DynamicList to eigen vector  -->works for BCs
        template<typename T>
        Eigen::VectorXd foamListToEigenVector(const DynamicList<T>& foamList, bool debug = false, word key = "")
        {
            label cols = foamList.size();
            Eigen::VectorXd eigenVector(cols);
            for(label i = 0; i < cols; ++i)
            {
                eigenVector(i) = foamList[i];
            }
            debugMsg("eigen Vector *" + key + "* has dimensions: " + std::to_string(eigenVector.rows()), debug);
            return eigenVector;
        }

        //works for fields inside list e.g. for the mean
        template<typename T>
        Eigen::VectorXd foamListFieldToEigenVector(const DynamicList<T>& foamList, bool debug = false, word key = "")
        {
            label cols = foamList[0].size();
            Eigen::VectorXd eigenVector(cols);
            for(label i = 0; i < cols; ++i)
            {
                eigenVector(i) = foamList[0][i];
            }
            debugMsg("eigen Vector *" + key + "* has dimensions: " + std::to_string(eigenVector.rows()), debug);
            return eigenVector;
        }


        //function to create the openfoam folder environment with 0 / constant / system --> this ensures that paraView can open the datas in the folder
        void createFolderEnvironment(const fileName& newPath, fvMesh& mesh, Time& runTime, bool debug = false)
        {
            // Reset target folder (remove if exists, then create)
            if (isDir(newPath))
            {
                rmDir(newPath);    // remove entire tree
            }
            mkDir(newPath);        // (re)create

            const fileName basePath = runTime.path();

            // the three standard OpenFOAM dirs we want to mirror
            static const wordList envDirs{ "0", "constant", "system"};

            for (const word& d : envDirs)
            {
                const fileName src = basePath / d;
                const fileName dst = newPath / d;

                if (isDir(src))
                {
                    // copy only if destination doesn't exist yet
                    if (!isDir(dst)) cp(src, dst);
                }
                else
                {
                    WarningInFunction
                        << "Source directory does not exist: " << src << nl << endl;
                }
            }
            debugMsg("folder environment *" + newPath + "* was created succesfully!", debug);
        }


        //helper function for flux BCs
        template<class T>
        inline void changeCalculatedBCsforSurfaceScalarFields(T& fld,  const typename T::value_type& BCvalue, bool debug=false)
        {

            auto& bfr = fld.boundaryFieldRef();
            forAll(bfr, patchI)
            {
                fvsPatchScalarField& pf = bfr[patchI];

                if (isA<calculatedFvsPatchScalarField>(pf))
                {
                    // replace calculated -> fixedValue
                    bfr.set
                    (
                        patchI,
                        new  fixedValueFvsPatchScalarField( pf.patch(), fld )
                    );
                    refCast<fixedValueFvsPatchScalarField>(bfr[patchI]) == BCvalue;
                    if (debug) Info<< "Patch " << bfr[patchI].patch().name()
                                << ": calculated -> fixedValue(0)\n";
                }
                else if (isA<fixedValueFvsPatchScalarField>(pf))
                {
                    // already fixedValue: make it homogeneous (zero)
                    refCast<fixedValueFvsPatchScalarField>(pf) == BCvalue;
                    if (debug) Info<< "Patch " << pf.patch().name()
                                << ": fixedValue set to 0\n";
                }
                // else: leave other BC types unchanged
            }
        }

        //helper function to take an eigenMatrix and save it to disk in foam format --> choose a relative path, can only save on that single path
        template<typename T>
        void EigenVecWriteToDisk(const Eigen::VectorXd& vec, const fvMesh& mesh, Time& runTime, word fieldName, word newFieldName, fileName relPath, bool debug = false, word key = "")
        {

            instantList timeFolders = getTimeDirs(runTime, true, debug);

            T initializer 
            (
                IOobject
                (
                    fieldName,
                    timeFolders[0].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            T newField
            (
                IOobject
                (
                    newFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                initializer
                //mesh, 
                //initializer.dimensions()     
                //initializer.boundaryField().types()
                //dimensioned<typename T::value_type>("zero", initializer.dimensions(), Zero)
            );

            auto& internal = newField.internalFieldRef();
            label internalCells = internal.size();
            const label nElems = internal.size();     // cells for vol*, internal faces for surface*

            // --- shape checks
            if constexpr (std::is_same<T, volScalarField>::value)
            {
                if (vec.rows() != nElems)
                {
                    FatalErrorInFunction
                        << "volScalarField size mismatch: vec.rows()=" << vec.rows()
                        << " vs nCells=" << nElems << nl << abort(FatalError);
                }
                for (label i=0; i<nElems; ++i) internal[i] = vec(i);
            }
            else if constexpr (std::is_same<T, volVectorField>::value)
            {
                if (vec.rows() % nElems != 0)
                {
                    FatalErrorInFunction
                        << "volVectorField size mismatch: vec.rows()=" << vec.rows()
                        << " not divisible by nCells=" << nElems << nl << abort(FatalError);
                }
                const label nComponents = vec.rows() / nElems; // expect 3
                if (nComponents != 3)
                {
                    WarningInFunction
                        << "volVectorField with nComponents=" << nComponents
                        << " (expected 3); proceeding anyway." << nl;
                }
                for (label i=0; i<nElems; ++i)
                    for (label k=0; k<nComponents; ++k)
                        internal[i][k] = vec(i*nComponents + k);
            }
            else if constexpr (std::is_same<T, surfaceScalarField>::value)
            {
                // internal size is nInternalFaces
                if (vec.rows() != nElems)
                {
                    FatalErrorInFunction
                        << "surfaceScalarField size mismatch: vec.rows()=" << vec.rows()
                        << " vs nInternalFaces=" << nElems << nl << abort(FatalError);
                }
                for (label i=0; i<nElems; ++i) internal[i] = vec(i);
            }


            const label oldIndex = runTime.timeIndex();
            const instant oldInst(runTime.time().value(), runTime.timeName());
            //trick runtime to go in relpath
            runTime.setTime(instant(runTime.time().value(), relPath), oldIndex);
            newField.write();
            debugMsg("field *" + key + "* successfully wrote to folder: " + runTime.path().name() + "/" + runTime.timeName(), debug);
            // Restore runtime I had
            runTime.setTime(oldInst, oldIndex);
        }

        //writes to disk but also changes BC to wanted value
        template<typename T>
        void EigenVecWriteToDiskBCchange(const Eigen::VectorXd& vec, typename T::value_type BCvalue, const fvMesh& mesh, Time& runTime, word fieldName, word newFieldName, fileName relPath, bool debug = false, word key = "")
        {

            instantList timeFolders = getTimeDirs(runTime, true, debug);

            T initializer 
            (
                IOobject
                (
                    fieldName,
                    timeFolders[0].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            T newField
            (
                IOobject
                (
                    newFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh, 
                initializer.dimensions(),           
                initializer.boundaryField().types()
            );

            auto& internal = newField.internalFieldRef();
            label internalCells = internal.size();
            label nRows = vec.rows();
            if constexpr (std::is_same<T,volVectorField>::value)
            {
                label nComponents = nRows / internalCells;
                for(int i = 0; i < internalCells; ++i)
                {
                    for(int k = 0; k < nComponents; ++k)
                    {
                        internal[i][k] = vec(i * nComponents + k);
                    }
                }
            }
            else if constexpr (std::is_same<T,volScalarField>::value || std::is_same<T,surfaceScalarField>::value) 
            {
                for(int i = 0; i < internalCells; ++i)
                {
                    internal[i] = vec(i);
                }
            }

            //change BC to wanted value
            const fvBoundaryMesh& boundary = mesh.boundary();
            forAll(boundary,patchI)
            {
                const word patchName = boundary[patchI].name();
                auto& pf = newField.boundaryFieldRef()[patchI];
                if(isA<fixedValueFvPatchField<typename T::value_type>>(pf))
                {
                    refCast<fixedValueFvPatchField<typename T::value_type>>(pf) == BCvalue;
                }
            }

            if constexpr (std::is_same<T,surfaceScalarField>::value)
                {
                    changeCalculatedBCsforSurfaceScalarFields(newField, BCvalue, debug);
                }

            newField.correctBoundaryConditions();

            const label oldIndex = runTime.timeIndex();
            const instant oldInst(runTime.time().value(), runTime.timeName());
            //trick runtime to go in relpath
            runTime.setTime(instant(runTime.time().value(), relPath), oldIndex);
            newField.write();
            debugMsg("field *" + key + "* successfully wrote to folder: " + runTime.path().name() + "/" + runTime.timeName(), debug);
            // Restore runtime I had
            runTime.setTime(oldInst, oldIndex);
        }

        //writes matrix to disk to relative path/time --> start at first time
        template<typename T> 
        void EigenMatrixWriteToDisk(const Eigen::MatrixXd& A, const fvMesh& mesh, Time& runTime, word fieldName, word newFieldName, fileName relPath, scalar firstTime = 0, scalar separation = 0, bool debug = false, word key = "")
        {
            //initializing the foam field 
            instantList timeFolders = getTimeDirs(runTime, true, debug);

            T initializer 
            (
                IOobject
                (
                    fieldName,
                    timeFolders[0].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            T newField
            (
                IOobject
                (
                    newFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh, 
                initializer.dimensions(),           
                initializer.boundaryField().types()
                //dimensioned<typename T::value_type>("zero", initializer.dimensions(), Zero)
            );

            auto& internal = newField.internalFieldRef();
            label internalCells = internal.size();
            label nRows = A.rows();
            //creating the saving place
            scalar lastTime = A.cols();
            instantList customSaving = createCustomTime(runTime, firstTime, lastTime, separation, debug);
            label nComponents = nRows / internalCells;

            const label oldIndex = runTime.timeIndex();
            const instant oldInst(runTime.time().value(), runTime.timeName());

            //loop over columns of A
            forAll(customSaving,ti)
            {
                runTime.setTime(customSaving[ti],ti);
                //auto& internal = newField.internalFieldRef();
                const scalar t = customSaving[ti].value();
                const label colI = lround((t - firstTime)/separation);
                if constexpr (std::is_same<T,volVectorField>::value)
                {
                    for(int i = 0; i < internalCells; ++i)
                    {
                        for(int k = 0; k < nComponents; ++k)
                        {
                            internal[i][k] = A(i * nComponents + k, colI);
                        }
                    }
                }
                else if constexpr (std::is_same<T,volScalarField>::value || std::is_same<T,surfaceScalarField>::value) 
                {
                    for(int i = 0; i < internalCells; ++i)
                    {
                        internal[i] = A(i, colI);
                    }
                }
                scalar saveTime = ti + firstTime;
                runTime.setTime(instant(runTime.time().value(), relPath/name(saveTime)), oldIndex);
                newField.write();
            }
            runTime.setTime(oldInst, oldIndex);
        }


        //writes matrix to disk to relative path/time --> start at first time, and changes BCs fixed to what I want (probably could do 1 single function for vec and matrix but not for now..)
        template<typename T> 
        void EigenMatrixWriteToDiskBCchange(const Eigen::MatrixXd& A, typename T::value_type BCvalue, const fvMesh& mesh, Time& runTime, word fieldName, word newFieldName, fileName relPath, scalar firstTime = 0, scalar separation = 0, bool debug = false, word key = "")
        {
            //initializing the foam field 
            instantList timeFolders = getTimeDirs(runTime, true, debug);

            T initializer 
            (
                IOobject
                (
                    fieldName,
                    timeFolders[0].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            T newField
            (
                IOobject
                (
                    newFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh, 
                initializer.dimensions(),           
                initializer.boundaryField().types()
                //dimensioned<typename T::value_type>("zero", initializer.dimensions(), Zero)
            );

            auto& internal = newField.internalFieldRef();
            label internalCells = internal.size();
            label nRows = A.rows();
            //creating the saving place
            scalar lastTime = A.cols();
            instantList customSaving = createCustomTime(runTime, firstTime, lastTime, separation, debug);
            label nComponents = nRows / internalCells;

            const label oldIndex = runTime.timeIndex();
            const instant oldInst(runTime.time().value(), runTime.timeName());

            //loop over columns of A
            forAll(customSaving,ti)
            {
                runTime.setTime(customSaving[ti],ti);
                //auto& internal = newField.internalFieldRef();
                const scalar t = customSaving[ti].value();
                const label colI = lround((t - firstTime)/separation);
                if constexpr (std::is_same<T,volVectorField>::value)
                {
                    for(int i = 0; i < internalCells; ++i)
                    {
                        for(int k = 0; k < nComponents; ++k)
                        {
                            internal[i][k] = A(i * nComponents + k, colI);
                        }
                    }
                }
                else if constexpr (std::is_same<T,volScalarField>::value || std::is_same<T,surfaceScalarField>::value) 
                {
                    for(int i = 0; i < internalCells; ++i)
                    {
                        internal[i] = A(i, colI);
                    }
                }

                //change BC to wanted value
                const fvBoundaryMesh& boundary = mesh.boundary();
                forAll(boundary,patchI)
                {
                    const word patchName = boundary[patchI].name();
                    auto& pf = newField.boundaryFieldRef()[patchI];
                    if(isA<fixedValueFvPatchField<typename T::value_type>>(pf))
                    {
                        refCast<fixedValueFvPatchField<typename T::value_type>>(pf) == BCvalue;
                    }
                }

                if constexpr (std::is_same<T,surfaceScalarField>::value)
                {
                    changeCalculatedBCsforSurfaceScalarFields(newField, BCvalue, debug);
                }

                newField.correctBoundaryConditions();

                scalar saveTime = ti + firstTime;
                runTime.setTime(instant(runTime.time().value(), relPath/name(saveTime)), oldIndex);
                newField.write();
            }
            runTime.setTime(oldInst, oldIndex);
        }
    
         template<typename T> 
        void EigenMatrixWriteToDiskRunTime(const Eigen::MatrixXd& A, const fvMesh& mesh, Time& runTime, scalar startTime, word fieldName, word newFieldName, bool debug = false, word key = "")
        {
            T initializer 
            (
                IOobject
                (
                    fieldName,
                    runTime.path()/"0",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            T newField
            (
                IOobject
                (
                    newFieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                initializer
            );

            auto& internal = newField.internalFieldRef();
            label internalCells = internal.size();
            label nRows = A.rows();
            label nComponents = nRows / internalCells;
            //std::cout << "nComponents: " << nComponents << "\n nRows: " << nRows << "\n internalCells " << internalCells << std::endl; 
            //std::cout << "cols of A are: " <<  A.cols() << std::endl;
            //loop over columns of A
            runTime.setTime(startTime,0);
            label colI = 0;
            while(runTime.loop())
            {
                if(runTime.outputTime())
                {
                    if (colI >= A.cols())
                    {
                        if (debug) Info<< "Writer: no more columns in A ("
                                    << colI << " >= " << A.cols() << "), stopping.\n";
                        break;
                    }
                    if constexpr (std::is_same<T,volVectorField>::value)
                    {
                        for(int i = 0; i < internalCells; ++i)
                        {
                            for(int k = 0; k < nComponents; ++k)
                            {
                                internal[i][k] = A(i * nComponents + k, colI);
                            }
                        }
                    }
                    else if constexpr (std::is_same<T,volScalarField>::value || std::is_same<T,surfaceScalarField>::value) 
                    {
                        for(int i = 0; i < internalCells; ++i)
                        {
                            internal[i] = A(i, colI);
                        }
                    }
                    ++colI;
                    newField.write();
                }            
            }
        }
    
}; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
