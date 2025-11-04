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

        //openfoam DynamicList to eigen vector
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
                mesh, 
                initializer.dimensions(),           
                initializer.boundaryField().types()
                //dimensioned<typename T::value_type>("zero", initializer.dimensions(), Zero)
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

        //reads field from disk and stores it in a ptrlist --> not stored in time folders, but custom folders e.g. for modes
        template<typename T>
        PtrList<T> readPtrField_noTime(word fieldName, const Time& runTime, const fvMesh& mesh, scalar start, scalar end, fileName relPath, bool debug = false, word key = "")
        {
            fileName path = runTime.path()/relPath;
            PtrList<T> fieldPtr;
            for(int i = start; i <= end; ++i)
            {
                T* ptr = new T
                {
                    IOobject
                    (
                        fieldName,
                        path/name(i),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh
                };
                fieldPtr.append(ptr);
            }
            debugMsg(debug, "fieldPtr ", key, " read successfully from ", path, " and has size: ", fieldPtr.size(), ", " , fieldPtr[0].size(), ")");
            return fieldPtr;
        }

        //reads field from disk and stores it in a ptrlist --> stored in timefolders (runtime)       
        template<typename T>
        PtrList<T> readPtrField(word fieldName, Time& runTime, const fvMesh& mesh, bool excludeZero = false, bool debug = false, word key = "")
        {
            PtrList<T> fieldPtr;
            instantList timeFolderList = getTimeDirs(runTime,excludeZero);

            forAll(timeFolderList, ti)
            {
                runTime.setTime(timeFolderList[ti],ti);
                T* ptr = new T
                {
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::READ_IF_PRESENT,
                        IOobject::NO_WRITE
                    ),
                    mesh
                };
                fieldPtr.append(ptr);
            }
            debugMsg(debug, "fieldPtr ", key, " read successfully from timeFolders and has size: (", fieldPtr.size(), ", " , fieldPtr[0].size(), ")");
            return fieldPtr;
        }
    
        template<typename T>
        void saveFoamMatrixToDisk(RectangularMatrix<T> M, fileName path, bool debug = false, word key = "")
        {
            OFstream os(path);
            os << M;
            debugMsg(debug, key, " successfully saved to ", path);
        }


        template<typename T>
        RectangularMatrix<T> readFoamMatrixFromDisk(fileName path, bool debug = false)
        {
            IFstream is(path);
            if (!is.good())
            {
                FatalErrorInFunction << "Cannot open " << path << exit(FatalError);
            }
            RectangularMatrix<T> A;
            is >> A;
            debugMsg(debug,"field in *", path, "* read successfully");
            return A;
        }
}; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
