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
class math: public utilities
{
    private: const fvMesh& mesh_; //reference to mesh
    public:
        explicit math(const fvMesh& mesh) 
        : mesh_(mesh)
        {}                         // Constructor
        math(const math&) = delete;          // no copy constructor
        void operator=(const math&) = delete;     // no assignment
        virtual ~math() = default;                // default destructor

        /****************************************/
        /*********** MEMBER FUNCTIONS ***********/
        /****************************************/
        
        //perform arithmetic mean of each column of a eigen matrix
        Eigen::VectorXd arithmeticMeanColWise(const Eigen::MatrixXd& eigenMatrix, bool debug = false, word key = "")        
        {
            Eigen::VectorXd mean = eigenMatrix.rowwise().mean();
            debugMsg("arithmetic mean of **" + key + "** computed without scaling --> dimension is: (" + std::to_string(mean.rows()) + ", " + std::to_string(mean.cols()) + ")", debug);
            return mean;
        }

        Eigen::VectorXd arithmeticMeanColWise(const Eigen::MatrixXd& eigenMatrix, const Eigen::VectorXd& scalingValues, bool debug = false, word key = "")        
        {
            Eigen::VectorXd mean = eigenMatrix.rowwise().mean();
            double meanScaling = scalingValues.mean();
            double scalingFactor = 1.0 / stabilise(meanScaling, VSMALL);
            debugMsg(debug,"scaling factor computed");
            debugMsg("arithmetic mean of *" + key + "* computed with scaling :" + std::to_string(scalingFactor) + "--> dimension is: (" + std::to_string(mean.rows()) + ", " + std::to_string(mean.cols()) + ")", debug);
            return mean * scalingFactor;
        }

        Eigen::MatrixXd computeFluct(const Eigen::MatrixXd& S, const Eigen::VectorXd& Smean, const Eigen::VectorXd& scaleBc, bool debug = false, word key = "")
        {
            assert(S.cols() == scaleBc.size());
            assert(S.rows() == Smean.size());
            Eigen::MatrixXd Sfluct = S;
            for(int i = 0; i < S.cols(); ++i)
            {
                Sfluct.col(i) -= Smean * scaleBc(i);
            }
            debugMsg("fluctuation field of *" + key + "* successfully created", debug);
            return Sfluct;
        }

        Eigen::MatrixXd computeFluct(const Eigen::MatrixXd& S, const Eigen::VectorXd& Smean, bool debug = false, word key = "")
        {
            assert(S.rows() == Smean.size());
            Eigen::MatrixXd Sfluct = S;
            for(int i = 0; i < S.cols(); ++i)
            {
                Sfluct.col(i) -= Smean;
            }
            debugMsg("fluctuation field of *" + key + "* successfully created", debug);
            return Sfluct;
        }
    
        Eigen::MatrixXd correlationMatrix(const Eigen::MatrixXd& S, bool debug =false, word key = "")
        {
            label nCells = mesh_.nCells(); //number of cell centers
            label Ns = S.cols();     //number of snapshots
            label nDof = S.rows();   //number of rows of snapshots matrix
            Eigen::MatrixXd C(Ns,Ns);
            label nComponents = nDof/nCells;
            debugMsg("This is the number of components of **" + key + "** for the correlation matrix " + std::to_string(nComponents), debug);
            
            Eigen::VectorXd intermediate(nCells);
            for(label i = 0; i < Ns; ++i)
            {
                for(label j = 0; j < Ns; ++j)
                {
                    for(label k = 0; k < nCells; ++k)
                    {
                        scalar value = 0.0;
                        for(label l = 0; l < nComponents; ++l)
                        {
                            label index = k*nComponents + l;
                            value += S(index, i) * S(index, j);
                        }
                        intermediate(k) = value;
                    }
                    scalar integral = 0.0;
                    for(label m = 0; m < nCells; ++m)
                    {
                        integral += intermediate(m) * mesh_.V()[m];
                    }
                    //C(i,j) = integral / Ns;
                    C(i,j) = integral;
                }
            }

            if(debug)
            {
                double relAsym = (C - C.transpose()).norm() / std::max(1e-16, C.norm());
                if (relAsym > 1e-10) std::cout<< "Warning: noticeable asymmetry: " << relAsym << nl;
            }
            //ensures that C is symmetric also numerically (just numerical artifact to make sure everything works)
            C = C.selfadjointView<Eigen::Lower>();
            if(debug)
            {
                double relAsym = (C - C.transpose()).norm() / std::max(1e-16, C.norm());
                if (relAsym > 1e-10) std::cout<< "Warning: noticeable asymmetry still in place after correction: " << relAsym << nl;
            }
            debugMsg("correlation matrix for *" +key+ "* successfully computed --> this is the first entry: " + std::to_string(C(0,0)), debug);
            return C;
        }

        //it should be identical to correlationMatrix1 but more efficient
        Eigen::MatrixXd correlationMatrix2(const Eigen::MatrixXd& S, bool debug =false, word key = "")
        {
            label nCells = mesh_.nCells(); //number of cell centers
            label Ns = S.cols();     //number of snapshots
            label nDof = S.rows();   //number of rows of snapshots matrix
            Eigen::MatrixXd C(Ns,Ns);
            label nComponents = nDof/nCells;
            debugMsg("This is the number of components of " + key + " for the correlation matrix " + std::to_string(nComponents), debug);
            //creates the volume weight for the inner product --> C = 1/Ns * S^TWS
            Eigen::VectorXd W(nComponents * nCells);
            for (int c = 0; c < nCells; ++c) {
                for (int j = 0; j < nComponents; ++j)
                    W[nComponents*c + j] = mesh_.V()[c];
            }

            C.noalias() = S.transpose() * (S.array().colwise() * W.array()).matrix();
            //C.array() /= double(Ns);

            if(debug)
            {
                double relAsym = (C - C.transpose()).norm() / std::max(1e-16, C.norm());
                if (relAsym > 1e-10) std::cout<< "Warning: noticeable asymmetry: " << relAsym << nl;
            }
            //ensures that C is symmetric also numerically (just numerical artifact to make sure everything works)
            C = C.selfadjointView<Eigen::Lower>();
            if(debug)
            {
                double relAsym = (C - C.transpose()).norm() / std::max(1e-16, C.norm());
                if (relAsym > 1e-10) std::cout<< "Warning: noticeable asymmetry still in place after correction: " << relAsym << nl;
            }

            debugMsg("correlation matrix for " +key+ " successfully computed --> this is the first entry: " + std::to_string(C(0,0)), debug);
            return C;
        }

        //computes SVD and returns eigenvectors as well as singular values already truncated by modes
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> SVDTruncation(const Eigen::MatrixXd& C, label nModes, bool debug = false, word key = "")
        {

            label Ns = C.rows();
            Eigen::VectorXd EW(Ns);     //eigenValues
            Eigen::MatrixXd EV(Ns, Ns);     //eigenVectors
            
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
            es.compute(C, Eigen::ComputeEigenvectors);
            if (es.info() != Eigen::Success) { FatalError << "Eigen Problem could not be properly computed..." << exit(FatalError);}
            //extrapolates eigenvalue and stores it from largest to smallest 
            EW = es.eigenvalues().reverse().eval();
            // extrapolates eigenvectors and stores them --> reverses columns to match eigenvalues reversing
            EV = es.eigenvectors().rowwise().reverse().eval();
            debugMsg("eigen problem of *" + key + "* successfully computed!", debug);
            plotVector(EW, "index", "eigenValues", "EigenValues Plot (log scale)");

            //corrects the numerical errors of having EW that are very close to 0 but negative
            double relTol = 1e-12;
            double tol = relTol * std::max(1.0, std::abs(EW[0]));
            Eigen::ArrayXd la = EW.array();
            EW = ((la > -tol).select(la.max(0.0), la) ).matrix();

            //truncates to decided number of modes
            Eigen::VectorXd EWTrunc = EW.head(nModes).eval();
            Eigen::MatrixXd EVTrunc = EV.leftCols(nModes).eval();
            //compute singular values and plot them
            //Eigen::VectorXd sigma = EWTrunc.cwiseSqrt();

            debugMsg("returning " + key + " truncated EW with dimensions: (" + std::to_string(EWTrunc.rows()) + ", " + std::to_string(EWTrunc.cols()) 
                    + ") and truncated V with dimensions: (" + std::to_string(EVTrunc.rows()) + ", " + std::to_string(EVTrunc.cols()) + ")", debug);
            return {EWTrunc, EVTrunc}; 
        }

        //computes SVD without truncation
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> SVD(const Eigen::MatrixXd& C, label nModes, bool debug = false, word key = "")
        {

            label Ns = C.rows();
            Eigen::VectorXd EW(Ns);     //eigenValues
            Eigen::MatrixXd EV(Ns, Ns);     //eigenVectors
            
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
            es.compute(C, Eigen::ComputeEigenvectors);
            if (es.info() != Eigen::Success) { FatalError << "Eigen Problem could not be properly computed..." << exit(FatalError);}
            //extrapolates eigenvalue and stores it from largest to smallest 
            EW = es.eigenvalues().reverse().eval();
            // extrapolates eigenvectors and stores them --> reverses columns to match eigenvalues reversing
            EV = es.eigenvectors().rowwise().reverse().eval();
            debugMsg("eigen problem of *" + key + "* successfully computed!", debug);

            //corrects the numerical errors of having EW that are very close to 0 but negative
            double relTol = 1e-12;
            double tol = relTol * std::max(1.0, std::abs(EW[0]));
            Eigen::ArrayXd la = EW.array();
            EW = ((la > -tol).select(la.max(0.0), la) ).matrix();
            plotVector(EW, "index", "eigenValues", "EigenValues Plot of " + key + " (log scale)");
            return {EW, EV};
        }

        //computing modes with phi_i = 1/lamda*Ns * (sum(uj*Qij)) with lambda being the eigenvalues and Q the eigenvectors and u the snapshots of matrix S
        Eigen::MatrixXd getModes(const Eigen::MatrixXd& Q, const Eigen::VectorXd& lambda, Eigen::MatrixXd S, bool debug = false, word key = "")
        {
            label nRows = S.rows();
            label Ns = Q.rows();
            Eigen::MatrixXd modes(nRows, Ns);
            //scalar factor1 = 1/double(Ns);

            //necessary to handle division by 0 of factor2
            const double relTol = 1e-12;
            const double scale  = std::max(1.0, lambda.cwiseAbs().maxCoeff());
            const double tol    = relTol * scale;

            for(label i = 0; i < Ns; ++i)
            {
                if (std::abs(lambda[i]) <= tol) 
                {
                    // λ_i ≈ 0 → skip/zero this mode to avoid 1/0
                    modes.col(i).setZero();
                    debugMsg("Mode " + name(i) + " zeroed (lambda≈0)", debug);
                    continue;
                }
                //scalar factor2 = factor1/std::sqrt(lambda[i]);
                scalar factor2 = 1.0/std::sqrt(lambda[i]); //if Ns is not in the correlation matrix no need factor 1
                for(label rowI = 0; rowI < nRows; ++rowI)
                {
                    scalar value = 0.0;
                    for(label j = 0; j < Ns; ++j)
                    {
                        value += S(rowI,j) * Q(j,i);
                    }
                    modes(rowI,i) = factor2 * value;
                }
            }
            debugMsg("mode of " + key + " computed. First entry has value: " + std::to_string(modes(0,0)), debug);
            return modes;
        }

        //computing modes with phi_i = 1/lamda*Ns * (sum(uj*Qij)) with lambda being the eigenvalues and Q the eigenvectors and u the snapshots of matrix S --> should be more efficient
        /******************************************************************************************/
        /******************************************************************************************/
        Eigen::MatrixXd getModes2(const Eigen::MatrixXd& Q, const Eigen::VectorXd& lambda, Eigen::MatrixXd S, bool debug = false, word key = "")
        {
            //const int Ns = Q.rows();
            const double relTol = 1e-12;
            const double tol = relTol * std::max(1.0, lambda.cwiseAbs().maxCoeff());

            // build per-mode scale: 1/(Ns*lambda_i) if |lambda_i|>tol, else 0
            //Eigen::VectorXd scale = lambda.unaryExpr([&](double li){return (std::abs(li) > tol) ? 1.0 / (double(Ns) * li) : 0.0;});
            Eigen::VectorXd scale = lambda.unaryExpr([&](double li){if (li <= tol) return 0.0; return 1.0 / std::sqrt(li); }); // --> this should be correct with sqrt(lambdai)
            //Eigen::MatrixXd modes = S * Q.transpose() * scale.asDiagonal();
            Eigen::MatrixXd modes = S * Q * scale.asDiagonal(); // --> this should be correct since it gives orthonormality
            //debugMsg("mode of " + key + " computed with method 2. First entry has value: " + std::to_string(modes(0,0)), debug);
            
            //check if modes are orthonormal
            if(debug && (key == "velocity" || key == "pressure"))
            {
                debugMsg("Checking if modes of " + key + " are orthonormal", debug);
                Eigen::MatrixXd ortCheck = L2InnerProductEigen(modes, modes, false, key);
                label nRows = std::min<int>(ortCheck.rows(), 5);
                label nCols = std::min<int>(ortCheck.cols(), 5);
                std::cout << "this are the first few entries of <" << key << "_i, " << key << "_j>_l2 " << nl 
                          << ortCheck.topLeftCorner(nRows, nCols) << nl <<nl;
            }
            return modes;
        }

        //L2 inner product eigen S^TWS --> can't use for flux for now (should not be necessary anyway)
        Eigen::MatrixXd L2InnerProductEigen(const Eigen::MatrixXd& M1, const Eigen::MatrixXd& M2, bool debug =false, word key = "")
        {
            label nCells = mesh_.nCells(); //number of cell centers
            label nCols = M1.cols();     //number of snapshots
            label nRows = M1.rows();   //number of rows of snapshots matrix
            Eigen::MatrixXd product(nCols, nCols);
            label nComponents = nRows/nCells;
            debugMsg("number of components is: " + name(nComponents), debug);
            //creates the volume weight for the inner product --> C = 1/Ns * S^TWS
            Eigen::VectorXd W(nComponents * nCells);
            for (int c = 0; c < nCells; ++c) {
                for (int j = 0; j < nComponents; ++j)
                    W[nComponents*c + j] = mesh_.V()[c];
            }

            product.noalias() = M2.transpose() * (M1.array().colwise() * W.array()).matrix();
            debugMsg("L2 inner product of " + key + " computed!", debug);
            return product;
        }

        

}; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
