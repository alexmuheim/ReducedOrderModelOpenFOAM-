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

#pragma once
#include "RectangularMatrix.H"
#include "List.H"
#include "IOstreams.H"      // Info, nl
#include "ResidualConstruction2eqs.H"        // your DSL header (Expr, residual, jacobian, timeTerm, timeJacobian)

//#include "ResidualConstruction2eqsTime.H"
#include <Eigen/Dense>      // used to solve J Δ = -R

namespace rom {

using Mat = Foam::RectangularMatrix<Foam::scalar>;
using Foam::label;
using Foam::scalar;

struct NewtonOptions
{
    label  maxNewtonIters = 50;
    scalar tolResidual    = 1e-8;
    scalar tolStep        = 1e-12;
    scalar relaxStart     = 1.0;
    label  maxBacktrack   = 8;
    bool   printInfo      = true;
};

// forward decls from your DSL (if not already in ResidualConstruction.H)
/* struct Expr;
Mat residual(const Expr&, const Mat&, label);
Mat jacobian(const Expr&, const Mat&, label);
Mat timeTerm(const Mat& a, const Mat& aPrev, scalar dt);
Mat timeJacobian(label n, scalar dt);
 */
// --------------------------------------------------------------------
class NewtonSolver
{
public:
    // -------- 1 equation constructor --------
    NewtonSolver(const Expr& exA,
                 const Mat&  a0,
                 scalar      dt,
                 scalar      tEnd,
                 NewtonOptions opt = {})
    :
        exA_(exA), exB_(), a_(a0), b_(0,1,0.0),
        aPrev_(a0), bPrev_(0,1,0.0),
        dt_(dt), tEnd_(tEnd), opt_(opt), twoEq_(false)
    {}

    // -------- 2 equations constructor --------
    NewtonSolver(const Expr& exA, const Expr& exB,
                 const Mat&  a0,   const Mat&  b0,
                 scalar      dt,   scalar      tEnd,
                 NewtonOptions opt = {})
    :
        exA_(exA), exB_(exB), a_(a0), b_(b0),
        aPrev_(a0), bPrev_(b0),
        dt_(dt), tEnd_(tEnd), opt_(opt), twoEq_(true)
    {}

    std::pair<Mat, Mat> solve()
    {
        const label nA = a_.m();
        const label nB = twoEq_ ? b_.m() : 0;
        const label n  = nA + nB;
        //Info << "starting newton loop" << endl;
        //Mat aPrev_ = a_;

        std::pair<Mat, scalar> Rnorm = assembleResidual(a_, b_, aPrev_, bPrev_);
        Mat R   = Rnorm.first;
        scalar norm = Rnorm.second;
        
        for (label it=0; it<opt_.maxNewtonIters; ++it)
        {
            if (opt_.printInfo) Foam::Info<< "  Newton iter " << (it+1) << Foam::nl;
            
            Mat J = assembleJacobian(a_, b_);
            Mat dX = solveLinear(J, R, /*negR=*/true);

            scalar resNorm  = l2(R) / norm;
            scalar stepNorm = l2(dX);

            if (opt_.printInfo)
                Foam::Info<< "    ||R|| = " << resNorm
                            << ", ||Δ|| = " << stepNorm << Foam::nl;

            if (resNorm < opt_.tolResidual) {
                if (opt_.printInfo) Foam::Info<< "  Converged by residual\n";
                break;
            }
            if (stepNorm < opt_.tolStep*(1.0 + l2(a_))) {
                if (opt_.printInfo) Foam::Info<< "  Converged by step\n";
                break;
            }

            // Damped update with simple backtracking
            scalar alpha = opt_.relaxStart;
            Mat trialA = a_, trialB = b_;
            bool accepted = false;

            for (label ls=0; ls<opt_.maxBacktrack; ++ls)
            {
                for (label i=0;i<nA;++i) trialA[i][0] = a_[i][0] + alpha * dX[i][0];
                if (twoEq_) for (label i=0;i<nB;++i) trialB[i][0] = b_[i][0] + alpha * dX[nA + i][0];


                std::pair<Mat, scalar> Rt_normt = assembleResidual(trialA, trialB, aPrev_, bPrev_);
                Mat   Rt    = Rt_normt.first;
                scalar nt   = Rt_normt.second;

                if (l2(Rt) < l2(R)) {
                    a_ = trialA; b_ = trialB;
                    R  = Rt;     norm = nt;
                    accepted = true; break;
                }
                alpha *= 0.5;
            }

            if (!accepted) {
                // tiny fallback step
                for (label i=0;i<nA;++i) a_[i][0] += 0.1 * dX[i][0];
                if (twoEq_) for (label i=0;i<nB;++i) b_[i][0] += 0.1 * dX[nA + i][0];
                Rnorm = assembleResidual(a_, b_, aPrev_, bPrev_);
                R     = Rnorm.first;
                norm  = Rnorm.second;
            }
        }
        aPrev_ = a_; bPrev_ = b_;
        return {a_, b_};
    }

    // Access current coefficients (after solve or mid-iteration if you like)
    const Mat& a() const { return a_; }
    const Mat& b() const { return b_; }  // only meaningful if twoEq_==true

private:
    Expr exA_, exB_;
    Mat  a_, b_;
    Mat  aPrev_, bPrev_;
    scalar dt_, tEnd_;
    NewtonOptions opt_;
    bool twoEq_;

    // ---- helpers ----------------------------------------------------------

    inline static scalar l2(const Mat& v)
    {
        scalar s=0; for (label i=0;i<v.m();++i) s += v[i][0]*v[i][0]; return std::sqrt(s);
    }
    inline static scalar l2Stack(const Mat& a, const Mat& b){
        scalar s=0; for (label i=0;i<a.m();++i) s += a[i][0]*a[i][0];
        for (label i=0;i<b.m();++i) s += b[i][0]*b[i][0]; return std::sqrt(s);
    }

    inline std::pair<Mat, scalar> assembleResidual(const Mat& aCur, const Mat& bCur,
                                            const Mat& aPrev, const Mat& bPrev) const
    {
        if (!twoEq_) {
            // 1 eq: R = time(a) + residual1(exA,a); norm from residual1
            //Info << "solving with 1 eq " << endl;
            std::pair<Mat, scalar> Rx = rom::residual1(exA_, aCur); // (Rexpr, scale)
            const label nA = aCur.m();
            Mat R(nA,1,0.0);
            Mat Rtime = rom::timeTerm(aCur, aPrev, dt_);
            for (label i=0;i<nA;++i) R[i][0] = Rtime[i][0] + Rx.first[i][0];
            return std::make_pair(R, Rx.second);
        } else {
            // 2 eqs: stack, then add time to the TOP block (a-equation).
            //Info << "solving with 2 eq " << endl;
            std::pair<Mat, Mat> RJ = rom::assemble2eq(exA_, exB_, aCur, bCur, /*rowScale=*/true);
            Mat R = RJ.first;
            Mat& Jdummy = RJ.second; (void)Jdummy; // not used here

            const label nA = aCur.m();
            Mat RtimeA = rom::timeTerm(aCur, aPrev, dt_);
            for (label i=0;i<nA;++i) R[i][0] += RtimeA[i][0];
            return std::make_pair(R, scalar(1.0)); // already row-scaled
        }
    }

    inline Mat assembleJacobian(const Mat& aCur, const Mat& bCur) const
    {
        if (!twoEq_) {
            Mat Jexpr = rom::jacobian1(exA_, aCur);
            Mat Jtime = rom::timeJacobian(aCur.m(), dt_);
            const label nA = aCur.m();
            Mat J(nA,nA,0.0);
            for (label i=0;i<nA;++i)
                for (label j=0;j<nA;++j)
                    J[i][j] = Jexpr[i][j] + Jtime[i][j];
            return J;
        } else {
            std::pair<Mat, Mat> RJ = rom::assemble2eq(exA_, exB_, aCur, bCur, /*rowScale=*/true);
            Mat J = RJ.second;

            // add time Jacobian to top-left block (a-equation)
            const label nA = aCur.m(), nB = bCur.m();
            Mat JtimeA = rom::timeJacobian(nA, dt_);
            for (label i=0;i<nA;++i)
                for (label j=0;j<nA;++j)
                    J[i][j] += JtimeA[i][j];
            return J;
        }
    }


    // Solve J * x = (negR ? -R : R) with Eigen (robust and simple)
    inline static Mat solveLinear(const Mat& J, const Mat& R, bool negR)
    {
        const label m = J.m(), n = J.n(); // here J is square (m==n)
        Eigen::MatrixXd Je(m,n);
        Eigen::VectorXd Re(m);
        for (label i=0;i<m;++i){
            Re[i] = (negR ? -R[i][0] : R[i][0]);
            for (label j=0;j<n;++j) Je(i,j) = J[i][j];
        }
        Eigen::VectorXd Xe = Je.fullPivLu().solve(Re);
        Mat X(n,1,0.0);
        for (label i=0;i<n;++i) X[i][0] = Xe[i];
        return X;
    }
};

} // namespace rom
