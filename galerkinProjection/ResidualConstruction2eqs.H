/*---------------------------------------------------------------------------*\
  Auto-dim ROM expression system (2 vars, 2 eqs, quad only in 'a')
  - Unknowns: a ∈ R^{nA}, b ∈ R^{nB}
  - One Expr = one equation (its residual length is inferred)
  - Linear terms can use 'a' or 'b' with rectangular operators
  - Quadratic terms are ONLY a ⊗ a (convective-like)
  - Jacobian row-blocks: ∂R/∂a (nOut×nA), ∂R/∂b (nOut×nB)
  - Assembly stacks two equations into [Ra; Rb] and block J
\*---------------------------------------------------------------------------*/
#pragma once
#include "RectangularMatrix.H"
// (Optional) #include "error.H"  // For FatalErrorInFunction, if desired

namespace rom {

using Mat = Foam::RectangularMatrix<Foam::scalar>;
using Foam::label; using Foam::scalar;

// ------------------- variable tags -------------------
struct Avar{}; inline constexpr Avar a{};
struct Bvar{}; inline constexpr Bvar b{};

enum class Var { a, b };

// ------------------- term types ----------------------
struct Lin  {
    const Mat* M = nullptr;   // nOut × nVar
    scalar     k = 1.0;
    bool       inNorm = false;
    Var        v = Var::a;    // which vector: a or b
};

struct Quad {
    const Mat* C = nullptr;   // nOut × (nA * nA)  (ONLY a⊗a)
    scalar     k = 1.0;
    bool       inNorm = false;
};

struct Cst  {
    const Mat* b = nullptr;   // nOut × 1
    scalar     k = 1.0;
};

struct Expr {
    Foam::List<Lin>  lins;
    Foam::List<Quad> quads;
    Foam::List<Cst>  csts;

    inline void add(const Lin&  t){ lins.append(t);  }
    inline void add(const Quad& t){ quads.append(t); }
    inline void add(const Cst&  t){ csts.append(t);  }
};

// ------------------- helpers & sugar -----------------
inline Cst constv(const Mat& b) { return Cst{&b, 1.0}; }

// linear tagging
inline Lin operator*(const Mat& M, Avar){ return Lin{&M, 1.0, false, Var::a}; }
inline Lin operator*(const Mat& M, Bvar){ return Lin{&M, 1.0, false, Var::b}; }

// quadratic ONLY in a: (a * C * a)
struct ACa { const Mat* C; scalar k = 1.0; };
inline ACa  operator*(Avar, const Mat& C){ return ACa{&C, 1.0}; }
inline Quad operator*(ACa ac, Avar)      { return Quad{ac.C, ac.k, false}; }

// flags + scalar scaling
inline Lin  normed(Lin  t){ t.inNorm = true; return t; }
inline Quad normed(Quad t){ t.inNorm = true; return t; }

inline Lin  operator*(scalar s, Lin  t){ t.k *= s; return t; }
inline Quad operator*(scalar s, Quad t){ t.k *= s; return t; }
inline Cst  operator*(scalar s, Cst  t){ t.k *= s; return t; }

// merging (Expr + Term, Term + Expr, Term + Term)
inline void appendInto(Expr& ex, const Lin&  t){ ex.add(t); }
inline void appendInto(Expr& ex, const Quad& t){ ex.add(t); }
inline void appendInto(Expr& ex, const Cst&  t){ ex.add(t); }

template<class Term>
inline Expr operator+(Expr ex, const Term& t){ appendInto(ex,t); return ex; }

template<class Term>
inline Expr operator+(const Term& t, Expr ex){ appendInto(ex,t); return ex; }

template<class T1, class T2>
inline Expr operator+(const T1& t1, const T2& t2){ Expr ex; appendInto(ex,t1); appendInto(ex,t2); return ex; }

template<class Term>
inline Expr operator-(Expr ex, const Term& t){ return ex + (Foam::scalar(-1)*t); }

// ------------------- time helpers (optional) ----------
inline Mat timeTerm(const Mat& a, const Mat& aPrev, scalar dt)
{
    const scalar invDt = 1.0/dt;
    const label n = a.m();
    Mat r(n,1,0.0);
    for (label i=0;i<n;++i) r[i][0] = invDt*(a[i][0]-aPrev[i][0]);
    return r;
}
inline Mat timeJacobian(label n, scalar dt)
{
    Mat J(n,n,0.0);
    const scalar invDt = 1.0/dt;
    for (label i=0;i<n;++i) J[i][i] = invDt;
    return J;
}

// ------------------- norms & checks -------------------
inline scalar l2(const Mat& v)
{
    scalar s = 0.0;
    for (label i=0;i<v.m();++i){ const scalar vi = v[i][0]; s += vi*vi; }
    return std::sqrt(s);
}

inline label deduceRows(const Expr& ex, label fallback)
{
    if (ex.csts.size()) return ex.csts[0].b->m();
    if (ex.lins.size()) return ex.lins[0].M->m();
    if (ex.quads.size()) return ex.quads[0].C->m();
    return fallback; // e.g. a.m() or b.m()
}

inline void assertShape(bool ok, const char* msgM,
                        Foam::label Mr, Foam::label Mc,
                        Foam::label Er, Foam::label Ec)
{
    if (!ok)
    {
        FatalErrorInFunction                 // <-- no Foam::
            << msgM << " shape " << Mr << "×" << Mc
            << " but expected " << Er << "×" << Ec << Foam::nl
            << Foam::abort(Foam::FatalError);
    }
}

// ------------------- one-equation residual ------------
/** Evaluate residual of ONE equation:
 *  - Automatically deduces nOut (equation length) from the first term found.
 *  - Linear terms may be rectangular (nOut×nA) or (nOut×nB).
 *  - Quadratic terms are only a⊗a with C of size nOut×(nA*nA).
 *  Returns: (R, scale) where scale = ||conv|| + ||diff|| for row scaling.
 */
inline std::pair<Mat, scalar>
residual(const Expr& ex, const Mat& avec, const Mat& bvec)
{
    const label nA = avec.m();
    const label nB = bvec.m();

    const label nOut = deduceRows(ex, /*fallback*/ nA ? nA : nB);
    Mat r(nOut, 1, 0.0), convVec(nOut,1,0.0), diffVec(nOut,1,0.0);

    // constants
    for (const auto& c: ex.csts){
        const Mat& B = *c.b;
        assertShape(B.n()==1, "Cst (b) must be column", B.m(), B.n(), nOut, 1);
        assertShape(B.m()==nOut, "Cst row mismatch", B.m(), B.n(), nOut, 1);
        for (label i=0;i<nOut;++i) r[i][0] += c.k * B[i][0];
    }

    // linear: r += k * M * (a or b)
    for (const auto& L: ex.lins){
        const Mat& M = *L.M;
        const bool usesA = (L.v == Var::a);
        const label nVar = usesA ? nA : nB;

        assertShape(M.m()==nOut, "Lin row mismatch", M.m(), M.n(), nOut, M.n());
        assertShape(M.n()==nVar, "Lin col mismatch", M.m(), M.n(), M.m(), nVar);

        for (label i=0;i<nOut;++i){
            scalar acc = 0.0;
            for (label j=0;j<nVar;++j)
                acc += M[i][j] * (usesA ? avec[j][0] : bvec[j][0]);
            const scalar c = L.k * acc;
            r[i][0] += c;
            if (L.inNorm) diffVec[i][0] += c;
        }
    }

    // quadratic ONLY in a: r_l += k * Σ_{i,j} C[l, i*nA + j] a_i a_j
    for (const auto& Q: ex.quads){
        const Mat& C = *Q.C;
        assertShape(C.m()==nOut, "Quad row mismatch", C.m(), C.n(), nOut, C.n());
        assertShape(C.n()==nA*nA, "Quad col mismatch (need nA*nA)", C.m(), C.n(), C.m(), nA*nA);

        for (label l=0;l<nOut;++l){
            scalar acc = 0.0;
            for (label i=0;i<nA;++i)
                for (label j=0;j<nA;++j)
                    acc += C[l][i*nA + j] * avec[i][0] * avec[j][0];
            const scalar c = Q.k * acc;
            r[l][0] += c;
            if (Q.inNorm) convVec[l][0] += c;
        }
    }

    const scalar eps = 1e-300;
    const scalar scale = Foam::max(l2(convVec) + l2(diffVec), eps);
    return {r, scale};
}

// ------------------- one-equation Jacobian row blocks -
/** Build Jacobian row-blocks for ONE equation:
 *  J_wrtA: nOut × nA,   J_wrtB: nOut × nB
 *  Quadratic (a⊗a) contributes only to J_wrtA.
 */
inline std::pair<Mat, Mat>
jacobianRow(const Expr& ex, const Mat& avec, const Mat& bvec)
{
    const label nA = avec.m();
    const label nB = bvec.m();
    const label nOut = deduceRows(ex, /*fallback*/ nA ? nA : nB);

    Mat JwrtA(nOut, nA, 0.0), JwrtB(nOut, nB, 0.0);

    // linear contributions
    for (const auto& L: ex.lins){
        const Mat& M = *L.M;
        const bool usesA = (L.v == Var::a);
        const label nVar = usesA ? nA : nB;

        assertShape(M.m()==nOut, "Lin row mismatch (J)", M.m(), M.n(), nOut, M.n());
        assertShape(M.n()==nVar, "Lin col mismatch (J)", M.m(), M.n(), M.m(), nVar);

        if (usesA){
            for (label i=0;i<nOut;++i)
                for (label j=0;j<nA;++j)
                    JwrtA[i][j] += L.k * M[i][j];
        } else {
            for (label i=0;i<nOut;++i)
                for (label j=0;j<nB;++j)
                    JwrtB[i][j] += L.k * M[i][j];
        }
    }

    // quadratic a⊗a: ∂/∂a_m sum_{i,j} C[l, i*nA + j] a_i a_j
    // = sum_j C[l, m*nA + j] a_j  +  sum_i C[l, i*nA + m] a_i
    for (const auto& Q: ex.quads){
        const Mat& C = *Q.C;
        assertShape(C.m()==nOut, "Quad row mismatch (J)", C.m(), C.n(), nOut, C.n());
        assertShape(C.n()==nA*nA, "Quad col mismatch (J)", C.m(), C.n(), C.m(), nA*nA);

        for (label l=0;l<nOut;++l){
            for (label m=0;m<nA;++m){
                scalar term1 = 0.0, term2 = 0.0;
                for (label j=0;j<nA;++j) term1 += C[l][m*nA + j] * avec[j][0];
                for (label i=0;i<nA;++i) term2 += C[l][i*nA + m] * avec[i][0];
                JwrtA[l][m] += Q.k * (term1 + term2);
            }
        }
    }

    return {JwrtA, JwrtB};
}

// ------------------- two-equation assembly ------------
/** Assemble stacked residual and Jacobian for TWO equations:
 *  exA: builds R_a (momentum-like)
 *  exB: builds R_b (continuity-like)
 *  If rowScale=true, each equation row is divided by its own scale.
 *
 *  Returns: { R (nA+nB × 1), J (nA+nB × nA+nB) }
 */
inline std::pair<Mat, Mat>
assemble2eq(const Expr& exA, const Expr& exB,
            const Mat& avec, const Mat& bvec,
            bool rowScale = true)
{
    const label nA = avec.m();
    const label nB = bvec.m();

    // Per-equation residuals + scales (auto nOut)
    auto [Ra, sA] = residual(exA, avec, bvec);
    auto [Rb, sB] = residual(exB, avec, bvec);

    const label mA = Ra.m();   // may be != nA if you design it that way
    const label mB = Rb.m();   // may be != nB likewise

    // Per-equation Jacobian row-blocks
    auto [Jaa, Jab] = jacobianRow(exA, avec, bvec); // sizes mA×nA, mA×nB
    auto [Jba, Jbb] = jacobianRow(exB, avec, bvec); // sizes mB×nA, mB×nB

    /* // Optional row scaling
    if (rowScale){
        const scalar eps = 1e-300;
        sA = Foam::max(sA, eps);
        sB = Foam::max(sB, eps);

        for (label i=0;i<mA;++i) Ra[i][0] /= sA;
        for (label i=0;i<mB;++i) Rb[i][0] /= sB;

        for (label i=0;i<mA;++i){ for (label j=0;j<nA;++j) Jaa[i][j] /= sA; for (label j=0;j<nB;++j) Jab[i][j] /= sA; }
        for (label i=0;i<mB;++i){ for (label j=0;j<nA;++j) Jba[i][j] /= sB; for (label j=0;j<nB;++j) Jbb[i][j] /= sB; }
    } */

    // Stack
    Mat R(mA + mB, 1, 0.0);
    Mat J(mA + mB, nA + nB, 0.0);

    // R
    for (label i=0;i<mA;++i) R[i][0]     = Ra[i][0];
    for (label i=0;i<mB;++i) R[mA + i][0]= Rb[i][0];

    // J top (A-equation rows)
    for (label i=0;i<mA;++i){
        for (label j=0;j<nA;++j) J[i][j]      = Jaa[i][j];
        for (label j=0;j<nB;++j) J[i][nA + j] = Jab[i][j];
    }
    // J bottom (B-equation rows)
    for (label i=0;i<mB;++i){
        for (label j=0;j<nA;++j) J[mA + i][j]      = Jba[i][j];
        for (label j=0;j<nB;++j) J[mA + i][nA + j] = Jbb[i][j];
    }

    return {R, J};
}

// ------------------- convenience (1-var fallback) -----
inline std::pair<rom::Mat, Foam::scalar>
residual1(const rom::Expr& ex, const rom::Mat& a)
{ rom::Mat b(0,1,0.0); return rom::residual(ex, a, b); }

inline rom::Mat jacobian1(const rom::Expr& ex, const rom::Mat& a)
{ rom::Mat b(0,1,0.0); return rom::jacobianRow(ex, a, b).first; }

// Function to check the Jacobian using finite differences
void checkJacobian2(const rom::Expr& ex2, 
                    const rom::Mat& a, 
                    const rom::Mat& b,
                    double epsilon = 1e-7) 
{
    const label nA = a.m();
    const label nB = b.m();
    
    // Get analytical Jacobian from your ROM system
    auto [Jba_analytical, Jbb_analytical] = rom::jacobianRow(ex2, a, b);
    
    // Compute numerical Jacobian using finite differences
    auto [R0, scale] = rom::residual(ex2, a, b);
    const label nOut = R0.m();
    
    // Check ∂R/∂a (Jacobian w.r.t. 'a')
    rom::Mat Jba_numerical(nOut, nA, 0.0);
    for (label j = 0; j < nA; ++j) {
        rom::Mat a_plus = a;
        a_plus[j][0] += epsilon;
        
        auto [R_plus, _] = rom::residual(ex2, a_plus, b);
        
        for (label i = 0; i < nOut; ++i) {
            Jba_numerical[i][j] = (R_plus[i][0] - R0[i][0]) / epsilon;
        }
    }
    
    // Check ∂R/∂b (Jacobian w.r.t. 'b')
    rom::Mat Jbb_numerical(nOut, nB, 0.0);
    for (label j = 0; j < nB; ++j) {
        rom::Mat b_plus = b;
        b_plus[j][0] += epsilon;
        
        auto [R_plus, _] = rom::residual(ex2, a, b_plus);
        
        for (label i = 0; i < nOut; ++i) {
            Jbb_numerical[i][j] = (R_plus[i][0] - R0[i][0]) / epsilon;
        }
    }
    
    // Compute errors
    scalar maxErrorA = 0.0, maxErrorB = 0.0;
    for (label i = 0; i < nOut; ++i) {
        for (label j = 0; j < nA; ++j) {
            scalar err = std::abs(Jba_analytical[i][j] - Jba_numerical[i][j]);
            maxErrorA = std::max(maxErrorA, err);
        }
        for (label j = 0; j < nB; ++j) {
            scalar err = std::abs(Jbb_analytical[i][j] - Jbb_numerical[i][j]);
            maxErrorB = std::max(maxErrorB, err);
        }
    }
    
    Info << "Jacobian check for ex2:" << nl
         << "  Max error in ∂R/∂a: " << maxErrorA << nl
         << "  Max error in ∂R/∂b: " << maxErrorB << nl;
    
    if (maxErrorA < 1e-5 && maxErrorB < 1e-5) {
        Info << "  ✓ Jacobian is CORRECT (errors < 1e-5)" << nl;
    } else {
        Info << "  ✗ Jacobian may be INCORRECT (errors > 1e-5)" << nl;
    }
}


} // namespace rom
