
//read fields stored in the offline stage
/* RectangularMatrix<scalar> B = io.readFoamMatrixFromDisk<scalar>(offlinePath/"diffusiveTerm", debug);
RectangularMatrix<scalar> C = io.readFoamMatrixFromDisk<scalar>(offlinePath/"convectiveTerm", debug);
RectangularMatrix<scalar> extra1 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra1", debug);
RectangularMatrix<scalar> extra2 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra2", debug);
RectangularMatrix<scalar> extra3 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra3", debug);
RectangularMatrix<scalar> extra4 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra4", debug);
RectangularMatrix<scalar> extra5 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra5", debug);
RectangularMatrix<scalar> A = io.readFoamMatrixFromDisk<scalar>(offlinePath/"pressureTerm", debug); */

B = io.readFoamMatrixFromDisk<scalar>(offlinePath/"diffusiveTerm", debug);
C = io.readFoamMatrixFromDisk<scalar>(offlinePath/"convectiveTerm", debug);
extra1 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra1", debug);
extra2 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra2", debug);
extra3 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra3", debug);
extra4 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra4", debug);
extra5 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra5", debug);
A = io.readFoamMatrixFromDisk<scalar>(offlinePath/"pressureTerm", debug);
D = io.readFoamMatrixFromDisk<scalar>(offlinePath/"D", debug);
F = io.readFoamMatrixFromDisk<scalar>(offlinePath/"F", debug);
extra6 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra6", debug);
extra7 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra7", debug);
extra8 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra8", debug);
extra9 = io.readFoamMatrixFromDisk<scalar>(offlinePath/"extra9", debug);
// If you saved them:
E1fo = io.readFoamMatrixFromDisk<scalar>(offlinePath/"E1_firstOrder", debug); // Np x 1
F1fo = io.readFoamMatrixFromDisk<scalar>(offlinePath/"F1_firstOrder", debug); // Np x Nf (or Np x Nu if you built that way)
F2fo = io.readFoamMatrixFromDisk<scalar>(offlinePath/"F2_firstOrder", debug); // Np x Nu


//create residual (no time term)
ex = - nu * normed((B*rom::a)) + normed((rom::a * C * rom::a)) + rom::constv(extra1) + (extra2 * rom::a) + 
    (extra3 * rom::a) - nu * rom::constv(extra4) + rom::constv(extra5) + (A*rom::b);

utils.debugMsg(debug, "first residual correct");

ex2 =
    (-1*(D*rom::b)        )  // -D b
  + (rom::a * F * rom::a  )  // a^T F a
  + (F2fo * rom::a        )  // + F2 a
  + (F1fo * rom::a        )  // + F1 a
  +  rom::constv(E1fo     )  // + E1
  -  rom::constv(extra6   ); // - ∫ χ_i Δ p̄  dΩ   (your extra6)


Eigen::MatrixXd De(D.m(), D.n());
for (label i = 0; i < D.m(); ++i)
    for (label j = 0; j < D.n(); ++j)
        De(i,j) = D[i][j];

Eigen::JacobiSVD<Eigen::MatrixXd> svd(De);
Info << "D matrix singular values:" << nl;
for (label i = 0; i < min(5, (label)svd.singularValues().size()); ++i) {
    Info << "  σ[" << i << "] = " << svd.singularValues()(i) << nl;
}
Info << "  σ[last] = " << svd.singularValues()(svd.singularValues().size()-1) << nl;
Info << "Condition number: " << svd.singularValues()(0) / 
        svd.singularValues()(svd.singularValues().size()-1) << nl;
