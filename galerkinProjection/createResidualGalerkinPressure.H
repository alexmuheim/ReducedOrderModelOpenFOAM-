
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

//create residual (no time term)
ex = - nu * normed((B*rom::a)) + normed((rom::a * C * rom::a)) + rom::constv(extra1) + (extra2 * rom::a) + 
    (extra3 * rom::a) - nu * rom::constv(extra4) + rom::constv(extra5) + (A*rom::a);

