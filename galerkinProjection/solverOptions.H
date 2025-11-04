//solver options

/* rom::NewtonOptionsTime opt2;
opt2.maxNewtonIters = maxIter;
opt2.tolResidual = tol; */
//scalar tolStep = 1e-12;
//scalar relaxStart = 1.0;
//label maxBacktrack = 8;
//opt2.printInfo = printInfos;



rom::NewtonOptions opt;

opt.maxNewtonIters = maxIter;
opt.tolResidual = tol;
//scalar tolStep = 1e-12;
//scalar relaxStart = 1.0;
//label maxBacktrack = 8;
opt.printInfo = printInfos;

scalar tEnd = endTime;
label Nt = std::floor(endTime / dt);
label nSaves = Nt / writeInterval + 1;


/* //solvertimenewton
RectangularMatrix<scalar> timeCoeff(nModesU, nSaves, 0.0); //solution will be stored in this matrix togheter withh initial condition

rom::TrajectorySink trj{ &timeCoeff, writeInterval, nSaves};
 */


