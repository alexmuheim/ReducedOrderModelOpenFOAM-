#include "MathOffline.H"

/**********************************************/
/******** Residual for same coeff ROM  ********/
/**********************************************/

// matrix for the residual --> equation dependent --> wrap into header file
PtrList<volVectorField> laplUModes = mathOff.fieldLaplacian(uModes,"laplacian of U modes");
RectangularMatrix<scalar> B = mathOff.secondOrder(uModes, laplUModes, "Diffusive Term");
IO.saveFoamMatrixToDisk<scalar>(B,offlinePath/"diffusiveTerm");

//non linear term with flux modes
RectangularMatrix<scalar> C = mathOff.thirdOrder<PtrList<volVectorField>, PtrList<surfaceScalarField>>(uModes, uModes, fModes, "Convective Term");
IO.saveFoamMatrixToDisk<scalar>(C,offlinePath/"convectiveTerm");

// <phi_j, umean.grad(u mean)>
PtrList<volVectorField> convMean = mathOff.convectiveTermFlux(uMean, fMean, "mean convective term");
RectangularMatrix<scalar> extra1 = mathOff.firstOrder(uModes, convMean, "extraTerm 1");
IO.saveFoamMatrixToDisk<scalar>(extra1,offlinePath/"extra1");

// <phii, umean grad phi_j>
PtrList<volVectorField> convFlux = mathOff.convectiveTermFlux(uMean, fModes, "flux convective term");
RectangularMatrix<scalar> extra2 = mathOff.secondOrder(uModes, convFlux, "extraTerm 2");
IO.saveFoamMatrixToDisk<scalar>(extra2,offlinePath/"extra2");

// <phi_i,phi_i.gradumean>
PtrList<volVectorField> convVel = mathOff.convectiveTermVec(uModes, fMean, "velocity convective term");
RectangularMatrix<scalar> extra3 = mathOff.secondOrder(uModes, convVel, "extraTerm 3");
IO.saveFoamMatrixToDisk<scalar>(extra3,offlinePath/"extra3");

// <phi_i, lapl umean>
PtrList<volVectorField> laplUMean = mathOff.fieldLaplacian(uMean, "laplacian of U mean");
RectangularMatrix<scalar> extra4 = mathOff.firstOrder(uModes, laplUMean, "extraTerm 4");
IO.saveFoamMatrixToDisk<scalar>(extra4, offlinePath/"extra4");


// <phi_i, grad pmean>
PtrList<volVectorField> gradPMean = mathOff.fieldGradient(pMean, "gradient of p mean");
RectangularMatrix<scalar> extra5 = mathOff.firstOrder(uModes, gradPMean, "extraTerm 5");
IO.saveFoamMatrixToDisk<scalar>(extra5, offlinePath/"extra5");

//pressure term
PtrList<volVectorField> gradPModes = mathOff.fieldGradient(pModes, "gradient of p modes");
RectangularMatrix<scalar> pressTerm = mathOff.secondOrder(uModes, gradPModes, "pressure Term");
IO.saveFoamMatrixToDisk<scalar>(pressTerm, offlinePath/"pressureTerm");




