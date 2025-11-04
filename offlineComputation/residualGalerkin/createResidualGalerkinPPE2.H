#include "MathOffline2.H"

/**********************************************/
/******** ITHACA-FV Style - Direct Approach ***/
/**********************************************/

// Diffusive term: B[i][j] = ∫ φ_i · ∇²φ_j dV
RectangularMatrix<scalar> B = mathOff2.diffusiveMatrix(uModes, "Diffusive Term");
IO.saveFoamMatrixToDisk<scalar>(B, offlinePath/"diffusiveTerm");

// Convective term with flux modes
RectangularMatrix<scalar> C = mathOff2.convectionTensor(uModes, uModes, fModes, "Convective Term");
IO.saveFoamMatrixToDisk<scalar>(C, offlinePath/"convectiveTerm");

// Extra terms - still need pre-computed fields for these
PtrList<volVectorField> convMean = mathOff2.convectiveTermFlux(uMean, fMean, "mean convective term");
RectangularMatrix<scalar> extra1 = mathOff2.firstOrder(uModes, convMean, "extraTerm 1");
IO.saveFoamMatrixToDisk<scalar>(extra1, offlinePath/"extra1");

PtrList<volVectorField> convFlux = mathOff2.convectiveTermFlux(uMean, fModes, "flux convective term");
RectangularMatrix<scalar> extra2 = mathOff2.secondOrder(uModes, convFlux, "extraTerm 2");
IO.saveFoamMatrixToDisk<scalar>(extra2, offlinePath/"extra2");

PtrList<volVectorField> convVel = mathOff2.convectiveTermVec(uModes, fMean, "velocity convective term");
RectangularMatrix<scalar> extra3 = mathOff2.secondOrder(uModes, convVel, "extraTerm 3");
IO.saveFoamMatrixToDisk<scalar>(extra3, offlinePath/"extra3");

PtrList<volVectorField> laplUMean = mathOff2.fieldLaplacian(uMean, "laplacian of U mean");
RectangularMatrix<scalar> extra4 = mathOff2.firstOrder(uModes, laplUMean, "extraTerm 4");
IO.saveFoamMatrixToDisk<scalar>(extra4, offlinePath/"extra4");

PtrList<volVectorField> gradPMean = mathOff2.fieldGradient(pMean, "gradient of p mean");
RectangularMatrix<scalar> extra5 = mathOff2.firstOrder(uModes, gradPMean, "extraTerm 5");
IO.saveFoamMatrixToDisk<scalar>(extra5, offlinePath/"extra5");

// Pressure term: direct computation
RectangularMatrix<scalar> pressTerm = mathOff2.pressureGradientMatrix(uModes, pModes, "pressure Term");
IO.saveFoamMatrixToDisk<scalar>(pressTerm, offlinePath/"pressureTerm");

// Extra6 - keep pre-computed for single field
PtrList<volScalarField> laplPMean = mathOff2.fieldLaplacian(pMean, "laplacian of p mean");
RectangularMatrix<scalar> extra6 = mathOff2.firstOrder(pModes, laplPMean, "extraTerm 6");
IO.saveFoamMatrixToDisk<scalar>(extra6, offlinePath/"extra6");

// Extra7
PtrList<volScalarField> divConvMean = mathOff2.fieldDiv<volScalarField>(convMean, "divergence of conv mean");
RectangularMatrix<scalar> extra7 = mathOff2.firstOrder(pModes, divConvMean, "extraTerm 7");
IO.saveFoamMatrixToDisk<scalar>(extra7, offlinePath/"extra7");

// Extra8
PtrList<volScalarField> divConvVel = mathOff2.fieldDiv<volScalarField>(convVel, "divergence of conv velocity");
RectangularMatrix<scalar> extra8 = mathOff2.secondOrder(pModes, divConvVel, "extraTerm 8");
IO.saveFoamMatrixToDisk<scalar>(extra8, offlinePath/"extra8");

// Extra9
PtrList<volScalarField> divConvFlux = mathOff2.fieldDiv<volScalarField>(convFlux, "divergence of conv flux");
RectangularMatrix<scalar> extra9 = mathOff2.secondOrder(pModes, divConvFlux, "extraTerm 9");
IO.saveFoamMatrixToDisk<scalar>(extra9, offlinePath/"extra9");

// Pressure Laplacian: D[i][j] = ∫ ∇χ_i · ∇χ_j dV
RectangularMatrix<scalar> D = mathOff2.pressureLaplacian(pModes, "D Term");
IO.saveFoamMatrixToDisk<scalar>(D, offlinePath/"D");

// PPE projection tensor
RectangularMatrix<scalar> F = mathOff2.ppeConvectionTensor(pModes, uModes, fModes, "F Term");
IO.saveFoamMatrixToDisk<scalar>(F, offlinePath/"F");


//check of mass matrix
RectangularMatrix<scalar> M = mathOff2.massMatrix(uModes, "Mass Term");
IO.saveFoamMatrixToDisk<scalar>(M,offlinePath/"Mu");


/***************************************************************** */
//try boundary term 

// All boundary patches:
RectangularMatrix<scalar> N = assembleNij(pModes, uModes);
IO.saveFoamMatrixToDisk<scalar>(N,offlinePath/"N");
/* // Only on a subset, e.g. {"wall","inlet"}:
List<word> subset(2); subset[0]="wall"; subset[1]="inlet";
RectangularMatrix<scalar> Nsub = assembleNij(pModes, uModes, subset); */

