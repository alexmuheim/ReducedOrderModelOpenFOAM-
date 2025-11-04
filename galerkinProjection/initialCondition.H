//creating the initial condition


RectangularMatrix<scalar> a0(nModesU, 1, 0.0);
RectangularMatrix<scalar> b0(nModesp, 1, 0.0);
//pMean = dimensionedScalar("zero", pMean.dimensions(), 0.0);
volVectorField uFluct0 = u0 - uMean;
volScalarField pFluct0 = p0 - pMean;


for(label l = 0; l < nModesU; ++l)
{
    a0[l][0] = fvc::domainIntegrate(uModes[l] & uFluct0).value();
}

for(label l = 0; l < nModesp; ++l)
{
    b0[l][0] = fvc::domainIntegrate(pModes[l] * pFluct0).value();
}

/* RectangularMatrix<scalar> rhs(nModesp,1,0.0);
RectangularMatrix<scalar> quad = (a0 * F * a0);   // size Np×1

for (label i=0;i<nModesp;++i)
    {
    //rhs[i][0] = - extra7[i][0] + quad[i][0] + (extra8 * a0)[i][0] + (extra9 * a0)[i][0];
        rhs[i][0] = - extra7[i][0] + extra6[i][0] - quad[i][0] - (extra8 * a0)[i][0] - (extra9 * a0)[i][0];
    }
b0 = solveSmallLinearSPD(D, rhs); */


//Info << "initial condition " << a0 << endl;

