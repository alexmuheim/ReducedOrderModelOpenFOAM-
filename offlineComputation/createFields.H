//create necessary fields --> need to be ptrlists
PtrList<volVectorField> uModes; //rows --> number of modes, cols--> number of dof
PtrList<surfaceScalarField> fModes;
PtrList<volScalarField> pModes;
PtrList<volVectorField> uMean;
PtrList<surfaceScalarField> fMean;
PtrList<volScalarField> pMean;

//read modes 
uModes = IO.readPtrField_noTime<volVectorField>("UMode", runTime, mesh, 1, nModesU, "POD", debug, "UMode");
//fModes = IO.readPtrField_noTime<surfaceScalarField>("fMode", runTime, mesh, 1, nModesU, "POD", debug, "fMode");
for(label m = 0; m < uModes.size(); ++m)
{
    tmp<surfaceScalarField> tFlux = fvc::interpolate(uModes[m]) & mesh.Sf();
    fModes.append(tFlux);
}

if(residual == "galerkin" || residual  == "galerkinPressure")
{
    pModes = IO.readPtrField_noTime<volScalarField>("pMode", runTime, mesh, 1, nModesp, "POD", debug, "pMode");
}
else if(residual == "PPE")
{
    pModes = IO.readPtrField_noTime<volScalarField>("pMode2", runTime, mesh, 1, nModesp, "POD", debug, "pMode2");
}

//read means
uMean = IO.readMeanPtrField<volVectorField>("uMean", runTime, mesh, "POD/Mean", debug, "uMean");
//fMean = IO.readMeanPtrField<surfaceScalarField>("fMean", runTime, mesh, "POD/Mean", debug, "fMean");
pMean = IO.readMeanPtrField<volScalarField>("pMean", runTime, mesh, "POD/Mean", debug, "pMean");
tmp<surfaceScalarField> tFluxMean = fvc::interpolate(uMean[0]) & mesh.Sf();
fMean.append(tFluxMean);

for (int i=0;i<pModes.size();++i)
  for (int j=0;j<pModes.size();++j)
    Info<<"<ϕ_i,ϕ_j>_V = "<<gSum(pModes[i]*pModes[j]*mesh.V())<<nl;

//prova fix per pressure

/* 
const dimensionedScalar Vsum("Vsum", dimVolume, gSum(mesh.V()));
for (label j=0; j<pModes.size(); ++j)
{
    const dimensionedScalar mj = fvc::domainIntegrate(pModes[j]) / Vsum;
    if (mag(mj.value()) > 1e-14) { pModes[j] -= mj; pModes[j].correctBoundaryConditions(); }
}

for (label j=0; j<pModes.size(); ++j)
{
    const dimensionedScalar Mjj = fvc::domainIntegrate(pModes[j]*pModes[j]) / Vsum; // ≈ 2.78994e-3
    const scalar s = 1.0/std::sqrt(std::max(SMALL, Mjj.value()));
    pModes[j] *= s;      // now <pj,pj>_W ≈ 1
}

for (label i = 0; i < pModes.size(); ++i)
for (label j = 0; j < pModes.size(); ++j)
{
    dimensionedScalar Mij = fvc::domainIntegrate( pModes[i]*pModes[j] ) / Vsum;
    if ( (i==j && mag(Mij.value()-1.0) > 1e-3) ||
         (i!=j && mag(Mij.value())      > 1e-3) )
        Info<< "Inner(p"<<i<<",p"<<j<<") = "<< Mij.value() << nl;
}


for (label j=0; j<pModes.size(); ++j)
{
    // subtract projections on previous modes
    for (label i=0; i<j; ++i)
    {
        const dimensionedScalar Mij = fvc::domainIntegrate(pModes[i]*pModes[j]) / Vsum;
        pModes[j] -= Mij.value()*pModes[i];
    }
    // renormalise
    const dimensionedScalar Mjj = fvc::domainIntegrate(pModes[j]*pModes[j]) / Vsum;
    pModes[j] *= 1.0/std::sqrt(std::max(SMALL, Mjj.value()));
}

for (label i=0; i<pModes.size(); ++i)
{
  for (label j=i; j<pModes.size(); ++j)
  {
    const dimensionedScalar Mij = fvc::domainIntegrate(pModes[i]*pModes[j]) / Vsum;
    if (i==j) Info<< "Inner(p"<<i<<",p"<<i<<")="<< Mij.value() << nl;
    else if (mag(Mij.value())>1e-3) Info<< "Offdiag(p"<<i<<",p"<<j<<")="<< Mij.value() << nl;
  }
} */