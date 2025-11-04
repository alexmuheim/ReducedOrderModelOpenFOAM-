/*********************************************/
/*********************************************/
/************Creates helper fields************/
/*********************************************/
/*********************************************/

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
scalar nu = readScalar(transportProperties.lookup("nu"));
Info << "viscosity is: nu = " << nu << endl;

volVectorField uMean
(
    IOobject
    (
        "uMean",
        runTime.path()/"POD"/"Mean",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField pMean
(
    IOobject
    (
        "pMean",
        runTime.path()/"POD"/"Mean",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

surfaceScalarField fMean
(
    IOobject
    (
        "fMean",
        runTime.path()/"POD"/"Mean",
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField u0
(
    IOobject
    (
        initialFieldNameU,
        runTime.path()/name(initialConditionTime),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField p0
(
    IOobject
    (
        initialFieldNamep,
        runTime.path()/name(initialConditionTime),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

PtrList<volVectorField> uModes;
uModes = io.readPtrField_noTime<volVectorField>("UMode", runTime, mesh, 1, nModesU, "POD", debug, "UMode");

PtrList<volScalarField> pModes;

if(residual == "galerkin" || residual  == "galerkinPressure")
    {
    pModes = io.readPtrField_noTime<volScalarField>("pMode", runTime, mesh, 1, nModesp, "POD", debug, "pMode");
    }
else if(residual == "PPE")
    {
        pModes = io.readPtrField_noTime<volScalarField>("pMode2", runTime, mesh, 1, nModesp, "POD", debug, "pMode2");
    }

volVectorField uRom
(
    IOobject
    (
        romNameU,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::AUTO_WRITE
        IOobject::NO_WRITE
    ),
    u0
);

volScalarField pRom
(
    IOobject
    (
        romNameP,
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        //IOobject::AUTO_WRITE
        IOobject::NO_WRITE
    ),
    p0
);

// Check how well pressure modes capture gradients
for (label i = 0; i < min(5, pModes.size()); ++i) {
    volVectorField gradMode = fvc::grad(pModes[i]);
    scalar maxGrad = max(mag(gradMode)).value();
    scalar meanGrad = average(mag(gradMode)).value();
    Info << "Mode " << i << ": max|∇p| = " << maxGrad 
         << ", mean|∇p| = " << meanGrad << nl;
}

// Check pressure mode energy distribution
scalar totalEnergy = 0;
for (label i = 0; i < pModes.size(); ++i) {
    scalar energy = fvc::domainIntegrate(pModes[i] * pModes[i]).value();
    totalEnergy += energy;
    if (i < 5) Info << "Mode " << i << " energy: " << energy << nl;
}


scalar pMeanMode0 = average(pModes[0]).value();
Info << "Mean of pressure mode 0: " << pMeanMode0 << nl;

/************************** */

