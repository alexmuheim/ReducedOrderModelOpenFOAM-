/*********************************************/
/*********************************************/
/************Creates helper fields************/
/*********************************************/
/*********************************************/

volVectorField uFom
(
    IOobject
    (
        FOMfieldU,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField pFom
(
    IOobject
    (
        FOMfieldP,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField uRom
(
    IOobject
    (
        ROMfieldU[i],
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField pRom
(
    IOobject
    (
        ROMfieldP[i],
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


//diff fields
volVectorField uDiff
(
    IOobject
    (
        "uDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (uFom - uRom)
);

volScalarField pDiff
(
    IOobject
    (
        "pDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (pFom - pRom)
);

volVectorField pGrad
(
    IOobject
    (
        "pGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(pFom)
);

volVectorField pDiffGrad
(
    IOobject
    (
        "pDiffGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(pDiff)
);






