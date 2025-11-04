/*********************************************/
/*********************************************/
/************Creates helper fields************/
/*********************************************/
/*********************************************/

volVectorField u0Fom
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

volScalarField p0Fom
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

volVectorField u0Rom
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

volScalarField p0Rom
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
volVectorField u0Diff
(
    IOobject
    (
        "uDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (u0Fom - u0Rom)
);

volScalarField p0Diff
(
    IOobject
    (
        "pDiff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (p0Fom - p0Rom)
);

volVectorField p0Grad
(
    IOobject
    (
        "pGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(p0Fom)
);

volVectorField p0DiffGrad
(
    IOobject
    (
        "pDiffGrad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(p0Diff)
);





