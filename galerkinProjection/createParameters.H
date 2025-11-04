/*********************************************/
/*********************************************/
/************Creates helper dicts************/
/*********************************************/
/*********************************************/

IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

scalar startTime = readScalar(controlDict.lookup("startTime"));
scalar endTime = readScalar(controlDict.lookup("endTime"));
scalar dt = readScalar(controlDict.lookup("deltaT"));
scalar writeInterval = readLabel(controlDict.lookup("writeInterval"));

IOdictionary podDict
(
    IOobject
    (
        "podDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

label nModesU = podDict.lookupOrDefault<label>("nModesU", 5);
label nModesp = podDict.lookupOrDefault<label>("nModesp", 5);
word romNameU = podDict.lookupOrDefault<word>("romNameU", "URom");
word romNameP = podDict.lookupOrDefault<word>("romNameP", "pRom");
word residual = podDict.lookupOrDefault<word>("stabilization", "galerkin");

IOdictionary debugDict
(
    IOobject
    (
        "debugDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const Switch debug = debugDict.lookupOrDefault<Switch>("debug",false);


IOdictionary parameterDict
(
    IOobject
    (
        "parameterDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
word patchName = parameterDict.lookupOrDefault<word>("patch","inlet");

IOdictionary newtonSolverDict
(
    IOobject
    (
        "newtonSolverDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
label maxIter = readLabel(newtonSolverDict.lookup("maxIter"));
scalar tol = readScalar(newtonSolverDict.lookup("tol"));
const Switch printInfos = newtonSolverDict.lookupOrDefault<Switch>("printInfo",false);
scalar initialConditionTime = newtonSolverDict.lookupOrDefault<scalar>("initialConditionTime", 0);
word initialFieldNameU = newtonSolverDict.lookupOrDefault<word>("fieldU", "U");
word initialFieldNamep = newtonSolverDict.lookupOrDefault<word>("fieldP", "p");

