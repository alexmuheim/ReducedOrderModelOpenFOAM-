
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
scalar startSaveSnap = podDict.lookupOrDefault<scalar>("startTime", startTime);
scalar endSaveSnap = podDict.lookupOrDefault<scalar>("endTime", endTime);
label nModesU = podDict.lookupOrDefault<label>("nModesU", 5);
label nModesp = podDict.lookupOrDefault<label>("nModesp", 5);
word residual = podDict.lookupOrDefault<word>("stabilization", "galerkin");

if(endSaveSnap > endTime)
{
    FatalError << "snapshot saving endTime can't be larger than simulation endTime " 
               << exit(FatalError); 
}
if(startSaveSnap >= endSaveSnap)
{
    Warning << "starTime can't be greater or equal endTime. " << nl
            << "changing start and end time to default values... startTime = " 
            << startTime << ", endTime = " << endTime << endl; 
            startSaveSnap = startTime;
            endSaveSnap = endTime;
}


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

