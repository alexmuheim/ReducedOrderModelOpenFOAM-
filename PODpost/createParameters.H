
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
scalar thresholdEnergy = podDict.lookupOrDefault<scalar>("threshold", 0.999);
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
const label infoLevel = debugDict.lookupOrDefault<label>("infoLevel",1);

IOdictionary parameterDict
(
    IOobject
    (
        "parameterDict",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<scalar> nuList;
List<vector> dirichletList;
word parameterSpace(parameterDict.lookupOrDefault<word>("parameter", "viscosity"));
word parPatch(parameterDict.lookupOrDefault<word>("patch","inlet"));
nuList = parameterDict.lookup("nu");
dirichletList = parameterDict.lookup("Ud"); 
label parSize = 0;
label nuSize = nuList.size();
utils.debugMsg("nuSize is " + std::to_string(nuSize), debug);
label dirichletSize = dirichletList.size();
if(parameterSpace == "dirichlet")
{
    parSize = dirichletList.size();
    if(nuSize >= 1)
    {
        Warning << "when dirichlet is chosen, nu should only have 1 value ..." << endl;
    }
}
else if(parameterSpace == "viscosity")
{
    parSize = nuList.size();
    if(dirichletSize != 1)
    {
        Warning << "when viscosity is chosen, Ud should only have 1 vector value ..." << endl;
    }
}
else
{
    FatalError << "parameter '" << parameterSpace << "' is not avaialable!" << nl
                << " you can choose between 'viscosity' or 'dirichlet' "  
                << exit(FatalError);
}

utils.debugMsg("parSize is: " + std::to_string(parSize),debug);
utils.debugMsg(debug, "start saving time: ", startSaveSnap, " end: ", endSaveSnap);