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


IOdictionary errorDict
(
    IOobject
    (
        "errorDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<word> ROMfieldU;
List<word> ROMfieldP;

word FOMfieldU = errorDict.lookupOrDefault<word>("FOMfieldU","U");
ROMfieldU = errorDict.lookup("ROMfieldU");

word FOMfieldP = errorDict.lookupOrDefault<word>("FOMfieldP","p");
ROMfieldP = errorDict.lookup("ROMfieldP");

std::vector<word> namesU(ROMfieldU.size());
forAll(ROMfieldU, i)
{
    namesU[i] = ROMfieldU[i];
}

std::vector<word> namesP(ROMfieldP.size());
forAll(ROMfieldP, i)
{
    namesP[i] = ROMfieldP[i];
}
/*******************************/
