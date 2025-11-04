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

