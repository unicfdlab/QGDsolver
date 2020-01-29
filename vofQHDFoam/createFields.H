Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();

const dimensionedScalar rhol (thermo.lookup("rhol"));
const dimensionedScalar rhog (thermo.lookup("rhog"));
const dimensionedScalar mul  (thermo.lookup("mul"));
const dimensionedScalar mug  (thermo.lookup("mug"));

volScalarField alphal
(
    IOobject
    (
        "alphal",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
const surfaceScalarField& tauQGDf = thermo.tauQGDf();
const surfaceScalarField& hQGDf = thermo.hQGDf();

Info << "Thermo corrected" << endl;

autoPtr<compressible::turbulenceModel> turbulence;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    (rhol - rhog)*alphal + rhog
);

volVectorField W
(
    IOobject
    (
        "W",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    U
);

volVectorField rhoU
(
    "rhoU",
    rho*U
);

Info<< "\nReading gravitationalProperties" << endl;
  
IOdictionary gravitationalProperties
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g"));

dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    thermo.subDict("mixture").subDict("transport")
);


surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(U)
);

surfaceScalarField  phiwo
(
    "phiwo",
    mesh.Sf() & linearInterpolate(W)
);

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & (linearInterpolate(U) - linearInterpolate(W))
);

volVectorField BdFrc
(
    "BdFrc",
    rho*g
);

Switch implicitDiffusion (thermo.implicitDiffusion());

surfaceScalarField phiRhof
(
    "phiRhof",
    linearInterpolate(rho)*phi
);

Info<< "Creating turbulence model\n" << endl;
turbulence.reset
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phiRhof,
        thermo
    ).ptr()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

//
//END-OF-FILE
//

