#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"



Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
volScalarField& T = const_cast<volScalarField&>(thermo.T());
const surfaceScalarField& hQGDf = thermo.hQGDf();
const surfaceScalarField& tauQGDf = thermo.tauQGDf();

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

volScalarField T0
(
    "T0",
    T
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
    thermo.rho()
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
    "phiwStar",
    mesh.Sf() & linearInterpolate(W)
);

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & (linearInterpolate(U) - linearInterpolate(W))
);

surfaceScalarField phiRhof
(
    "phiRhof",
    linearInterpolate(rho)*phi
);


volVectorField BdFrc
(
    "BdFrc",
    T*g*beta
);

Switch implicitDiffusion (thermo.implicitDiffusion());

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
rho.oldTime();

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

#include "createIcoZeroSources.H"

//
//END-OF-FILE
//

