Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiQGDThermo> pThermo
(
    psiQGDThermo::New(mesh)
);
psiQGDThermo& thermo = pThermo();
thermo.correct();

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();

const surfaceScalarField& tauQGDf = thermo.tauQGDf();
const surfaceScalarField& hQGDf   = thermo.hQGDf();
const volScalarField& gamma = thermo.gamma();
const volScalarField& psi = thermo.psi();

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
    thermo.rho()
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*e + rho*0.5*magSqr(U)
);

volScalarField H
(
    IOobject
    (
        "Ht",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (rhoE + p)/rho
);

surfaceScalarField phi
(
    "phi",
    linearInterpolate(rho*U) & mesh.Sf()
);

surfaceScalarField phiRhof
(
    "phiRhof",
    linearInterpolate(rho)*phi
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

#include "createZeroSources.H"

//
//END-OF-FILE
//

