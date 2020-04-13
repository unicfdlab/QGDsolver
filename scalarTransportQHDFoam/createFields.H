Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const surfaceScalarField& hQGDf = thermo.hQGDf();
const surfaceScalarField& tauQGDf = thermo.tauQGDf();

Info << "Thermo corrected" << endl;

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

surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(U)
);

surfaceScalarField  phiwo
(
    "phiwo",
    mesh.Sf() & linearInterpolate(U)
);

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & (linearInterpolate(U))
);

Switch implicitDiffusion (thermo.implicitDiffusion());

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

fvScalarMatrix TSu
(
    T,
    T.dimensions()*dimVolume/dimTime
);

//
//END-OF-FILE
//

