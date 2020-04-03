//#include "readGravitationalAcceleration.H"
#include "readhRef.H"
//#include "gh.H"


Info<< "Reading thermophysical properties\n" << endl;
autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();


Info<< "Reading field Urel\n" << endl;
volVectorField Urel
(
    IOobject
    (
        "Urel",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating SRF model\n" << endl;
autoPtr<SRF::SRFModel> SRF
(
    SRF::SRFModel::New(Urel)
);


//volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
volScalarField& T = const_cast<volScalarField&>(thermo.T());
//const surfaceScalarField& hQGDf = thermo.hQGDf();
//const surfaceScalarField& tauQGDf = thermo.tauQGDf();

Info << "Thermo corrected" << endl;

autoPtr<compressible::turbulenceModel> turbulence;
/*
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
*/

/*
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    thermo.subDict("mixture").subDict("transport")
);
*/
//volScalarField T0
//(
//    "T0",
//    T
//);
//T0.correctBoundaryConditions();

/*
volScalarField rhok0
(
    "rhok0",
    1.0 - beta*T
);

volScalarField rhok
(
    "rhok",
    1.0 - beta*T
);
*/

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
    Urel
);

surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(Urel)
);

surfaceScalarField  phiwo
(
    "phiwStar",
    mesh.Sf() & linearInterpolate(W)
);

//volScalarField p_rgh
//(
//    IOobject
//    (
//        "p_rgh",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & (linearInterpolate(Urel) - linearInterpolate(W))
);

surfaceScalarField phiRhof
(
    "phiRhof",
    linearInterpolate(rho)*phi
);


volVectorField BdFrc
(
    "BdFrc",
//    (rhok-rhok0)*g
    -2.0*(SRF->omega()^Urel)
);

Switch implicitDiffusion (thermo.implicitDiffusion());

Info<< "Creating turbulence model\n" << endl;
turbulence.reset
(
    compressible::turbulenceModel::New
    (
        rho,
        Urel,
        phiRhof,
        thermo
    ).ptr()
);
rho.oldTime();

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);

volScalarField alpha1(T);

// MULES Correction
tmp<surfaceScalarField> talphaPhiCorr0;


//
//END-OF-FILE
//

