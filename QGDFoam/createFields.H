Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoQGDThermo> pThermo
(
    rhoQGDThermo::New(mesh)
);
rhoQGDThermo& thermo = pThermo();
thermo.correct();

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& hQGD = thermo.hQGD();
const volScalarField& tauQGD = thermo.tauQGD();
const volScalarField& beta = thermo.beta();
//const volScalarField& gamma = thermo.gamma();
//const volScalarField& c = thermo.c();
//const volScalarField& psi = thermo.psi();

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

volVectorField W
(
    IOobject
    (
        "W",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    )
    U
)

surfaceScalarField phiu
(
    "phiu",
    linearInterpolate(U) & mesh.Sf()
);

surfaceScalarField  phiwo
(
    "phiwo"
    linearInterpolate(W) & mesh.Sf()
);

Switch implicitDiffusion (thermo.implicitDiffusion());

Info<< "Creating turbulence model\n" << endl;
turbulence.reset
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    ).ptr()
);

//
//END-OF-FILE
//

