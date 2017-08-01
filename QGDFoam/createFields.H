Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();
thermo.correct();

volScalarField& e = thermo.he();

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.primitiveField()) > 0.0)
{
    inviscid = false;
}

#warning "needs to be moved in dynamic library"
extendedFaceStencil faceStencil
(
    IOobject
    (
        "testFaceStencil",
        runTime.timeName(),
        mesh,
        regIOobject::NO_READ,
        regIOobject::NO_WRITE
    ),
    false
);
    
IOdictionary QGDDict
(
    IOobject
    (
        "QGDDict",              // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);
    
scalar ScQGD(readScalar(QGDDict.lookup("ScQGD")));
scalar PrQGD(readScalar(QGDDict.lookup("PrQGD")));

surfaceScalarField alphaQGD
(
    IOobject
    (
        "alphaQGD",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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
    rho*(e + 0.5*magSqr(U))
);

volScalarField rPsi
(
    "rPsi", 
    1.0/psi
);

volScalarField gamma
(
    "gamma",
    thermo.Cp()/thermo.Cv()
);

volScalarField c
(
    "c",
    sqrt(gamma * rPsi)
);

surfaceScalarField phi
(
    "phi",
    linearInterpolate(rho*U) & mesh.Sf()
);

