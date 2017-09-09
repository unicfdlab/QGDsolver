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
const volScalarField& hQGD = thermo.hQGD();
const volScalarField& tauQGD = thermo.tauQGD();
const volScalarField& gamma = thermo.gamma();
const volScalarField& c = thermo.c();
const volScalarField& psi = thermo.psi();

Info << "Thermo corrected" << endl;

autoPtr<compressible::turbulenceModel> turbulence;
autoPtr<volScalarField> muEffPtr;

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

volTensorField rhoUU
(
    IOobject
    (
        "rhoUU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U*U
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

surfaceScalarField phi
(
    "phi",
    linearInterpolate(rho*U) & mesh.Sf()
);

bool inviscid(true);
if (max(thermo.mu()().primitiveField()) > SMALL)
{
    inviscid = false;
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

    muEffPtr.reset
    (
        new volScalarField
        (
            "muEff",
            turbulence->muEff()
        )
    );
}
else
{
    muEffPtr.reset
    (
        new volScalarField
        (
            "muEff",
            thermo.mu()*0.0
        )
    );
}
//
//END-OF-FILE
//

