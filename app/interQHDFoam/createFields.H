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

Info<< "Reading thermophysical properties\n" << endl;
twoPhaseIcoQGDThermo thermo(mesh,U);
thermo.correct();


const dimensionedScalar& rho1 = thermo.rho1();
const dimensionedScalar& rho2 = thermo.rho2();
const dimensionedScalar& nu1 = thermo.nu1();
const dimensionedScalar& nu2 = thermo.nu2();
const dimensionedScalar& Tau1 = thermo.Tau1();
const dimensionedScalar& Tau2 = thermo.Tau2();

volScalarField& alpha1 = thermo.alpha1();
volScalarField& alpha2 = thermo.alpha2();
volScalarField& p = thermo.p();
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

surfaceScalarField phi
(
    "phi",
    mesh.Sf() & linearInterpolate(U)
);

volVectorField BdFrc
(
    "BdFrc",
    rho*g
);

volVectorField cFrc
(
    "cFrc",
    fvc::reconstruct(thermo.surfaceTensionForce()*mesh.magSf())
);

surfaceScalarField hQGDf = thermo.hQGDf();

Switch implicitDiffusion (thermo.implicitDiffusion());

Info<< "Creating turbulence model\n" << endl;
// Construct incompressible turbulence model
//autoPtr<incompressible::turbulenceModel> turbulence
//(
//    incompressible::turbulenceModel::New(U, phi, mixture)
//);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, thermo.subDict("QGD"), pRefCell, pRefValue);

//
//END-OF-FILE
//

