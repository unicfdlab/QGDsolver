//Gradients and divergence
//---------Start---------

surfaceVectorField gradPf
(

    "gradPf", fvsc::grad(p)
);
    Info<< "grP"<< nl << endl;
//surfaceVectorField gradPf
//(
//    "gradP_rghf", fvsc::grad(p)
//);

surfaceTensorField gradUrelf
(
    "gradUrelf",
    fvsc::grad(Urel)
);

//---------End---------
//P equation fluxes
//---------Start---------
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiUrelf
(
    "phiUrelf",
    phi * Urelf + (mesh.Sf() & (Wf * Wf))
);
/*
surfaceScalarField phiTf
(
    "phiTf",
    phi*Tf
);
*/
//surfaceVectorField phiPi
//(
//    "phiPi",
//    muf/rhof * mesh.Sf() & (gradUf+gradUf.T())
//);

//---------End---------

// Fluxes for T balance equation
//---------Start---------
surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);
Info << "Creating Face Fluxes is done"<< endl;

//---------End---------
