//Gradients and divergence
//---------Start---------

surfaceVectorField gradPf 
(
    "gradPf", fvsc::grad(p)
);

//surfaceVectorField gradPf
//(
//    "gradP_rghf", fvsc::grad(p)
//);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceTensorField gradWf
(
    "gradWf",
    fvsc::grad(W)
);

surfaceVectorField gradTf
(
    "gradTf",
    fvsc::grad(T)
);

//---------End---------
//P equation fluxes
//---------Start---------
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiUf
(
    "phiUf",
    phi * Uf + (mesh.Sf() & (Wf * Wf))
);

surfaceScalarField phiTf
(
    "phiTf",
    phi*Tf
);

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
