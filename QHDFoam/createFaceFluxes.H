
/*Field for debugging parallel execution*/
//cellNo.correctBoundaryConditions();
//surfaceVectorField cellNoGradf
//(
//    "cellNoGradf", fvsc::grad(cellNo)
//);
//
//cellNoGradf.write();

//Gradients and divergence
//---------Start---------


surfaceVectorField gradPf 
(
    "gradPf", fvsc::grad(p)
);

surfaceVectorField gradTf
(
    "gradTf", fvsc::grad(T)
);

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

surfaceVectorField gradRhof
(
    "gradRhof",
    fvsc::grad(rho)
);

surfaceScalarField divUf
(
    "divUf",
    tr(gradUf)
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

surfaceScalarField Pbyrho
(
    "Pbyrho",
    pf/rhof
);

surfaceScalarField phiTf
(
    "phiTf",
    phi*Tf
);

surfaceVectorField phiPi
(
    "phiPi",
    muf/rhof * mesh.Sf() & (gradUf+gradUf.T())
);

//---------End---------

// Fluxes for T balance equation
//---------Start---------
surfaceScalarField phiUmWT
(
    "phiUmWT",
    ((Uf-Wf)*Tf) & mesh.Sf()
);

surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);
Info << "Creating Face Fluxes is done"<< endl;

//---------End---------
