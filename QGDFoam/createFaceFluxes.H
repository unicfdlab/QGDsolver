
/*Field for debugging parallel execution*/
//surfaceVectorField cellNoGradf
//(
//    "cellNoGradf", faceStencil.faceScalarGrad(cellNo)
//);
//
//cellNoGradf.write();

//Gradients and divergence
//---------Start---------
surfaceVectorField gradPf 
(
    "gradPf", fvsc::grad(p)
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
surfaceVectorField phiUUmWUmWW
(
    "phiUUmWUmWW",
    ((Uf * Uf)-(Wf * Uf)-(Wf * Wf)) & mesh.Sf()
);

surfaceScalarField Pbyrho
(
    "Pbyrho",
    pf/rhof
);

surfaceVectorField phiPi
(
    "phiPi",
    muf*((gradUf+gradUf.T())) & mesh.Sf()/rhof
);

surfaceVectorField phiBdFrc
(
    "phiBdFrc",
    betaf*g*Tf
);
//---------End---------

// Fluxes for T balance equation
//---------Start---------
surfaceScalarField phiUmWT
(
    "phiUmWT",
    ((Uf-Wf)*Tf) & mesh.Sf()
);

surfaceScalarField HiT
(
    "HiT",
    alphaf/rhof*Tf
);
//---------End---------
