//Gradients and divergence
//---------Start---------
surfaceVectorField gradPf 
(
    "gradPf", faceStencil.faceScalarGrad(p)
);

surfaceTensorField gradUf
(
    "gradUf",
    faceStencil.faceVectorGrad(U)
);

surfaceScalarField divUf
(
    "divUf",
    faceStencil.faceVectorDiv(U)
);
//---------End---------

//Continuity equation fluxes
//---------Start---------
surfaceScalarField phivf
(
    "phivf",
    Uf & mesh.Sf()
);

phi = phivf*rhof;

surfaceVectorField rhoW1
(
    "rhoW1",
    tauQGD * faceStencil.faceTensorDiv(rho * (U * U))
);

surfaceScalarField phiRhoW1
(
    "phiRhoW1", rhoW1 & mesh.Sf()
);

surfaceVectorField rhoW2
(
    "rhoW2",
    tauQGD * gradPf
);

surfaceScalarField phiRhoW2
(
    "phiRhoW2",
    rhoW2 & mesh.Sf()
);

surfaceVectorField jm
(
    "jm",
    Uf*rhof - rhoW1 - rhoW2
);

surfaceScalarField phiJm
(
    "phiJm",
    jm & mesh.Sf()
);
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiJmU
(
    "phiJmU",
    (jm * Uf) & mesh.Sf()
);

surfaceVectorField phiP
(
    "phiP",
    pf*mesh.Sf()
);
        
surfaceTensorField Pif
(
    "Pif",
    //QGD diffusive fluxes
    tauQGD * 
    (
        Uf * (rhof * (Uf & gradUf) + gradPf)
        +
        I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
    )
    +
    //NS diffusive fluxes
    (
        muQGD*(gradUf + gradUf.T())
        -
        muQGD*I*divUf
    )
);

surfaceVectorField phiPi
(
    "phiPi",
    Pif & mesh.Sf()
);
//---------End---------

// Fluxes for energy balance equation
//---------Start---------
surfaceScalarField phiJmH
(
    "phiJmH",
    phiJm * Hf
);

surfaceVectorField qf
(
    "qf",
    -
    kappaQGD*faceStencil.faceScalarGrad(T)
    -
    tauQGD* 
    ( 
        ((linearInterpolate(rho*(U*U))) & faceStencil.faceScalarGrad(e))
        +
        (pf * rhof * Uf * (Uf & faceStencil.faceScalarGrad(1/rho)))
    )
);

surfaceScalarField phiQ
(
    "phiQ",
    qf & mesh.Sf()
);

surfaceScalarField phiPiU
(
    "phiPiU",
    (Pif & Uf) & mesh.Sf()
);
//---------End---------
