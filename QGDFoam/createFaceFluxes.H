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

surfaceVectorField rhofUf
(
    "rhofUf",
    rhof*Uf
);

surfaceVectorField rhoW
(
    "rhoW",
    tauQGDf * 
    (
        faceStencil.faceTensorDiv(rho * (U * U))
        +
        gradPf
    )
);

surfaceVectorField jm
(
    "jm",
    rhofUf - rhoW
);

surfaceScalarField phiJm
(
    "phiJm",
    (Uf*rhof - rhoW) & mesh.Sf()
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
    tauQGDf * 
    (
        Uf * (rhof * (Uf & gradUf) + gradPf)
        +
        I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
    )
);

surfaceVectorField phiPi
(
    "phiPi",
    Pif & mesh.Sf()
);

autoPtr<surfaceTensorField> tauMCPtr;
tauMCPtr.reset
(
    new surfaceTensorField
    (
        "tauMC",
        //muEffPtr()*dev2(Foam::T(fvc::grad(U)))
        muf*
        (
            Foam::T(faceStencil.faceVectorGrad(U))//Don't forget to transpose
            -
            I*faceStencil.faceVectorDiv(U)
        )
    )
);

surfaceVectorField phiTauMC
(
    "phiTauMC",
    tauMCPtr() & mesh.Sf()
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
    tauQGDf* 
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

autoPtr<surfaceVectorField> sigmaDotUPtr;
sigmaDotUPtr.reset
(
    new surfaceVectorField
    (
        "sigmaDotU",
        tauMCPtr() & linearInterpolate(U) //to be updated later
    )
);

surfaceScalarField phiSigmaDotU
(
    "phiSigmaDotU",
    sigmaDotUPtr() & mesh.Sf()
);

//---------End---------
