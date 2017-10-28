
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

surfaceScalarField divUf
(
    "divUf",
    fvsc::div(U)
);

surfaceVectorField gradef
(
    "gradef",
    fvsc::grad(e)
);

surfaceVectorField gradInvRhof
(
    "gradInvRhof",
    fvsc::grad(1.0/rho)
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
        fvsc::div(rho * (U * U))
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
        0.0*muf*
        (
            Foam::T(gradUf)//Don't forget to transpose
            -
            I*divUf
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
        ((linearInterpolate(rho*(U*U))) & gradef)
        +
        (pf * rhof * Uf * (Uf & gradInvRhof))
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
