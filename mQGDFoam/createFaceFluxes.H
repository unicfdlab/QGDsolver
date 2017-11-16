
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

//surfaceTensorField gradUfT
//(
//    "gradUfT",
//    Foam::T(gradUf)
//);

surfaceScalarField divUf
(
    "divUf",
    Foam::tr(gradUf)
);

surfaceVectorField gradef
(
    "gradef",
    fvsc::grad(e)
);

surfaceVectorField gradRhof
(
    "gradRhof",
    fvsc::grad(rho)
);

//---------End---------

//Continuity equation fluxes
//---------Start---------

surfaceVectorField rhofUf
(
    "rhofUf",
    rhof*Uf
);

phi = mesh.Sf() & rhofUf;

surfaceVectorField rhoW
(
    "rhoW",
    tauQGDf * 
    (
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
    mesh.Sf() & (Uf*rhof - rhoW)
);

//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiJmU
(
    "phiJmU",
    mesh.Sf() & (jm * Uf)
);

surfaceVectorField phiP
(
    "phiP",
    mesh.Sf() * pf
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
    mesh.Sf() & Pif
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
    mesh.Sf() & tauMCPtr()
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
        -
        (pf * rhof * Uf * (Uf & gradRhof) / rhof / rhof)
    )
);

surfaceScalarField phiQ
(
    "phiQ",
    mesh.Sf() & qf
);

surfaceScalarField phiPiU
(
    "phiPiU",
    mesh.Sf() & (Pif & Uf)
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
    mesh.Sf() & sigmaDotUPtr()
);

//---------End---------
