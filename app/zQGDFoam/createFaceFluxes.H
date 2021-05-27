
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

surfaceScalarField divRhoUf
(
    "divRhoUf",
    fvsc::div(rhoU)
);

//---------End---------

//Continuity equation fluxes
//---------Start---------

phi = linearInterpolate(rhoU) & mesh.Sf();

surfaceVectorField jm
(
    "jm",
    //rhoUf - rhoW
    linearInterpolate(rhoU)
);

surfaceScalarField phiJm
(
    "phiJm",
    mesh.Sf() & jm
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
    pf*mesh.Sf()
);

surfaceTensorField Pif
(
    "Pif",
    //QGD diffusive fluxes
    tauQGDf * (Uf * gradPf)
    +
    tauQGDf * I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
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
    tauMCPtr() & mesh.Sf()
);

//---------End---------

// Fluxes for energy balance equation
//---------Start---------
surfaceScalarField phiJmH
(
    "phiJmH",
    phiJm * linearInterpolate(e)
);

surfaceScalarField phih2
(
    "phih2",
    (0.25*hQGDf*hQGDf*(fvc::snGrad(U) * fvc::snGrad(p)) & mesh.Sf())
);
phih2.setOriented(false);

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
