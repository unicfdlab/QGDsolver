
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

//---------End---------

//Continuity equation fluxes
//---------Start---------

surfaceVectorField rhoUf
(
    "rhoUf",
    linearInterpolate(rho*U)
);

surfaceVectorField rhofUf
(
    "rhofUf",
    rhof*Uf
);

phi = rhofUf & mesh.Sf();

surfaceVectorField rhoW
(
    "rhoW",
    tauQGDf * ((rhofUf & gradUf) + gradPf)
);

surfaceVectorField jm
(
    "jm",
    rhofUf - rhoW
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
    Uf * rhoW
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
    phiJm * Hf
);

surfaceVectorField qf
(
    "qf",
    - alphauf * gradef * 0.0 //make flux zero for implicit diffusion
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
