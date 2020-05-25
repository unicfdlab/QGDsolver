
surfaceVectorField gradTf
(
    "gradTf", fvsc::grad(T)
);

surfaceScalarField phiTf
(
    "phiTf",
    phiu*Tf
);

surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);

//---------End---------
