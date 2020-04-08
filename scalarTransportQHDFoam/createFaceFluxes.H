
surfaceVectorField gradTf
(
    "gradTf", fvsc::grad(T)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceScalarField phiTf
(
    "phiTf",
    phi*Tf
);

surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);

//---------End---------
