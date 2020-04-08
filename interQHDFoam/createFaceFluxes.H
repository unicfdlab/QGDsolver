surfaceVectorField gradAlpha1f
(
    "gradAlpha1f", fvsc::grad(alpha1)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceVectorField gradpf
(
    "gradpf",
    fvsc::grad(p)
);

surfaceScalarField rhoPhi
(
    "rhoPhi",
    rho1*phi
);

surfaceVectorField phiUfRhof
(
    "phiUf",
    phi * Uf * rho1
); 
phiUfRhof.setOriented(true);

surfaceVectorField phiRhofWf
(
    "phiUf",
    phi * Uf * rho1
);
phiRhofWf.setOriented(true);

surfaceScalarField phiAlpha1f
(
    "phiAlpha1f",
    phi*alpha1f
);

surfaceScalarField phiAlpha2f
(
    "phiAlpha2f",
    phi*alpha1f
);

surfaceScalarField phiwo1
(
    "phiwo1",
    phi*alpha1f
);
surfaceScalarField phiwo2
(
    "phiwo2",
    phi*alpha1f
);
surfaceScalarField phiw1
(
    "phiw1",
    phi*alpha1f
);
surfaceScalarField phiw2
(
    "phiw2",
    phi*alpha1f
);

surfaceScalarField phiwm
(
    "phiwm",
    phiwo1*0.0
);
surfaceScalarField coeffp
(
    "coeffp",
    linearInterpolate(Tau1 / rho)
);

surfaceScalarField phiWr
(
    "phiWr",
    phiwo1*0.0
);

surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(U)
);
surfaceScalarField phi1
(
    "phi1",
    mesh.Sf() & linearInterpolate(U)
);
surfaceScalarField phi2
(
    "phi2",
    mesh.Sf() & linearInterpolate(U)
);


//
//END-OF-FILE
//


