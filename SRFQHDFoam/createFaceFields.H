    /*
     *
     * Linear interpolation of fields from volumes to face centers
     *
     */

    // Density
    surfaceScalarField rhof
    (
        "rhof",
        linearInterpolate(rho)
    );

    // Velocity
    surfaceVectorField Urelf
    (
        "Urelf",
        linearInterpolate(Urel)
    );
/*
    // Velocity to QGDCourant.H
    surfaceVectorField Uf
    (
        "Urelf",
        linearInterpolate(Urel)
    );
*/
    //
    surfaceVectorField Wf
    (
        "Wf",
        linearInterpolate(W)
    );

    surfaceVectorField BdFrcf
    (
        "BdFrcf",
        linearInterpolate(BdFrc)
    );

    surfaceScalarField Tf
    (
        "Tf",
        linearInterpolate(T)
    );

    // Heat capacity at constant pressure
    surfaceScalarField Cpf
    (
        "Cpf",
        linearInterpolate(thermo.Cp())
    );

    surfaceScalarField muf
    (
        "muf",
        linearInterpolate(thermo.mu()*0.0)
    );

    /*
     *
     * QGD coefficients
     *
     */
    // OpenFOAM delta coeffs can be negative for coupled patches
    const surfaceScalarField& tauQGDf = thermo.tauQGDf();
    surfaceScalarField taubyrhof
    (
        "taubyrhof",
        tauQGDf/rhof
    );

    surfaceScalarField alphaf
    (
        "alphaf",
        linearInterpolate(thermo.alpha())
    );
