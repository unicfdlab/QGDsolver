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
    surfaceVectorField Uf
    (
        "Uf",
        linearInterpolate(U)
    );

    surfaceScalarField Tf
    (
        "Tf",
        linearInterpolate(T)
    );
    
    /*
     *
     * QGD coefficients
     *
     */
    // OpenFOAM delta coeffs can be negative for coupled patches
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


