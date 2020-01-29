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

    // Liquid fraction
    surfaceScalarField alphalf
    (
        "rhof",
        linearInterpolate(alphal)
    );
    
    // Velocity
    surfaceVectorField Uf
    (
        "Uf",
        linearInterpolate(U)
    );

    surfaceScalarField pf
    (
        "pf",
	linearInterpolate(p)
    );
    
    //
    surfaceVectorField rhoWf
    (
        "rhoWf",
        linearInterpolate(rho*W)
    );
    
    surfaceVectorField rhoUf
    (
        "rhoUf",
        linearInterpolate(rho*U)
    );
    
    surfaceVectorField BdFrcf
    (
        "BdFrcf",
        linearInterpolate(BdFrc)
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


