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
    
    // Pressure
    surfaceScalarField pf
    (
        "pf",
        linearInterpolate(p)
    );
    
    // Heat capacities ratio
    surfaceScalarField gammaf 
    (
        "gammaf",
        linearInterpolate(gamma)
    );
    
    surfaceScalarField gammam1f
    (
        "gammam1",
        gammaf - 1
    );
    
    // Speed of sound
    surfaceScalarField cf
    (
        "cf",linearInterpolate(c)
    );

    // Heat capacity at constant pressure
    surfaceScalarField Cpf
    (
        "Cpf",
        linearInterpolate(thermo.Cp())
    );
    
    surfaceScalarField Hf
    (
        "Hf",
        linearInterpolate((rhoE + p)/rho)
        //or (rhoEf + pf)/rhof?
    );
    
    /*
     *
     * QGD coefficients
     *
     */
    // OpenFOAM delta coeffs can be negative for coupled patches
    surfaceScalarField hQGD
    (
        "hQGD",
        1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs())
    );
    
    surfaceScalarField tauQGD
    (
        "tauQGD",
        alphaQGD * hQGD / cf
    );
    
    surfaceScalarField muQGD
    (
        "muQGD",
        tauQGD*pf*ScQGD
    );
    
    surfaceScalarField kappaQGD
    (
        "kappaQGDf",
        muQGD*Cpf / PrQGD
    );
    