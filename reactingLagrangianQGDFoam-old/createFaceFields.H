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

    //
    surfaceScalarField yf
    (
        "yf",
        linearInterpolate(y)
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
    
    surfaceScalarField muf
    (
        "muf",
        linearInterpolate(thermo.mu()*0.0)
    );
    
    surfaceScalarField alphauf
    (
        "alphauf",
        muf*0.0
    );
    
    surfaceScalarField cf
    (
        "cf",
        linearInterpolate(thermo.c())
    );
    