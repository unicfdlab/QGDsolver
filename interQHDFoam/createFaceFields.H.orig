    /*
     *
     * Linear interpolation of fields from volumes to face centers
     *
     */
    
    // Density
    
    // Liquid fraction
    surfaceScalarField alpha1f
    (
        "alpha1f",
        linearInterpolate(alpha1)
    );
    surfaceScalarField alpha2f
    (
        "alpha2f",
        1.0 - alpha1f
    );
    
    // Velocity
    surfaceVectorField Uf
    (
        "Uf",
        linearInterpolate(U)
    );
    surfaceVectorField W1
    (
        "Wf",
        linearInterpolate(U)*0.0
    );
    surfaceVectorField W2
    (
        "Wf",
        linearInterpolate(U)*0.0
    );

    // Forces
    surfaceVectorField BdFrcf
    (
        "BdFrcf",
        linearInterpolate(BdFrc)
    );
    
    surfaceVectorField cFrcf
    (
        "cFrcf",
        linearInterpolate(cFrc)
    );
    
    // fluid viscosity
    surfaceScalarField muf
    (
        "muf",
        linearInterpolate(thermo.mu()*0.0)
    );
    




