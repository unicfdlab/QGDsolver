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
    
    surfaceScalarField muf
    (
        "muf",
        linearInterpolate(thermo.mu()*0.0)
    );
    if (!inviscid)
    {
        muf = linearInterpolate(turbulence->muEff());
    }

    /*
     *
     * QGD coefficients
     *
     */
    // OpenFOAM delta coeffs can be negative for coupled patches
    surfaceScalarField hQGDf
    (
        "hQGD",
        1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs())
    );
    forAll(mesh.boundaryMesh(), iPatch)
    {
        if (!isA<coupledFvPatch>(mesh.boundaryMesh()[iPatch]))
        {
            hQGDf.boundaryFieldRef()[iPatch] *= 2.0;
        }
    }
    
    surfaceScalarField tauQGDf
    (
        "tauQGD",
        linearInterpolate(tauQGD)
    );
    