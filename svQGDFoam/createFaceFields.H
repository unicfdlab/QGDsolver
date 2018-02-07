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
    
    surfaceScalarField zetaf
    (
        "zetaf",
        muf*0.0
    );
    
    /*
     *
     * QGD coefficients
     *
     */
    surfaceScalarField tauQGDf
    (
        "tauQGD",
        linearInterpolate(tauQGD)
    );
    
    surfaceScalarField hQGDf
    (
        "hQGDf",
        1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs())
    );
    
    {
        scalar hown = 0.0;
        scalar hnei = 0.0;
        forAll(hQGDf.primitiveField(), iFace)
        {
            hown = mag(mesh.C()[mesh.owner()[iFace]] - mesh.Cf()[iFace]);
            hnei = mag(mesh.C()[mesh.neighbour()[iFace]] - mesh.Cf()[iFace]);
            hQGDf.primitiveFieldRef()[iFace] = 2.0*min(hown, hnei);
        }
        
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (!fvp.coupled())
            {
                hQGDf.boundaryFieldRef()[patchi] *= 2.0;
            }
        }
    }
    
    surfaceScalarField cf
    (
        "cf",
        linearInterpolate(thermo.c())
    );
    