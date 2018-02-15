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
    
    //
    surfaceVectorField Wf
    (
        "Wf",
        linearInterpolate(W)
    );

    // Pressure
    surfaceScalarField pf
    (
        "pf",
        linearInterpolate(p)
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
    
    surfaceScalarField taubyrhof
    (
        "taubyrhof",
	linearInterpolate(tauQGD/rho)
    );


    surfaceScalarField tauQGDf
    (
        "tauQGD",
        linearInterpolate(tauQGD)
    );
    

    surfaceScalarField betaf
    (
        "betaf",
	linearInterpolate(beta)
    );

    surfaceScalarField alphaf
    (
        "alphaf",
	linearInterpolate(thermo.alpha())
    );

    surfaceVectorField BdFrcf
    (
        "BdFrcf",
	linearInterpolate(BdFrc)
    );
