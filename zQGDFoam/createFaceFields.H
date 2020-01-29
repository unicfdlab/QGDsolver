    /*
     *
     * Linear interpolation of fields from volumes to face centers
     *
     */

    surfaceScalarField pos
    (
        IOobject
        (
            "pos",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("pos", dimless, 1.0)
    );
    
    surfaceScalarField neg
    (
        IOobject
        (
            "neg",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("neg", dimless, -1.0)
    );

    // Density
    surfaceScalarField rhof
    (
        "rhof",
        linearInterpolate(rho)
    );
    surfaceScalarField rho_pos
    (
        interpolate(rho,pos)
    );
    surfaceScalarField rho_neg
    (
        interpolate(rho,neg)
    );
    
    surfaceScalarField rhoLnf = 1.0 / logMean(rho_pos,rho_neg);
    
    surfaceScalarField ef
    (
        "ef",
        linearInterpolate(e)
    );
    surfaceScalarField e_pos
    (
        interpolate(e,pos)
    );
    surfaceScalarField e_neg
    (
        interpolate(e,neg)
    );
    surfaceScalarField eLnf= e_pos*e_neg*logMean(e_pos,e_neg);

    // Velocity
    surfaceVectorField Uf
    (
        "Uf",
        linearInterpolate(U)
    );
    surfaceVectorField U_pos
    (
        interpolate(U,pos)
    );
    surfaceVectorField U_neg
    (
        interpolate(U,neg)
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
    surfaceScalarField E1f
    (
        "E1f",
        linearInterpolate(rhoE + p)
        //or (rhoEf + pf)/rhof?
    );
    
    surfaceScalarField p1f
    (
        "p1f",
        gammam1f*rhof*ef
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


        