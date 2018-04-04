/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);

// Velocity
Uf = linearInterpolate(U);
    
// Pressure
pf = linearInterpolate(p);

// T
Tf = linearInterpolate(T);

alphaf = linearInterpolate(thermo.alpha());

//Hi
Hif = alphaf/rhof;

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
tauQGDf    = linearInterpolate(tauQGD);
muf        = linearInterpolate(turbulence->muEff());


