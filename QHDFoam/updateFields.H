/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
muf        = linearInterpolate(turbulence->muEff());
alphaf = linearInterpolate(turbulence->alphaEff());

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

//Hi
Hif = alphaf/rhof;

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());




