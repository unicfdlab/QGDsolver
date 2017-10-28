/*
 *
 * Viscous part of N-S equations
 *
 */

rhoUU = rhoU*U;

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);

// Momentum
rhofUf = linearInterpolate(rhoU);

// Velocity
Uf = linearInterpolate(U);
    
// Pressure
pf = linearInterpolate(p);

// Heat capacities ratio
gammaf = linearInterpolate(gamma);
gammam1f = gammaf - 1;

// Speed of sound
cf = linearInterpolate(c);

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());
    
// Total specific enthalpy
Hf = linearInterpolate(rhoE + p)/rhof;

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
tauQGDf    = linearInterpolate(tauQGD);
alphauf    = linearInterpolate(turbulence->muEff());
muf        = linearInterpolate(turbulence->alphaEff());


