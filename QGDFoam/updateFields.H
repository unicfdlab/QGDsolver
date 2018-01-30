/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);

// Velocity
Uf = linearInterpolate(U);

// Momentum
//rhofUf = rhof*Uf; //linearInterpolate(rhoU);
rhoUf = linearInterpolate(rhoU);
rhofUf = rhof*Uf;
//rhofUf = linearInterpolate(rhoU);

// Pressure
pf = linearInterpolate(p);

//Speed of sound
cf = linearInterpolate(thermo.c());

// Heat capacities ratio
gammaf = linearInterpolate(gamma);
gammam1f = gammaf - 1;

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());
    
// Total specific enthalpy
Hf = (linearInterpolate(rhoE) + pf)/rhof;

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
alphauf    = linearInterpolate(turbulence->alphaEff());
muf        = linearInterpolate(turbulence->muEff());



