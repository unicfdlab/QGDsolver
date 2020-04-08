/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = qgdInterpolate(rho);

// Velocity
Uf = qgdInterpolate(U);

// Momentum
rhoUf = qgdInterpolate(rhoU);
rhofUf= rhof*Uf;

// Nonlinear Tensor
UrhoUf = qgdInterpolate(U*rhoU);

// Pressure
pf = qgdInterpolate(p);

//Speed of sound
cf = qgdInterpolate(thermo.c());

// Heat capacities ratio
gammaf = qgdInterpolate(gamma);
gammam1f = gammaf - 1;

// Heat capacity at constant pressure
Cpf = qgdInterpolate(thermo.Cp());
    
// Total specific enthalpy
H = (rhoE + p)/rho;
Hf = qgdInterpolate(H);

/*
 *
 * QGD coefficients are now updated in the thermo library
 *
 */
alphauf    = qgdInterpolate(turbulence->alphaEff());
muf        = qgdInterpolate(turbulence->muEff());




