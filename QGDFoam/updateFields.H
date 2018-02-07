/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);

// Velocity
Uf = linearInterpolate(U);
Wf = tauQGD*(Uf & gradUf + gradPf/rho + betaf*g*Tf);
    
// Pressure
pf = linearInterpolate(p);

// Heat capacities ratio
gammaf = linearInterpolate(gamma);
gammam1f = gammaf - 1;

// T
Tf = lineatInterpolate(T);

//Hi
Hif = linearInterpolate(Hi);

// Speed of sound
cf = linearInterpolate(c);

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
tauQGDf    = linearInterpolate(tauQGD);
alphauf    = linearInterpolate(turbulence->alphaEff());
muf        = linearInterpolate(turbulence->muEff());
betaf      = linearInterpolate(beta);

