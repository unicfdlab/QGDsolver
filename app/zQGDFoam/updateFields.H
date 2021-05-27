/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);
rho_pos = interpolate(rho,pos);
rho_neg = interpolate(rho,neg);
rhoLnf  = 1.0 / logMean(rho_pos,rho_neg);

ef = linearInterpolate(e);
e_pos = interpolate(e,pos);
e_neg = interpolate(e,neg);
eLnf  = e_pos*e_neg*logMean(e_pos,e_neg);

// Velocity
Uf = linearInterpolate(U);
U_pos = interpolate(U,pos);
U_neg = interpolate(U,neg);

// Pressure
pf = linearInterpolate(p);

//Speed of sound
cf = linearInterpolate(thermo.c());

// Heat capacities ratio
gammaf = linearInterpolate(gamma);
gammam1f = gammaf - 1;

// Total specific enthalpy
E1f = (0.5*rhoLnf*(U_pos&U_neg)+rhoLnf*eLnf+pf);
p1f = gammam1f*rhof*ef;

/*
 *
 * QGD coefficients are now updated in the thermo library
 *
 */
alphauf    = linearInterpolate(turbulence->alphaEff());
muf        = linearInterpolate(turbulence->muEff());

