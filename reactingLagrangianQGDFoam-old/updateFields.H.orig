/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);

//
yf = linearInterpolate(y);

// Velocity
Uf = linearInterpolate(U);

// Momentum
rhoUf = linearInterpolate(rhoU);
rhofUf= rhof*Uf;

// Nonlinear Tensor
UrhoUf = linearInterpolate(U*rhoU);

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
//Hf = (linearInterpolate(rhoE) + pf)/rhof;
Hf = linearInterpolate((rhoE + p)/rho);

/*
 *
 * QGD coefficients are now updated in the thermo library
 *
 */
alphauf    = linearInterpolate(turbulence->alphaEff());
muf        = linearInterpolate(turbulence->muEff());

//surfaceScalarField tauGammaPf
//(
//    "tauGammaPf",
//    linearInterpolate(tauQGD*gamma*p)
//);
//
//surfaceVectorField tauUf
//(
//    "tauUf",
//    linearInterpolate(tauQGD*U)
//);
//
//surfaceTensorField tauUrhoUf
//(
//    "tauUrhoUf",
//    linearInterpolate(tauQGD*U*rhoU)
//);

//tauQGDf    = 0.5* hQGDf / cf;
//muf        = pf * tauQGDf;
//alphauf    = muf*1.0;


