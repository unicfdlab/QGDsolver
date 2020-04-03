/*
 *
 * Gradients at faces
 *
 */
Urel.correctBoundaryConditions();
gradUrelf = fvsc::grad(Urel);

/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
muf    = linearInterpolate(turbulence->mu());
alphaf = linearInterpolate(turbulence->alpha());

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rho.correctBoundaryConditions();
rhof = linearInterpolate(rho);

// Velocity
Urelf = linearInterpolate(Urel);
    
// T
Tf = linearInterpolate(T);


//rhok   = 1.0 - beta*T;

BdFrc  = -2.0*(SRF->omega()^Urel);

//BdFrc  = (rhok - rhok0)*g;
BdFrc.correctBoundaryConditions();
BdFrcf = fvc::interpolate(BdFrc);

//Hi
Hif = alphaf/rhof;

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());




