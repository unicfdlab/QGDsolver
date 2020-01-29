/*
 *
 * Gradients at faces
 *
 */
gradUf = fvsc::grad(U);

gradAlphalf = fvsc::grad(alphal);
/*
 *
 * QGD coefficients now updated in thermo library
 *
 */
muf = linearInterpolate
(
    alphal*mul + (1.0 - alphal)*mug
);
alphaf = linearInterpolate
(
    alphal*mul + (1.0 - alphal)*mug
);

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Density
rhof = linearInterpolate(rho);

// Velocity
rhoUf = linearInterpolate(rho*U);
Uf = rhoUf/rhof;
    
// alphal
alphalf = linearInterpolate(alphal);

volScalarField rhoStar =
    rho + thermo.tauQGD()*(U & fvc::grad(alphal))*(rhog-rhol);

BdFrc  = -fvc::grad(rhoStar)*(g & mesh.C()); //rho*g; //(rho + thermo.tauQGD()*(U & fvc::grad(alphal))*(rhog-rhol))*g;
BdFrcf = -fvsc::grad(rhoStar)*(g & mesh.Cf()); //rhof*g;
//BdFrc    = -fvc::grad(rho)*(g & mesh.C());// + thermo.tauQGD()*(U & fvc::grad(alphal))*(rhog-rhol)*g;
//BdFrcf   = -fvsc::grad(rho)*(g & mesh.Cf());// + tauQGDf*(Uf & gradAlphalf)*(rhog-rhol)*g;
BdFrcf.setOriented(true);

//Hi
Hif = alphaf/rhof;

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());




