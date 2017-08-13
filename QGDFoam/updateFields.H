/*
 *
 * Auxiliary fields
 *
 */

// Inversed compressibility
rPsi = 1.0 / psi;

// Heat capacities ratio
gamma = thermo.Cp()/thermo.Cv();

// Speed of sound
c = sqrt(gamma * rPsi);

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

// Heat capacities ratio
gammaf = linearInterpolate(gamma);
gammam1f = gammaf - 1;

// Speed of sound
cf = linearInterpolate(c);

// Heat capacity at constant pressure
Cpf = linearInterpolate(thermo.Cp());
    
// Total specific enthalpy
Hf = linearInterpolate((rhoE + p)/rho);

/*
 *
 * QGD coefficients
 *
 */
tauQGD = alphaQGD * hQGD / cf;
    
muQGD = tauQGD*pf*ScQGD;
    
kappaQGD = muQGD*Cpf / PrQGD;

/*
 *
 * Viscous part of N-S equations
 *
 */
if (!inviscid)
{
    muEffPtr() = turbulence->muEff()();
    tauMCPtr() = muEffPtr()*dev2(Foam::T(fvc::grad(U)));
    sigmaDotUPtr() = 
        (
            linearInterpolate(muEffPtr())*mesh.magSf()*fvc::snGrad(U)
          + fvc::dotInterpolate(mesh.Sf(), tauMCPtr())
        )
        & Uf;
}
