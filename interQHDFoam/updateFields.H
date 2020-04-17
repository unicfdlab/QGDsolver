/*
 *
 * Gradients at faces
 *
 */
gradUf = fvsc::grad(U);

gradAlpha1f = fvsc::grad(alpha1);

muf = qgdInterpolate
(
    alpha1*(rho1*nu1-rho2*nu2) + rho2*nu2
);

Info << "max/min muf:" << max(muf).value() << "/" << min(muf).value() << endl;

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */
// Velocity
Uf = qgdInterpolate(U);
    
// alphal
alpha1f = qgdInterpolate(alpha1);
alpha2f = 1.0 - alpha1f;

surfaceScalarField da1dtf = -(Uf & gradAlpha1f);
volScalarField da1dt = -(U & fvc::grad(alpha1));

cFrc   = fvc::reconstruct(thermo.surfaceTensionForce()*mesh.magSf());
cFrcf  = qgdInterpolate(cFrc);
BdFrc  = g * (rho + da1dt*(Tau1*rho1 - Tau2*rho2));


