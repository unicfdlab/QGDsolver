//Gradients and divergence
//---------Start---------
gradUf = fvsc::grad(U);

gradWf = fvsc::grad(W);

//gradUfT= Foam::T(gradUf);

divUf  = tr(gradUf);

gradRhof = fvsc::grad(rho);

//---------End---------

//Continuity equation fluxes
//---------Start---------

//rhoW = tauQGDf * (fvsc::div(rhoUU) + gradPf); //rhoW should be zero at wall

phiu  = mesh.Sf() & U;
phiwo = mesh.Sf() & tauQGD*(Uf & gradUf + betaf*g*Tf);

taubyrho = tauQGD/rho;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiPi = muf * (mesh.Sf() & ((gradUf+gradUf.T()))/rhof);

BdFrc = betaf*g*Tf;

//---------End---------

// Fluxes for T
//---------Start---------

//---------End---------
