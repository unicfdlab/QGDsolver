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

BdFrcf = beta*g*Tf;

phiu  = mesh.Sf() & Uf;

phiwo = mesh.Sf() & (tauQGDf*((Uf & gradUf) + BdFrcf));

taubyrhof = tauQGDf/rhof;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------


phiPi = muf/rhof * mesh.Sf() & (gradUf+gradUf.T());


//---------End---------

// Fluxes for T
//---------Start---------

//---------End---------
