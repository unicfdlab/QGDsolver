//Gradients and divergence
//---------Start---------
gradPf = fvsc::grad(p);

gradUf = fvsc::grad(U);

gradWf = fvc::grad(W);

//gradUfT= Foam::T(gradUf);

divUf  = tr(gradUf);

gradef = fvsc::grad(e);

gradRhof = fvsc::grad(rho);

//---------End---------

//Continuity equation fluxes
//---------Start---------

//rhoW = tauQGDf * (fvsc::div(rhoUU) + gradPf); //rhoW should be zero at wall

phiu  = mesh.Sf() & Uf;
phiwo = mesh.Sf() & Wf;
taubyrho = tauQGD/rho;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiUUmWUmWW = ((Uf * Uf)-(Wf * Uf)-(Wf * Wf)) & mesh.Sf();

Pbyrho = pf/rhof;

phiPi = muf*((gradUf+gradUf.T())) & mesh.Sf()/rhof;

phiBdFrc = beta*g*T*mesh.Sf();

//---------End---------

// Fluxes for T
//---------Start---------
phiUmWT = ((Uf - Wf)*Tf) & mesh.Sf();
HiT = Hif *Tf;

//---------End---------
