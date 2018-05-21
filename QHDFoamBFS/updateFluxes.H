//Continuity equation fluxes
//---------Start---------
phiu  = mesh.Sf() & Uf;

phiwo = mesh.Sf() & (tauQGDf*((Uf & gradUf) + BdFrcf));

taubyrhof = tauQGDf/rhof;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

//phiPi = muf/rhof * mesh.Sf() & (gradUf+gradUf.T());

//---------End---------

// Fluxes for T
//---------Start---------

//---------End---------
