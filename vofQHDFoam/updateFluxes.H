//Continuity equation fluxes
//---------Start---------
phiu  = mesh.Sf() & Uf;
phiu.setOriented(true);
phiwo = mesh.Sf() & (tauQGDf*((Uf & gradUf) - BdFrcf/rhof));
phiwo.setOriented(true);

taubyrhof = tauQGDf/rhof;
//---------End---------



// Fluxes for momentum balance equation
//---------Start---------

//phiPi = muf/rhof * mesh.Sf() & (gradUf+gradUf.T());

//---------End---------

// Fluxes for T
//---------Start---------

//---------End---------
