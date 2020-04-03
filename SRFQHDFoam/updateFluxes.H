//Continuity equation fluxes
//---------Start---------
phiu  = mesh.Sf() & Urelf;
phiu.setOriented(true);

phiwo = mesh.Sf() & (tauQGDf*((Urelf & gradUrelf) - BdFrcf));
phiwo.setOriented(true);

taubyrhof = tauQGDf/rhof;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

//phiPi = muf/rhof * mesh.Sf() & (gradUf+gradUf.T());

//---------End---------

// Fluxes for T
//---------Start---------

//phiTf = phi*Tf;

//---------End---------
