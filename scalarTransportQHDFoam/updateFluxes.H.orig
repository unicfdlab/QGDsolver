//Continuity equation fluxes
//---------Start---------
Uf    = qgdInterpolate(U);
gradUf= fvsc::grad(U);
phiu  = mesh.Sf() & Uf;
phiu.setOriented(true);
phiwo = mesh.Sf() & (tauQGDf*(Uf & gradUf));
phiwo.setOriented(true);

taubyrhof = tauQGDf/rhof;
//---------End---------
