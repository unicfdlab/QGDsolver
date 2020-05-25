//Continuity equation fluxes
//---------Start---------
Uf    = qgdInterpolate(U);
phiu  = mesh.Sf() & Uf;
phiu.setOriented(true);

taubyrhof = tauQGDf/rhof;
//---------End---------
