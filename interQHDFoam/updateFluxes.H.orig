//Continuity equation fluxes
//---------Start---------
phiu  = mesh.Sf() & Uf;
phiu.setOriented(true);
{
    surfaceScalarField phiwon = 
        mesh.Sf() & ((Uf & gradUf) - g);
    phiwon.setOriented(true);
    surfaceScalarField phicf = mesh.Sf() & cFrcf;
    phicf.setOriented(true);
    
    phiwo1 = (phiwon - phicf/rho1) * Tau1;
    phiwo2 = (phiwon - phicf/rho2) * Tau2;
}

//---------End---------
