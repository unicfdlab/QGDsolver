//Gradients and divergence
//---------Start---------
gradPf = faceStencil.faceScalarGrad(p);

gradUf = faceStencil.faceVectorGrad(U);

divUf = faceStencil.faceVectorDiv(U);
//---------End---------

//Continuity equation fluxes
//---------Start---------
phivf = Uf & mesh.Sf();
phi = phivf*rhof;

rhoW1 = tauQGD * faceStencil.faceTensorDiv(rho * (U * U));

phiRhoW1 = rhoW1 & mesh.Sf();

rhoW2 = tauQGD * gradPf;

phiRhoW2 = rhoW2 & mesh.Sf();

jm = Uf*rhof - rhoW1 - rhoW2;

phiJm = jm & mesh.Sf();
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
phiJmU = (jm * Uf) & mesh.Sf();

phiP = pf*mesh.Sf();

Pif = 
    //QGD diffusive fluxes
    tauQGD * 
    (
        Uf * (rhof * (Uf & gradUf) + gradPf)
        +
        I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
    )
    +
    //NS diffusive fluxes
    (
        muQGD*(gradUf + gradUf.T())
        -
        muQGD*I*divUf
    );

phiPi = Pif & mesh.Sf();
//---------End---------

// Fluxes for energy balance equation
//---------Start---------
phiJmH = phiJm * Hf;

qf = 
    -
    kappaQGD*faceStencil.faceScalarGrad(T)
    -
    tauQGD* 
    ( 
        ((linearInterpolate(rho*(U*U))) & faceStencil.faceScalarGrad(e))
        +
        (pf * rhof * Uf * (Uf & faceStencil.faceScalarGrad(1/rho)))
    );

phiQ = qf & mesh.Sf();

phiPiU = (Pif & Uf) & mesh.Sf();
//---------End---------
