//Gradients and divergence
//---------Start---------
gradPf = fvsc::grad(p);

gradUf = fvsc::grad(U);

divUf  = fvsc::div(U);

gradef = fvsc::grad(e);

gradInvRhof = fvsc::grad(1/rho);
//---------End---------

//Continuity equation fluxes
//---------Start---------
phivf = Uf & mesh.Sf();

rhoW = tauQGDf * (fvsc::div(rhoUU) + gradPf);

//force rhoW to be zero at wall patches
forAll(mesh.boundary(), patchi)
{
    if (isA<wallFvPatch>(mesh.boundary()[patchi]))
    {
        rhoW.boundaryFieldRef()[patchi] = vector::zero;
    }
}

//divergent scheme
//rhoW = tauQGDf * 
//        (
//            (faceStencil.faceVectorGrad(rhoU) &  Uf)
//            + 
//            (rhofUf * divUf)
//            + gradPf 
//        );


jm = rhofUf - rhoW;

phiJm = jm & mesh.Sf();
phi = phiJm;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = (jm * Uf) & mesh.Sf();

phiP = pf*mesh.Sf();

Pif = 
    //QGD diffusive fluxes
    tauQGDf *
    (
        Uf * (rhof * (Uf & gradUf) + gradPf)
        +
        I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
    );
    
if (!implicitDiffusion)
{
    Pif +=
        muf*
        (
            gradUf
            +
            Foam::T(gradUf)
            -
            (2.0/3.0)*I*divUf
        );
}
else
{
    tauMCPtr() = muf*linearInterpolate(Foam::T(fvc::grad(U)) - (2.0/3.0)*I*fvc::div(U));
    phiTauMC = tauMCPtr() & mesh.Sf();
}

phiPi = Pif & mesh.Sf();

//---------End---------

// Fluxes for energy balance equation
//---------Start---------
phiJmH = phiJm * Hf;

qf =
    -
    tauQGDf*rhofUf*
    (
        (Uf & gradef)
        +
        (pf * (Uf & gradInvRhof))
    );
if (!implicitDiffusion)
{
    qf -=
        alphauf*gradef;
}

phiQ = qf & mesh.Sf();

phiPiU = (Pif & Uf) & mesh.Sf();

//---------End---------
