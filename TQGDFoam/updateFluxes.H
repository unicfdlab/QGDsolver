//Gradients and divergence
//---------Start---------
gradPf = fvsc::grad(p);

gradUf = fvsc::grad(U);

if (!implicitDiffusion)
{
    gradef = fvsc::grad(e);
    divUf  = tr(gradUf);
}

//---------End---------

//Continuity equation fluxes
//---------Start---------
rhoW = tauQGDf * ((rhofUf & gradUf) + gradPf);

//jm = rhoUf - rhoW; //was rhoUf
jm = rhofUf - rhoW; //was rhoUf

phiJm = mesh.Sf() & jm;
phi   = mesh.Sf() & rhofUf;
//---------End---------

//Momentum balance equation fluxes
//---------Start---------

phiJmU = mesh.Sf() & (jm * Uf);

phiP = mesh.Sf() * pf;

Pif = 
    //QGD diffusive fluxes
    Uf * rhoW;
    
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
    //tauMCPtr() = muf*linearInterpolate(Foam::T(fvc::grad(U)) - (2.0/3.0)*I*fvc::div(U));
    //tauMCPtr() = muf*(Foam::T(gradUf) - (2.0/3.0)*I*divUf);
    tauMCPtr() = linearInterpolate(turbulence->muEff() * dev2(Foam::T(fvc::grad(U))));
    phiTauMC = mesh.Sf() & tauMCPtr();
}

phiPi = mesh.Sf() & Pif;

//---------End---------

// Fluxes for energy balance equation
//---------Start---------
phiJmH = phiJm * Hf;

if (!implicitDiffusion)
{
    qf = -alphauf*gradef;
}

phiQ = mesh.Sf() & qf;

phiPiU = mesh.Sf() & (Pif & Uf);

//---------End---------
