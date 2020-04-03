//Gradients and divergence
//---------Start---------
gradUf = fvsc::grad(U);

divUf  = tr(gradUf);

gradef = fvsc::grad(e);

gradRhof = fvsc::grad(rho);

//---------End---------

//Continuity equation fluxes
//---------Start---------

rhoW = tauQGDf *
        (
            ((Uf * gradRhof) & Uf)
            +
            (rhoUf * divUf)
            +
            (rhoUf & gradUf)
        );

phiw = mesh.Sf() & rhoW;

gradPf = fvsc::grad(p);

rhoW += tauQGDf * gradPf;

jm = rhoUf - rhoW;

phiJm = mesh.Sf() & jm;

phi = mesh.Sf() & rhoUf;

phiJmyf = phiJm * yf;

surfaceScalarField dydt = 
    (-fvsc::div(rhoU*y) + yf*fvsc::div(rhoU))/rhof;
phiJmyf += (mesh.Sf() & (rhoUf*tauQGDf*dydt));

//phiJma  = 
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = mesh.Sf() & (jm * Uf);
phiP = mesh.Sf() * pf;

Pif = 
    //QGD diffusive fluxes
    tauQGDf*
    (
        (UrhoUf & gradUf)
        +
        Uf*gradPf
    )
    +
    tauQGDf *
    (
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

qf =
    -tauQGDf*
    (
        UrhoUf & 
        (
            gradef
            - (pf / rhof / rhof) * gradRhof
            //+ pf * fvsc::grad(1.0 / rho)
        )
    );
    
if (!implicitDiffusion)
{
    qf -=
        alphauf*gradef;
}

phiQ = mesh.Sf() & qf;

phiPiU = mesh.Sf() & (Pif & Uf);

//---------End---------
