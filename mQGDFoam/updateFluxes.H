//Gradients and divergence
//---------Start---------
gradPf = fvsc::grad(p);

gradUf = fvsc::grad(U);

//gradUfT= Foam::T(gradUf);

divUf  = tr(gradUf);

gradef = fvsc::grad(e);

gradRhof = fvsc::grad(rho);

//---------End---------

//Continuity equation fluxes
//---------Start---------

//rhoW = tauQGDf * (fvsc::div(rhoUU) + gradPf); //rhoW should be zero at wall

rhoW = tauQGDf *
        (
            //(gradUfT & rhofUf)
            (rhofUf & gradUf)
            +
            //Uf * (Uf & gradRhof)
            ((Uf * gradRhof) & Uf)
            +
            (rhofUf * divUf)
            +
            gradPf
        );

jm = rhofUf - rhoW;

phiJm = mesh.Sf() & jm;
phi = phiJm;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = mesh.Sf() & (jm * Uf);
//phiJmU = Foam::T(jm * Uf) & mesh.Sf();

phiP = mesh.Sf() * pf;

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

phiPi = mesh.Sf() & Pif;

//---------End---------

// Fluxes for energy balance equation
//---------Start---------
phiJmH = phiJm * Hf;

qf =
    -
    tauQGDf*rhofUf*
    (
        (Uf & gradef)
        -
        (pf * (Uf & gradRhof) / rhof / rhof)
    );
if (!implicitDiffusion)
{
    qf -=
        alphauf*gradef;
}

phiQ = mesh.Sf() & qf;

phiPiU = mesh.Sf() & (Pif & Uf);

//---------End---------
