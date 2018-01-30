//Gradients and divergence
//---------Start---------
gradPf = fvsc::grad(p);

gradUf = fvsc::grad(U);

//gradRhoUf = fvsc::grad(rhoU);
//gradUfT= Foam::T(gradUf);

divUf  = tr(gradUf);

gradef = fvsc::grad(e);

gradRhof = fvsc::grad(rho);

//---------End---------

//Continuity equation fluxes
//---------Start---------

//rhoW = tauQGDf * (fvsc::div(rhoU*U) + gradPf); //rhoW should be zero at wall
//rhoW = tauQGDf * 
//        (
//            (rhofUf * divUf)
//            +
//            (Uf & gradRhoUf)
//            +
//            gradPf
//        );

rhoWStar = tauQGDf * ((rhofUf & gradUf) + gradPf);

rhoW = tauQGDf *
        (
            ((Uf * gradRhof) & Uf)
            +
            (rhoUf * divUf) //was rhoUf
        )
        +
        rhoWStar;

jm = rhoUf - rhoW; //was rhoUf

phiJm = mesh.Sf() & jm;
phi   = mesh.Sf() & rhofUf;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = mesh.Sf() & (jm * Uf);
//phiJmU = Foam::T(jm * Uf) & mesh.Sf();

phiP = mesh.Sf() * pf;

Pif = 
    //QGD diffusive fluxes
    Uf * rhoWStar
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
