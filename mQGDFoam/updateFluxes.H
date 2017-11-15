//Gradients and divergence
//---------Start---------
gradPf = fvsc::grad(p);

gradUf = fvsc::grad(U);

gradUfT= Foam::T(gradUf);

divUf  = tr(gradUf);

gradef = fvsc::grad(e);

gradRhof = fvsc::grad(rho);

//---------End---------

//Continuity equation fluxes
//---------Start---------

//rhoW = tauQGDf * (fvsc::div(rhoUU) + gradPf); //rhoW should be zero at wall

rhoW = tauQGDf *
        (
            (gradUfT & rhofUf)
            +
            //(Uf * fvsc::div(rhoU))
            Uf * (Uf & gradRhof)
            +
            rhofUf * divUf
            +
            gradPf
        );

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
            gradUfT
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
        -
        (pf * (Uf & gradRhof) / rhof / rhof)
    );
if (!implicitDiffusion)
{
    qf -=
        alphauf*gradef;
}

phiQ = qf & mesh.Sf();

phiPiU = (Pif & Uf) & mesh.Sf();

//---------End---------
