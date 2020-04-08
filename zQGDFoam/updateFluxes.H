//Gradients and divergence
//---------Start---------
gradUf = fvsc::grad(U);

divUf  = tr(gradUf);

gradef = fvsc::grad(e);

gradRhof = fvsc::grad(rho);

gradPf = fvsc::grad(p);

divRhoUf = fvsc::div(rhoU);

//---------End---------

//Continuity equation fluxes
//---------Start---------

surfaceVectorField wHatf = 
    (tauQGDf / rhof)
    *
    (
        ((rhof*Uf) & gradUf) + gradPf
    );
surfaceVectorField wf    =
  wHatf + (tauQGDf/rhof)*(Uf * divRhoUf);

jm = rhoLnf*(Uf - wf);
phiJm = mesh.Sf() & jm;
phi = mesh.Sf() & linearInterpolate(rhoU);
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = mesh.Sf() & (jm * Uf);
phiP = mesh.Sf() * pf;

Info << "gammaf=" << gammaf << endl;

//QGD diffusive fluxes
Pif = 
    rhof*Uf*wHatf
    +
    tauQGDf *
    (
        I*((Uf & gradPf) + (gammaf*p1f*divUf))
    );
    
if (!implicitDiffusion)
{
    Pif +=
        muf*
        (
            gradUf
//            +
//            Foam::T(gradUf)
//            -
//            (2.0/3.0)*I*divUf
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
phih2 = (0.25*hQGDf*hQGDf*(fvc::snGrad(U) * fvc::snGrad(p)) & mesh.Sf());

phiJmH = (E1f*(Uf - wf)) & mesh.Sf();
phih2.setOriented(phiJmH.oriented()());
phiJmH -= phih2;

qf =
    -tauQGDf*rhof*Uf*
    (
        (
            (Uf&gradef)
            -
            (
                Uf&
                (
                    (p1f*gradRhof)/rhof/rhof
                )
            )
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
