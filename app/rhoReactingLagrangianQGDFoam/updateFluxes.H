//Gradients and divergence
//---------Start---------
gradUf = fvsc::grad(U);

divUf  = tr(gradUf);

gradef = fvsc::grad(e); //*Cpf/gammaf;

gradRhof = fvsc::grad(rho);

gradPf = fvsc::grad(p);

//---------End---------

//Continuity equation fluxes
//---------Start---------
rhoW = tauQGDf*
        (
            ((Uf * gradRhof) & Uf)
            +
            (rhoUf * divUf)
            +
            (rhoUf & gradUf)
        );

rhoW += tauQGDf * gradPf;

jm = rhoUf - rhoW;

phiJm = mesh.Sf() & jm;
phi = mesh.Sf() & rhoUf;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = qgdFlux(phiJm, U, Uf);
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
    tauMCPtr() = qgdInterpolate(turbulence->muEff() * dev2(Foam::T(fvc::grad(U))));
    phiTauMC = mesh.Sf() & tauMCPtr();
}

phiPi = mesh.Sf() & Pif;

//---------End---------

// Fluxes for energy balance equation
//---------Start---------
phiJmH = qgdFlux(phiJm,H,Hf);

qf =
    -tauQGDf*
    (
        UrhoUf & 
        (
            gradef
            - (pf / rhof / rhof) * gradRhof
        )
    );
    
if (!implicitDiffusion)
{
    qf -=
        alphauf*gradef;
}

phiQ = mesh.Sf() & qf;

// QGD-141/142 bug hotfix
phiQ.setOriented(true);

phiPiU = mesh.Sf() & (Pif & Uf);

//clean diffusive flux for inert specie
{
    forAll(diffusiveFlux[inertIndex].internalField(), iFace)
    {
        diffusiveFlux[inertIndex].primitiveFieldRef()[iFace] = 0.0;
    }
    forAll(diffusiveFlux[inertIndex].boundaryFieldRef(), iPatch)
    {
        forAll(diffusiveFlux[inertIndex].boundaryFieldRef()[iPatch], iFace)
        {
            diffusiveFlux[inertIndex].boundaryFieldRef()[iPatch][iFace] = 0.0;
        }
    }
}
forAll(Y, i)
{
    
    if ( i != inertIndex )
    {
        gradYf = fvsc::grad(Y[i]);
        phiJmY[i] = qgdFlux(phiJm,Y[i],Yf[i]);
        surfaceScalarField dydtflux = - phi * tauQGDf
            * (Uf & gradYf);
        phiJmY[i] += dydtflux;
        diffusiveFlux[i] = 0.0*dydtflux;
        //inert specie flux will be updated in QGDYEqn.H
    }
    // QGD-141/142 bug hotfix
    diffusiveFlux[i].setOriented(true);
}

//---------End---------
