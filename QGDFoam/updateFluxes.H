/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of QGDsolver library, based on OpenFOAM+.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    updateFluxes

Description
    Updates fluxes for continuity equation, momentum balance equation,
    energy balance equation.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = qgdFlux(phiJm,U,Uf);
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

phiPiU = mesh.Sf() & (Pif & Uf);

//---------End---------


// ************************************************************************* //
