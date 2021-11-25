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
    createFaceFluxes

Description
    Creates the face-flux fields.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Gradients and divergence
//---------Start---------
surfaceVectorField gradPf 
(
    "gradPf", fvsc::grad(p)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceScalarField divUf
(
    "divUf",
    Foam::tr(gradUf)
);

surfaceVectorField gradef
(
    "gradef",
    fvsc::grad(e)
);

surfaceVectorField gradRhof
(
    "gradRhof",
    fvsc::grad(rho)
);

//---------End---------

//Continuity equation fluxes
//---------Start---------

surfaceVectorField rhoUf
(
    "rhoUf",
    linearInterpolate(rho*U)
);

surfaceVectorField rhofUf
(
    "rhofUf",
    rhof*Uf
);

surfaceTensorField UrhoUf
(
    "UrhoUf",
    linearInterpolate(U*rhoU)
);

phi = rhoUf & mesh.Sf();

surfaceVectorField rhoW
(
    "rhoW",
    tauQGDf*gradPf
);

surfaceScalarField phiw
(
    "phiwStar",
    mesh.Sf() & rhoW
);

surfaceVectorField jm
(
    "jm",
    rhoUf - rhoW
);

surfaceScalarField phiJm
(
    "phiJm",
    mesh.Sf() & jm
);

//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiJmU
(
    "phiJmU",
    mesh.Sf() & (jm * Uf)
);

surfaceVectorField phiP
(
    "phiP",
    pf*mesh.Sf()
);

surfaceTensorField Pif
(
    "Pif",
    //QGD diffusive fluxes
    tauQGDf * ((UrhoUf & gradUf) + Uf * gradPf)
    +
    tauQGDf * I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
);

surfaceVectorField phiPi
(
    "phiPi",
    mesh.Sf() & Pif
);

autoPtr<surfaceTensorField> tauMCPtr;
tauMCPtr.reset
(
    new surfaceTensorField
    (
        "tauMC",
        0.0*muf*
        (
            Foam::T(gradUf)//Don't forget to transpose
            -
            I*divUf
        )
    )
);

surfaceVectorField phiTauMC
(
    "phiTauMC",
    tauMCPtr() & mesh.Sf()
);

//---------End---------

// Fluxes for energy balance equation
//---------Start---------
surfaceScalarField phiJmH
(
    "phiJmH",
    phiJm * Hf
);

surfaceVectorField qf
(
    "qf",
    -
    tauQGDf* 
    ( 
        ((linearInterpolate(rho*(U*U))) & gradef)
        -
        (pf * rhof * Uf * (Uf & gradRhof) / rhof / rhof)
    )
);

surfaceScalarField phiQ
(
    "phiQ",
    qf & mesh.Sf()
);

surfaceScalarField phiPiU
(
    "phiPiU",
    (Pif & Uf) & mesh.Sf()
);

autoPtr<surfaceVectorField> sigmaDotUPtr;
sigmaDotUPtr.reset
(
    new surfaceVectorField
    (
        "sigmaDotU",
        tauMCPtr() & linearInterpolate(U) //to be updated later
    )
);

surfaceScalarField phiSigmaDotU
(
    "phiSigmaDotU",
    sigmaDotUPtr() & mesh.Sf()
);

//---------End---------


// ************************************************************************* //