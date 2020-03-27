/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
                QGDsolver   | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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


\*---------------------------------------------------------------------------*/

    gradPf = fvsc::grad(p);
    Wf = tauQGDf*((Uf & gradUf) + gradPf/rhof - BdFrcf);
    
    surfaceVectorField phiUfWf = mesh.Sf() & (Uf * Wf);
    phiUfWf.setOriented(true);
    phiUf = qgdFlux(phi,U,Uf);
    phiUf.setOriented(true);
    phiUf -= phiUfWf;

    // --- Solve U
    if (implicitDiffusion)
    {
        solve
        (
            fvm::ddt(U)
            +
            fvc::div(phiUf)
            -
            fvm::laplacian(muf/rhof,U)
            -
            fvc::div(muf/rhof * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))))
            ==
            -
            fvc::grad(p)/rho
            +
            BdFrc
            +
            USu
        );
    }
    else
    {
        solve
        (
            fvm::ddt(U)
            +
            fvc::div(phiUf)
            -
            fvc::laplacian(muf/rhof,U)
            -
            fvc::div(muf/rhof * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))))
            ==
            -
            fvc::grad(p)/rho
            +
            BdFrc
            +
            USu
        );
    }

// ************************************************************************* //
