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
    updateFields

Description
    Updates fields.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
 *
 * Gradients at faces
 *
 */
gradUf = fvsc::grad(U);

gradAlpha1f = fvsc::grad(alpha1);

muf = qgdInterpolate
(
    alpha1*(rho1*nu1-rho2*nu2) + rho2*nu2
);

Info << "max/min muf:" << max(muf).value() << "/" << min(muf).value() << endl;

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */
// Velocity
Uf = qgdInterpolate(U);
    
// alphal
alpha1f = qgdInterpolate(alpha1);
alpha2f = 1.0 - alpha1f;

surfaceScalarField da1dtf = -(Uf & gradAlpha1f);
volScalarField da1dt = -(U & fvc::grad(alpha1));

cFrc   = fvc::reconstruct(thermo.surfaceTensionForce()*mesh.magSf());
cFrcf  = thermo.surfaceTensionForce()*mesh.Sf() / mesh.magSf();

BdFrc  = g * (rho + da1dt*(Tau1*rho1 - Tau2*rho2));


// ************************************************************************* //
