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

surfaceVectorField gradAlpha1f
(
    "gradAlpha1f", fvsc::grad(alpha1)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceVectorField gradpf
(
    "gradpf",
    fvsc::grad(p)
);

surfaceScalarField rhoPhi
(
    "rhoPhi",
    rho1*phi
);

surfaceVectorField phiUfRhof
(
    "phiUf",
    phi * Uf * rho1
); 
phiUfRhof.setOriented(true);

surfaceVectorField phiRhofWf
(
    "phiUf",
    phi * Uf * rho1
);
phiRhofWf.setOriented(true);

surfaceScalarField phiAlpha1f
(
    "phiAlpha1f",
    phi*alpha1f
);

surfaceScalarField phiAlpha2f
(
    "phiAlpha2f",
    phi*alpha1f
);

surfaceScalarField phiwo1
(
    "phiwo1",
    phi*alpha1f
);
surfaceScalarField phiwo2
(
    "phiwo2",
    phi*alpha1f
);
surfaceScalarField phiw1
(
    "phiw1",
    phi*alpha1f
);
surfaceScalarField phiw2
(
    "phiw2",
    phi*alpha1f
);

surfaceScalarField phiwm
(
    "phiwm",
    phiwo1*0.0
);
surfaceScalarField coeffp
(
    "coeffp",
    linearInterpolate(Tau1 / rho)
);

surfaceScalarField phiWr
(
    "phiWr",
    phiwo1*0.0
);

surfaceScalarField phiu
(
    "phiu",
    mesh.Sf() & linearInterpolate(U)
);
surfaceScalarField phi1
(
    "phi1",
    mesh.Sf() & linearInterpolate(U)
);
surfaceScalarField phi2
(
    "phi2",
    mesh.Sf() & linearInterpolate(U)
);


// ************************************************************************* //