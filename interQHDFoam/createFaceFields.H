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
    createFaceFields

Description
    Creates the face fields: linear interpolation of fields from volumes 
    to face centers.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Density

// Liquid fraction
surfaceScalarField alpha1f
(
    "alpha1f",
    linearInterpolate(alpha1)
);
surfaceScalarField alpha2f
(
    "alpha2f",
    1.0 - alpha1f
);

// Velocity
surfaceVectorField Uf
(
    "Uf",
    linearInterpolate(U)
);
surfaceVectorField W1
(
    "Wf",
    linearInterpolate(U)*0.0
);
surfaceVectorField W2
(
    "Wf",
    linearInterpolate(U)*0.0
);

// Forces
surfaceVectorField BdFrcf
(
    "BdFrcf",
    linearInterpolate(BdFrc)
);

surfaceVectorField cFrcf
(
    "cFrcf",
    linearInterpolate(cFrc)
);

// fluid viscosity
surfaceScalarField muf
(
    "muf",
    linearInterpolate(thermo.mu()*0.0)
);
    

// ************************************************************************* //