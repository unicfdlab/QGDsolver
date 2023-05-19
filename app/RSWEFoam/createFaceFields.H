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

// Velocity
surfaceVectorField Uf
(
    IOobject
    (
        "Uf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U)
);

// Height

surfaceScalarField hf
(
    IOobject
    (
        "hf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(h)
);

// Height^2
surfaceScalarField h2f
(
    "h2f",
    linearInterpolate(h*h)
);

// bathymetry
surfaceScalarField bf
(
    "bf",
    linearInterpolate(b)
);

surfaceScalarField tauf
(
    IOobject
    (
        "tauf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(tau)
);

surfaceScalarField nuf
(
    "nuf",
    tauf*magg*hf
);

surfaceVectorField phih2
(
    "phih2",
    hf*hf*mesh.Sf()
);

surfaceVectorField phib
(
    "phib",
    bf*mesh.Sf()
);

const faceList& faces = mesh.faces();
const labelList& faceOwner = mesh.faceOwner();          // Face to owner cell
const labelList& faceNeighbour = mesh.faceNeighbour();  // Face to neighbour cell 

const fvPatchList& patches = mesh.boundary();

// ************************************************************************* //