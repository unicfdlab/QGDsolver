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
    This file is part of QGDsolver, based on OpenFOAM library.

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
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "tmp.H"

#ifndef fvsc_H
#define fvsc_H

namespace Foam
{

namespace fvsc
{

//returns gradient of volScalarField at the surface mesh
tmp<surfaceVectorField> grad(const volScalarField& vF);

//-
tmp<surfaceVectorField> grad(const tmp<volScalarField>& vF);

//returns gradient of volVectorField at the surface mesh
tmp<surfaceTensorField> grad(const volVectorField& vF);

//-
tmp<surfaceTensorField> grad(const tmp<volVectorField>& vF);

//returns divergence of volVectorField at the surface mesh
tmp<surfaceScalarField> div(const volVectorField& vF);

//-
tmp<surfaceScalarField> div(const tmp<volVectorField>& vF);

//returns divergence of volTensorField at the surface mesh
tmp<surfaceVectorField> div(const volTensorField& vF);

//-
tmp<surfaceVectorField> div(const tmp<volTensorField>& vF);

} //namespace fvsc

} //namespace Foam

#endif

//END-OF-FILE

