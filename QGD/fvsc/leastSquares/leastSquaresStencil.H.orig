
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
Class
   leastSquares 
Discription 
    Realisation least squares method for calculationg of defferential operators
Source file
    leastSquaresStencil.C 
    


\*---------------------------------------------------------------------------*/


#ifndef leastSquaresStencil_H
#define leastSquaresStencil_H

#include "fvscStencil.H"
#include "regIOobject.H"
#include "labelList.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "vector.H"
#include "List.H"

#include "processorFvPatch.H"
#include "leastSquaresBase.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class leastSquares
:
    public fvscStencil, public leastSquaresBase
{

public:

        TypeName ("leastSquares");

// Constructors

        //- Construct from IOobject. Optional flag for if IOobject is the
        //  top level regIOobject.
        leastSquares(const IOobject&);
        
        //- Destructor
        ~leastSquares();

        tmp<surfaceVectorField> Grad(const volScalarField& iF);

        tmp<surfaceTensorField> Grad(const volVectorField& iVF);

        tmp<surfaceScalarField> Div(const volVectorField& iVF);

        tmp<surfaceVectorField> Div(const volTensorField& iTF);
};

} //fvsc

} //Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

