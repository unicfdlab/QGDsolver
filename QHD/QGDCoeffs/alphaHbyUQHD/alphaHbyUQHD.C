/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                QGDsolver   | Copyright (C) 2016-2017 ISP RAS (www.unicfd.ru)
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
    alphaHbyUQHD

Description
    Class for one of possible ways of tau calculating.

SourceFiles
    alphaHbyUQHD.C

\*---------------------------------------------------------------------------*/



#include "alphaHbyUQHD.H"
#include "rhoQGDThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(alphaHbyUQHD,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        alphaHbyUQHD,
        dictionary
    );
}
}

Foam::qgd::
alphaHbyUQHD::alphaHbyUQHD
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict),
    uQGD_(0.0)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;
    
//    
    dict.lookup("uQGD") >> uQGD_;
    
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    muQGD_.primitiveFieldRef() = 0.0;
    alphauQGD_.primitiveFieldRef() = 0.0;
    
    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    muQGD_.boundaryFieldRef() = 0.0;
    alphauQGD_.boundaryFieldRef() = 0.0;
}

Foam::qgd::
alphaHbyUQHD::~alphaHbyUQHD()
{
}

void Foam::qgd::
alphaHbyUQHD::correct(const Foam::rhoQGDThermo& qgdThermo)
{
    this->tauQGD_ = this->aQGD_ * this->hQGD_ /dimensionedScalar("uQGD", dimLength/dimTime, uQGD_);
}

//
//END-OF-FILE
//


