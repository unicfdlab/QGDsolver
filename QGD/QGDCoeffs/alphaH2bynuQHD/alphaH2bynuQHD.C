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
    alphaH2bynuQHD

Description
    Class for one of possible ways of tau calculating.

SourceFiles
    alphaH2bynuQHD.C

\*---------------------------------------------------------------------------*/

#include "alphaH2bynuQHD.H"
#include "QGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(alphaH2bynuQHD,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        alphaH2bynuQHD,
        dictionary
    );
}
}

Foam::qgd::
alphaH2bynuQHD::alphaH2bynuQHD
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;
//

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
alphaH2bynuQHD::~alphaH2bynuQHD()
{
}

void Foam::qgd::
alphaH2bynuQHD::correct(const QGDThermo& qgdThermo)
{
    const volScalarField  nu     = qgdThermo.mu()/qgdThermo.rho();
    this->tauQGD_ = this->aQGD_ * sqr(this->hQGD_) / nu;
    this->tauQGDf_= linearInterpolate(this->tauQGD_);
}

//
//END-OF-FILE
//
