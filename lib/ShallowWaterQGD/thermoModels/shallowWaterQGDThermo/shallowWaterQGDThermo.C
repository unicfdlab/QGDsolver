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
    
Group
    grpShallowWaterQGDThermo

\*---------------------------------------------------------------------------*/

#include "shallowWaterQGDThermo.H"
#include "QGDCoeffs.H"
#include "RSWETau.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(shallowWaterQGDThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shallowWaterQGDThermo::shallowWaterQGDThermo(const fvMesh& mesh)
:
    IOdictionary
    (
        IOobject
        (
            "shallowWaterProperties",
            mesh.time().constant(),
            mesh.time().db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    QGDThermo(mesh, *this)
{
    validateQGDCoeffs();
    this->read();
    this->lookup("alpha") >> alpha_;
    
    wellBalancedScheme_ = false;
    dryZoneCondition_ = false;
    NS_ = 1.0;
    eps0_ = 1e-6;
    tauU_ = 0;
    n_ = 0.0;

    if (this->found("wellBalancedScheme"))
    {
        this->lookup("wellBalancedScheme") >> wellBalancedScheme_;
    }

    if (this->found("NS"))
    {
        this->lookup("NS") >> NS_;
    }

    if (this->found("dryZoneCondition"))
    {   
        this->lookup("dryZoneCondition") >> dryZoneCondition_;
    }

    if (this->found("eps0"))
    {   
        this->lookup("eps0") >> eps0_;
    }

    if (this->found("tauU"))
    {   
        this->lookup("tauU") >> tauU_;
    }
    
    if (this->found("n"))
    {   
        this->lookup("n") >> n_;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shallowWaterQGDThermo::~shallowWaterQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shallowWaterQGDThermo::validateQGDCoeffs()
{
    if (!isA<qgd::RSWETau>(this->qgdCoeffs()))
    {
        FatalErrorInFunction
        << "shallowWaterQGDThermo works only with RSWETau QGD Coeffs Model"
        << exit(FatalError);
    }
}

bool Foam::shallowWaterQGDThermo::read()
{

    if (!QGDThermo::read())
    {
      return false;
    }

    return true;
}

const Foam::scalar& Foam::shallowWaterQGDThermo::alpha() const
{
    return alpha_;
}

const Foam::scalar& Foam::shallowWaterQGDThermo::NS() const
{
    return NS_;
}

const Foam::scalar& Foam::shallowWaterQGDThermo::eps0() const
{
    return eps0_;
}

const Foam::scalar& Foam::shallowWaterQGDThermo::tauU() const
{
    return tauU_;
}

const Foam::scalar& Foam::shallowWaterQGDThermo::n() const
{
    return n_;
}

const Foam::Switch& Foam::shallowWaterQGDThermo::dryZoneCondition() const
{
    return dryZoneCondition_;
}

const Foam::Switch& Foam::shallowWaterQGDThermo::wellBalancedScheme() const
{
    return wellBalancedScheme_;
}

const Foam::volScalarField& Foam::shallowWaterQGDThermo::c() const
{
    notImplemented ("const volScalarField& Foam::shallowWaterQGDThermo::c() const");
    return volScalarField::null();
}

Foam::tmp<Foam::volScalarField> Foam::shallowWaterQGDThermo::rho() const
{
    notImplemented ("const volScalarField& Foam::shallowWaterQGDThermo::rho() const");
    return volScalarField::null();
}

Foam::tmp<Foam::volScalarField> Foam::shallowWaterQGDThermo::mu() const
{
    notImplemented ("const volScalarField& Foam::shallowWaterQGDThermo::mu() const");
    return volScalarField::null();
}

const Foam::volScalarField& Foam::shallowWaterQGDThermo::p() const
{
    notImplemented ("const volScalarField& Foam::shallowWaterQGDThermo::p() const");
    return volScalarField::null();
}


void Foam::shallowWaterQGDThermo::correct()
{
    qgdCoeffs().correct(*this);
}



// ************************************************************************* //
