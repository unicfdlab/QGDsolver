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

#include "twoPhaseIcoQGDThermo.H"
#include "QGDCoeffs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseIcoQGDThermo, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseIcoQGDThermo::twoPhaseIcoQGDThermo(const fvMesh& mesh, const volVectorField& U)
:
    IOdictionary
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh.time().db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    twoPhaseMixture(mesh, *this),
    qInterfaceProperties(alpha1(),U,*this),
    QGDThermo(mesh, *this),
    p_
    (
        IOobject
        (
            "p",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho1_("rho"+this->phase1Name(),*this),
    rho2_("rho"+this->phase2Name(),*this),
    nu1_("nu"+this->phase1Name(), *this),
    nu2_("nu"+this->phase2Name(),*this),
    Tau1_("tau"+this->phase1Name(), *this),
    Tau2_("tau"+this->phase2Name(),*this)
{
    this->read();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseIcoQGDThermo::~twoPhaseIcoQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



bool Foam::twoPhaseIcoQGDThermo::read()
{

    if (!QGDThermo::read())
    {
      return false;
    }

    return true;
}

const Foam::volScalarField& Foam::twoPhaseIcoQGDThermo::p() const
{
    return p_;
}

Foam::volScalarField& Foam::twoPhaseIcoQGDThermo::p()
{
    return p_;
}

const Foam::volScalarField& Foam::twoPhaseIcoQGDThermo::c() const
{
    notImplemented ("const volScalarField& Foam::twoPhaseIcoQGDThermo::c() const");
    return volScalarField::null();
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseIcoQGDThermo::rho() const
{
  return (rho1_-rho2_)*this->alpha1() + rho2_;
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseIcoQGDThermo::mu() const
{
  return (rho1_*nu1_-rho2_*nu2_)*this->alpha1() + rho2_*nu2_;
}

void Foam::twoPhaseIcoQGDThermo::correct()
{
    qInterfaceProperties::correct();
}

const Foam::dimensionedScalar& Foam::twoPhaseIcoQGDThermo::nu1() const
{
    return nu1_;
}

const Foam::dimensionedScalar& Foam::twoPhaseIcoQGDThermo::nu2() const
{
    return nu2_;
}

const Foam::dimensionedScalar& Foam::twoPhaseIcoQGDThermo::rho1() const
{
    return rho1_;
}

const Foam::dimensionedScalar& Foam::twoPhaseIcoQGDThermo::rho2() const
{
    return rho2_;
}

const Foam::dimensionedScalar& Foam::twoPhaseIcoQGDThermo::Tau1() const
{
    return Tau1_;
}

const Foam::dimensionedScalar& Foam::twoPhaseIcoQGDThermo::Tau2() const
{
    return Tau2_;
}


// ************************************************************************* //
