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
    grpTwoPhaseIcoQGDThermo

\*---------------------------------------------------------------------------*/

#include "twoPhaseIcoQGDThermo.H"
#include "QGDCoeffs.H"
#include "twoPhaseConstTau.H"

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
    constTwoPhaseProperties(*this, phase1Name(), phase2Name()),
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
    )
{
    validateQGDCoeffs();
    this->read();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseIcoQGDThermo::~twoPhaseIcoQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseIcoQGDThermo::validateQGDCoeffs()
{
    if (!isA<qgd::twoPhaseConstTau>(this->qgdCoeffs()))
    {
        FatalErrorInFunction
        << "twoPhaseIcoQGDThermo works only with twoPhaseConstTau QGD Coeffs Model"
        << exit(FatalError);
    }
}

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
  return (rho1()-rho2())*this->alpha1() + rho2();
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseIcoQGDThermo::mu() const
{
  return (rho1()*nu1()-rho2()*nu2())*this->alpha1() + rho2()*nu2();
}

void Foam::twoPhaseIcoQGDThermo::correct()
{
    qgdCoeffs().correct(*this);
    qInterfaceProperties::correct();
}



// ************************************************************************* //
