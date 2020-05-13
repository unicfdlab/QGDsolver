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
#include "constTwoPhaseProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constTwoPhaseProperties::constTwoPhaseProperties
(
    const dictionary& dict,
    const word& phase1name,
    const word& phase2name
)
:
    nu1_("nu"+phase1name,dict),
    nu2_("nu"+phase2name,dict),
    rho1_("rho"+phase1name,dict),
    rho2_("rho"+phase2name,dict),
    Tau1_("tau"+phase1name,dict),
    Tau2_("tau"+phase2name,dict)
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constTwoPhaseProperties::~constTwoPhaseProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::nu1() const
{
    return nu1_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::nu2() const
{
    return nu2_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::rho1() const
{
    return rho1_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::rho2() const
{
    return rho2_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::Tau1() const
{
    return Tau1_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::Tau2() const
{
    return Tau2_;
}


// ************************************************************************* //
