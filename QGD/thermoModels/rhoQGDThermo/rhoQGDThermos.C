/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::qgd::psiQGDThermo

Description
    Possible properties

\*---------------------------------------------------------------------------*/

#include "rhoQGDThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "rhoConst.H"
#include "PengRobinsonGas.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "hPolynomialThermo.H"
#include "polynomialTransport.H"

#include "heRhoQGDThermo.H"
#include "pureMixture.H"
#include "powerLawTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */

/*
 *
 * Perfect gas EOS
 *
 */

makeThermo
(
    rhoQGDThermo,
    heRhoQGDThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoQGDThermo,
    heRhoQGDThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoQGDThermo,
    heRhoQGDThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoQGDThermo,
    heRhoQGDThermo,
    pureMixture,
    powerLawTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

/*
 *
 * const density gas EOS
 *
 */

makeThermo
(
    rhoQGDThermo,
    heRhoQGDThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    rhoConst,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
