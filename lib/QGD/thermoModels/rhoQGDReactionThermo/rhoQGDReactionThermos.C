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

Group
    grpPsiQGDReactionThermo

\*---------------------------------------------------------------------------*/

#include "makeReactionThermo.H"

#include "rhoQGDReactionThermo.H"
#include "heRhoQGDThermo.H"

#include "specie.H"
#include "perfectFluid.H"
#include "perfectGas.H"
#include "hConstThermo.H"
#include "janafThermo.H"
#include "thermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"
#include "singleStepReactingMixture.H"
#include "singleComponentMixture.H"

#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// perfectFluid

makeReactionThermos
(
    rhoQGDThermo,
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    reactingMixture,
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectFluid,
    specie
);

// perfectGas

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Multi-component thermo for internal energy

makeThermoPhysicsReactionThermos
(
    rhoQGDThermo,
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    multiComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoQGDThermo,
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    multiComponentMixture,
    gasEThermoPhysics
);


// Reaction thermo for internal energy

makeThermoPhysicsReactionThermos
(
    rhoQGDThermo,
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    reactingMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermos
(
    rhoQGDThermo,
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    reactingMixture,
    gasEThermoPhysics
);


// Single-step reaction thermo for internal energy

makeThermoPhysicsReactionThermos
(
    rhoQGDThermo,
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    singleStepReactingMixture,
    gasEThermoPhysics
);


// Single-component thermo for sensible enthalpy

makeThermoPhysicsReactionThermo
(
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    singleComponentMixture,
    constGasHThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    singleComponentMixture,
    gasHThermoPhysics
);


// Single-component thermo for internal energy

makeThermoPhysicsReactionThermo
(
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    singleComponentMixture,
    constGasEThermoPhysics
);

makeThermoPhysicsReactionThermo
(
    rhoQGDReactionThermo,
    heRhoQGDThermo,
    singleComponentMixture,
    gasEThermoPhysics
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
