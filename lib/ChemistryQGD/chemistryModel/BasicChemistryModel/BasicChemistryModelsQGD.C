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
    This file is part of QGDsolver library, based on OpenFOAM+.

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
    grpCheQGDModels

Description
    Creates chemistry model instances templated on the type of thermodynamics
    for QGD thermo

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"
#include "psiQGDReactionThermo.H"
#include "rhoQGDReactionThermo.H"
#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    makeChemistryModel(psiQGDReactionThermo);

    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constEThermoPhysics
    );

}

namespace Foam
{
    // Make base types
    makeChemistryModel(rhoQGDReactionThermo);

    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoQGDReactionThermo,
        constEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        rhoQGDReactionThermo,
        constEThermoPhysics
    );

}

// ************************************************************************* //
