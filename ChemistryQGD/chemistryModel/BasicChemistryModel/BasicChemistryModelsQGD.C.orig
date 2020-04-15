/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics
    for QGD thermo

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiQGDReactionThermo.H"
#include "rhoReactionThermo.H"

#include "StandardChemistryModel.H"
#include "TDACChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Make base types
    makeChemistryModel(psiQGDReactionThermo);
    
    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        psiQGDReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        TDACChemistryModel,
        psiQGDReactionThermo,
        constHThermoPhysics
    );



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
        StandardChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        StandardChemistryModel,
        rhoReactionThermo,
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

// ************************************************************************* //
