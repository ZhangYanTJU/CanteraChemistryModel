/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeChemistryModel.H"

#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"
#include "CanteraChemistryModel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Cantera+++

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constEThermoPhysics
    );


    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constGasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        gasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        icoPoly8EThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidEThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constEThermoPhysics
    );


    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        psiReactionThermo,
        constHThermoPhysics
    );


    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constGasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        gasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        icoPoly8HThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constAdiabaticFluidHThermoPhysics
    );

    makeChemistryModelType
    (
        CanteraChemistryModel,
        rhoReactionThermo,
        constHThermoPhysics
    );
// Cantera---
}

// ************************************************************************* //
