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

Class
    Foam::CanteraChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    CanteraChemistryModelI.H
    CanteraChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef CanteraChemistryModel_H
#define CanteraChemistryModel_H

#include "StandardChemistryModel.H"
#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                      Class CanteraChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template<class ReactionThermo, class ThermoType>
class CanteraChemistryModel
:
    public StandardChemistryModel<ReactionThermo, ThermoType>
{
    // Private Data
    const word CanteraMechanismFile_;
    Cantera::IdealGasMix CanteraMixture_;
    //- Relative tolerance to control CVode
    scalar relTol_;
    //- Absolute tolerance to control CVode
    scalar absTol_;

    void setNumerics(Cantera::ReactorNet &r);

    // Private Member Functions

        //- Disallow copy constructor
        CanteraChemistryModel(const CanteraChemistryModel&);

        //- Disallow default bitwise assignment
        void operator=(const CanteraChemistryModel&);

        // override for Cantera
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("Cantera");


    // Constructors

        //- Construct from thermo
        CanteraChemistryModel(ReactionThermo& thermo);


    //- Destructor
    virtual ~CanteraChemistryModel();


    // Member Functions

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalar deltaT); //outer API-1

        //- Solve the reaction system for the given time step
        //  and return the characteristic time
        virtual scalar solve(const scalarField& deltaT); //outer API-2

        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const = 0;

        /*
        //- Return reaction rate of the speciei in reactionI
        virtual tmp<volScalarField::Internal> calculateRR
        (
            const label reactionI,
            const label speciei
        ) const;
        */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CanteraChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
