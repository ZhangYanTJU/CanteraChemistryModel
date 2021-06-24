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

\*---------------------------------------------------------------------------*/

#include "CanteraChemistryModel.H"
#include "UniformField.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "atomicWeights.H"
#include "chemistryReader.H"
#include "turbulentFluidThermoModel.H"
#include "physicoChemicalConstants.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::CanteraChemistryModel
(
    ReactionThermo& thermo
)
:
    StandardChemistryModel<ReactionThermo, ThermoType>(thermo),
    CanteraMechanismFile_(this->subDict("Cantera").lookup("CanteraMechanismFile_")),
    CanteraMixture_(CanteraMechanismFile_, ""),
    relTol_(this->subDict("Cantera").lookupOrDefault("relTol",1e-9)),
    absTol_(this->subDict("Cantera").lookupOrDefault("absTol",1e-15))
{
    Info<<"--- I am here in Cantera-construct ---"<<endl;
    Info<<"CanteraMechanismFile_ === "<<CanteraMechanismFile_<<endl;
    Info<<"relTol_ === "<<relTol_<<endl;
    Info<<"absTol_ === "<<absTol_<<endl;

    const long unsigned int nSpecies_OF = this->nSpecie_;
    if (CanteraMixture_.nSpecies() != nSpecies_OF)
    {
        FatalErrorIn("CanteraChemistryModel::solve")
            << "nSpecies in Cantera is not consistent with that in OpenFOAM" << endl
            << abort(FatalError) ;
    }
    else
    {
        for (label i=0; i<this->nSpecie_; i++)
        {
            if (CanteraMixture_.speciesName(i) != this->Y_[i].name())
            {
                FatalErrorIn("CanteraChemistryModel::solve")
                    << "species name in Cantera is not consistent with that in OpenFOAM" << endl
                    << "in Cantera:" << CanteraMixture_.speciesName(i) << endl
                    << "in OpenFOAM:" << this->Y_[i].name() << endl
                    << abort(FatalError) ;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::
~CanteraChemistryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalar deltaT
)
{
    // Don't allow the time-step to change more than a factor of 2
    return min
    (
        this->solve<UniformField<scalar>>(UniformField<scalar>(deltaT)),
        2*deltaT
    );
}


template<class ReactionThermo, class ThermoType>
Foam::scalar Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::solve
(
    const scalarField& deltaT
)
{
    return this->solve<scalarField>(deltaT);
}


template<class ReactionThermo, class ThermoType>
template<class DeltaTType>
Foam::scalar Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::solve
(
    const DeltaTType& deltaT
)
{
    Info<<"I am here in Cantera-solve === "<<endl;

    BasicChemistryModel<ReactionThermo>::correct();

    scalar deltaTMin = great;

    if (!this->chemistry_)
    {
        return deltaTMin;
    }

    Cantera::Reactor react;
    //Cantera::IdealGasReactor react;  // Constant-UV, default, constant Volumn
    //Cantera::IdealGasConstPressureReactor react;  // Constant-HP, constant pressure

    scalarField c0(this->nSpecie_);

    tmp<volScalarField> trho(this->thermo().rho());
    const scalarField& rho = trho();

    const scalarField& T = this->thermo().T();
    const scalarField& p = this->thermo().p();

    forAll(rho, cellI)
    {
        scalar Ti = T[cellI];
        const scalar rhoi = rho[cellI];
        scalar pi = p[cellI];
        try
        {
            scalarList yTemp_(this->nSpecie_);
            for (label i=0; i<this->nSpecie_; i++)
            {
                yTemp_[i] = this->Y_[i][cellI];
            }

            CanteraMixture_.setState_TPY(Ti, pi, yTemp_.begin());
            CanteraMixture_.getConcentrations(c0.begin());


            react.insert(CanteraMixture_);
            // useless in single mesh (since volume is 1 m^3 in single mesh and the default value is 1 too in Cantera)
            // not sure for normal case
            //react.setInitialVolume(1.0);
            react.setEnergy(0); // keep T const before and after sim.advance. this will give you a little improvement
            Cantera::ReactorNet sim;
            sim.addReactor(react);
            setNumerics(sim);


            sim.advance(deltaT[cellI]);


            CanteraMixture_.getConcentrations(this->c_.begin());
            for (label i=0; i<this->nSpecie_; i++)
            {
                this->RR_[i][cellI] =
                    (this->c_[i] - c0[i])*this->specieThermo_[i].W()/deltaT[cellI];
            }
        }
        catch(Cantera::CanteraError& err)
        {
            // handle exceptions thrown by Cantera
            std::cout << err.what() << std::endl;

            FatalErrorIn("CanteraChemistryModel::solve")
                << " Cantera complained in cell " << cellI
                << " with a Cantera::CanteraError"  << endl
                << abort(FatalError) ;
        }
    }

    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
void Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::setNumerics(Cantera::ReactorNet &sim)
{
    sim.setTolerances(relTol_,absTol_);
}

// ************************************************************************* //
