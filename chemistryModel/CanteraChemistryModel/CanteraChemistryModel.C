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

    ////Cantera::Reactor react;
    //Cantera::IdealGasReactor react;  // Constant-UV, default, constant Volumn
    Cantera::IdealGasConstPressureReactor react;  // Constant-HP, constant pressure

    scalarField c0(this->nSpecie_);
    scalarField c_ctr(this->nSpecie_);

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
                c0[i] = rhoi*this->Y_[i][cellI]/this->specieThermo_[i].W();
                yTemp_[i] = this->Y_[i][cellI];
            }

            CanteraMixture_.setState_TPY(Ti,pi,yTemp_.begin());


            CanteraMixture_.getConcentrations(c_ctr.begin());
            if ( mag(max(c_ctr-c0))> 0.1 )
            {
                FatalErrorIn("CanteraChemistryModel::solve")
                    << "Before ODE, c0 in OpenFOAM is not consistent with that in Cantera" << endl
                    << "in Cantera:" << c_ctr << endl
                    << "in OpenFOAM:" << c0 << endl
                    << abort(FatalError);
            }

            react.insert(CanteraMixture_);
            Cantera::ReactorNet sim;
            sim.addReactor(react);
            sim.setInitialTime(0);
            setNumerics(sim);

            // method-1:
            //scalar timeLeft = deltaT[cellI];
            //scalar old = 0;
            //while(timeLeft > small)
            //{
            //    sim.setMaxTimeStep(timeLeft);
            //    scalar now = sim.step();
            //    scalar dt = now-old;
            //    old = now;
            //    if (old >0) deltaTMin = min(dt, deltaTMin); // discard the first timestep, since it's too small
            //    //timeLeft = deltaT[cellI]-now; // this is also true
            //    timeLeft -= dt;
            //}



            // method-2:
            sim.advance(deltaT[cellI]);


            // method-3:
            //scalar timeLeft = deltaT[cellI];
            //scalar old = 0;
            //while(timeLeft > small)
            //{
            //    scalar now = sim.advance(deltaT[cellI], /* applylimit = */ true);
            //    scalar dt = now-old;
            //    old = now;
            //    if (old >0) deltaTMin = min(dt, deltaTMin); // discard the first timestep, since it's too small
            //    timeLeft -= dt;
            //}

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

    Info << "deltaTMin-ChemistryModel === " << deltaTMin << endl;
    return deltaTMin;
}


template<class ReactionThermo, class ThermoType>
void Foam::CanteraChemistryModel<ReactionThermo, ThermoType>::setNumerics(Cantera::ReactorNet &sim)
{
    sim.setTolerances(relTol_,absTol_);
}

// ************************************************************************* //
