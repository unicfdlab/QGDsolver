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
    grpRhoQGDThermo

\*---------------------------------------------------------------------------*/

#include "heRhoQGDThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::calculate()
{
    if (this->correctPsiOnly_)
    {
        calculatePsi();
        return;
    }

    const scalarField& hCells = this->he();
    const scalarField& pCells = this->p_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        TCells[celli] = mixture_.THE
        (
            hCells[celli],
            pCells[celli],
            TCells[celli]
        );

        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
        muCells[celli] = mixture_.mu(pCells[celli], TCells[celli]);
        alphaCells[celli] = mixture_.alphah(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
	this->rho_.boundaryFieldRef();

    volScalarField::Boundary& heBf =
        this->he().boundaryFieldRef();

    volScalarField::Boundary& muBf =
        this->mu_.boundaryFieldRef();

    volScalarField::Boundary& alphaBf =
        this->alpha_.boundaryFieldRef();


    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
	fvPatchScalarField& prho = rhoBf[patchi];
        fvPatchScalarField& phe = heBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pT.fixesValue())
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);
                phe[facei] = mixture_.HE(pp[facei], pT[facei]);
		prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
        else
        {
            forAll(pT, facei)
            {
                const typename MixtureType::thermoType& mixture_ =
                    this->patchFaceMixture(patchi, facei);

                pT[facei] = mixture_.THE(phe[facei], pp[facei], pT[facei]);
		prho[facei] = mixture_.rho(pp[facei], pT[facei]);
                pmu[facei] = mixture_.mu(pp[facei], pT[facei]);
                palpha[facei] = mixture_.alphah(pp[facei], pT[facei]);
            }
        }
    }

    calculatePsi();

    if (!this->isochoric())
    {
        this->gamma_ == (this->Cp() / this->Cv());
        this->c_ = sqrt(this->gamma_ / this->psi());
    }
    this->correctQGD(this->mu_, this->alpha_);
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::calculatePsi()
{
    if (this->multicomponent_)
    {
        multicomponentCalculatePsi();
        return;
    }
    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& pCells = this->p_;

    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();

    forAll(TCells, celli)
    {
        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);

        psiCells[celli] = mixture_.psi(pCells[celli], TCells[celli]);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];
    
        forAll(pT, facei)
        {
            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);
    
            ppsi[facei] = mixture_.psi(pp[facei], pT[facei]);
            prho[facei] = mixture_.rho(pp[facei], pT[facei]);
        }
    }
}

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::multicomponentCalculatePsi()
{
    const scalar R1 = 7000.0;
    const scalar R2 = 287.0;

    const scalar rho01 = 99.9; //.0;
    const scalar rho02 = 0.0;

    const volScalarField &Y1 = this->T_.mesh().thisDb().template lookupObject<volScalarField>("water");
    const volScalarField &Y2 = this->T_.mesh().thisDb().template lookupObject<volScalarField>("air");

    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& pCells = this->p_;
    const scalarField& Y1Cells= Y1.primitiveField();
    const scalarField& Y2Cells= Y2.primitiveField();

    scalarField& psiCells = this->psi_.primitiveFieldRef();
    scalarField& rhoCells = this->rho_.primitiveFieldRef();

    scalar psi1 = 0.0;
    scalar psi2 = 0.0;
    scalar rho1 = 0.0;
    scalar rho2 = 0.0;
    scalar y1   = 0.0;
    scalar y2   = 0.0;
    scalar y1byrho1 = 0.0;
    scalar y2byrho2 = 0.0;

    forAll(TCells, celli)
    {
        psi1 = 1.0 / (R1 *TCells[celli]);
        psi2 = 1.0 / (R2 *TCells[celli]);
        rho1 = psi1*pCells[celli] + rho01;
        rho2 = psi2*pCells[celli] + rho02;
        y1   = Y1Cells[celli];
        y2   = Y2Cells[celli];
        y1byrho1 = y1 / rho1;
        y2byrho2 = y2 / rho2;
        psiCells[celli] = pow(y1byrho1 + y2byrho2,-2.0)*(psi1*y1byrho1/rho1+psi2*y2byrho2/rho2);

        const typename MixtureType::thermoType& mixture_ =
            this->cellMixture(celli);
        rhoCells[celli] = mixture_.rho(pCells[celli], TCells[celli]);
    }

    volScalarField::Boundary& pBf =
        this->p_.boundaryFieldRef();

    volScalarField::Boundary& TBf =
        this->T_.boundaryFieldRef();

    volScalarField::Boundary& psiBf =
        this->psi_.boundaryFieldRef();

    const volScalarField::Boundary& Y1Bf =
        Y1.boundaryField();

    const volScalarField::Boundary& Y2Bf =
        Y2.boundaryField();

    volScalarField::Boundary& rhoBf =
        this->rho_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pT = TBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        const fvPatchScalarField& pY1 = Y1Bf[patchi];
        const fvPatchScalarField& pY2 = Y2Bf[patchi];
        fvPatchScalarField& prho = rhoBf[patchi];

        forAll(pT, facei)
        {
            psi1 = 1.0 / (R1 *pT[facei]);
            psi2 = 1.0 / (R2 *pT[facei]);
            rho1 = psi1*pp[facei] + rho01;
            rho2 = psi2*pp[facei] + rho02;
            y1byrho1 = y1 / rho1;
            y2byrho2 = y2 / rho2;
            y1   = pY1[facei];
            y2   = pY2[facei];
            ppsi[facei] = pow(y1byrho1 + y2byrho2,-2.0)*(psi1*y1byrho1/rho1+psi2*y2byrho2/rho2);

            const typename MixtureType::thermoType& mixture_ =
                this->patchFaceMixture(patchi, facei);
            prho[facei] = mixture_.rho(pp[facei], pT[facei]);
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::heRhoQGDThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName)
{
    calculate();
    // Switch on saving old time
    this->psi_.oldTime();
}

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::heRhoQGDThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    heThermo<BasicPsiThermo, MixtureType>(mesh, phaseName, dictionaryName)
{
    calculate();
    // Switch on saving old time
    this->psi_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::~heRhoQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicPsiThermo, class MixtureType>
void Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::correct()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    // force the saving of the old-time values
    this->psi_.oldTime();

    calculate();

    if (debug)
    {
        Info<< "    Finished" << endl;
    }
}

template<class BasicPsiThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::heRhoQGDThermo<BasicPsiThermo, MixtureType>::gamma() const
{
    return this->gamma_;
}


// ************************************************************************* //
