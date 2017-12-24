/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "psiQGDThermo.H"
#include "surfaceFields.H"
#include "typeInfo.H"
#include "coupledPolyPatch.H"
#include "coupledFvPatch.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "symmetryPlaneFvPatch.H"
#include "symmetryFvPatch.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(psiQGDThermo, 0);
    defineRunTimeSelectionTable(psiQGDThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::psiQGDThermo::psiQGDThermo(const fvMesh& mesh, const word& phaseName)
:
    psiThermo(mesh, phaseName),
    qgdCoeffsPtr_
    (
        Foam::qgd::QGDCoeffs::New
        (
            this->subDict("QGD").lookup("QGDCoeffs"),
            mesh,
            this->subDict("QGD")
        )
    ),
    c_
    (
        "thermo:c",
        qgdCoeffsPtr_->hQGD()/mesh.time().deltaT()
    ),
    gamma_
    (
        "thermo:gamma",
        c_ / c_
        //this->Cp() / this->Cv()
    ),
    implicitDiffusion_(false)
{
    this->read();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiQGDThermo> Foam::psiQGDThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<psiQGDThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::psiQGDThermo::~psiQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::psiQGDThermo::correctQGD()
{
    this->gamma_ == (this->Cp() / this->Cv());
    c_ = sqrt(gamma_ / this->psi());
    
    qgdCoeffsPtr_->correct(*this);
    
    const volScalarField& muQGD = this->muQGD();
    const volScalarField& alphauQGD = this->alphauQGD();
    
    forAll(mu_.primitiveField(), celli)
    {
        mu_.primitiveFieldRef()[celli] +=
            muQGD.primitiveField()[celli];
        
        alpha_.primitiveFieldRef()[celli] +=
            alphauQGD.primitiveField()[celli];
    }
    
    forAll(mu_.boundaryField(), patchi)
    {
        forAll(p_.boundaryField()[patchi], facei)
        {
            mu_.boundaryFieldRef()[patchi][facei] +=
                muQGD.boundaryField()[patchi][facei];
            
            alpha_.boundaryFieldRef()[patchi][facei] +=
                alphauQGD.boundaryField()[patchi][facei];
        }
    }
}

bool Foam::psiQGDThermo::read()
{
    if (!basicThermo::read())
    {
        return false;
    }
    
    //this->subDict("QGD").lookup("ScQGD") >> ScQGD_;
    //this->subDict("QGD").lookup("PrQGD") >> PrQGD_;
    this->subDict("QGD").lookup("implicitDiffusion") >> implicitDiffusion_;
    
    return true;
}

const Foam::volScalarField& Foam::psiQGDThermo::c() const
{
    return this->c_;
}

const Foam::volScalarField& Foam::psiQGDThermo::hQGD() const
{
    return qgdCoeffsPtr_->hQGD();
}

const Foam::volScalarField& Foam::psiQGDThermo::tauQGD() const
{
    return qgdCoeffsPtr_->tauQGD();
}

const Foam::volScalarField& Foam::psiQGDThermo::muQGD() const
{
    return qgdCoeffsPtr_->muQGD();
}

const Foam::volScalarField& Foam::psiQGDThermo::alphauQGD() const
{
    return qgdCoeffsPtr_->alphauQGD();
}

Foam::Switch Foam::psiQGDThermo::implicitDiffusion() const
{
    return implicitDiffusion_;
}

// ************************************************************************* //
