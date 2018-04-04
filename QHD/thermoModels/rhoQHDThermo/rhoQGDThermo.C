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

Class
    Foam::qgd::rhoQGDThermo



Description
    Class describing thermophysical properties of perfect gas with motion
    governed by Quasi-Hydro dynamics equations.


SourceFiles

    rhoQGDThermo.C



\*---------------------------------------------------------------------------*/

#include "rhoQGDThermo.H"
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
    defineTypeNameAndDebug(rhoQGDThermo, 0);
    defineRunTimeSelectionTable(rhoQGDThermo, fvMesh);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rhoQGDThermo::rhoQGDThermo(const fvMesh& mesh, const word& phaseName)
:
    rhoThermo(mesh, phaseName),
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

Foam::autoPtr<Foam::rhoQGDThermo> Foam::rhoQGDThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{   
    return basicThermo::New<rhoQGDThermo>(mesh, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rhoQGDThermo::~rhoQGDThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rhoQGDThermo::correctQGD()
{
    this->gamma_ == (this->Cp() / this->Cv());
    c_ = sqrt(gamma_ / this->psi());
    
    
    qgdCoeffsPtr_->correct(*this);
    const volScalarField& muQGD = this->muQGD();
    const volScalarField& alphauQGD = this->alphauQGD();
    
//    this->tauQGD_ = this->aQGD_ * this->hQGD_  / this->c_;
    
    forAll(mu_.primitiveField(), celli)
    {
//        muQGD_.primitiveFieldRef()[celli] = 
//            p_.primitiveField()[celli] * 
//            ScQGD_ *
//            tauQGD_.primitiveField()[celli];
//            //aQGD_.primitiveField()[celli] *
//            //hQGD_.primitiveField()[celli] /
//            //c_.primitiveField()[celli];
//        
//        alphauQGD_.primitiveFieldRef()[celli] = muQGD_.primitiveField()[celli] / PrQGD_;
//        
        mu_.primitiveFieldRef()[celli] +=
            muQGD.primitiveField()[celli];
        alpha_.primitiveFieldRef()[celli] +=
            alphauQGD.primitiveField()[celli];
   }
    
    forAll(mu_.boundaryField(), patchi)
    {
        forAll(p_.boundaryField()[patchi], facei)
        {
//            muQGD_.boundaryFieldRef()[patchi][facei] = 
//                p_.boundaryField()[patchi][facei] * ScQGD_ *
//                tauQGD_.boundaryField()[patchi][facei];
//                //aQGD_.boundaryField()[patchi][facei] *
//                //hQGD_.boundaryField()[patchi][facei] /
//                //c_.boundaryField()[patchi][facei];
//
//            alphauQGD_.boundaryFieldRef()[patchi][facei] = 
//                muQGD_.boundaryField()[patchi][facei] / PrQGD_;
        }
        
        mu_.boundaryFieldRef()[patchi] +=
            muQGD.boundaryField()[patchi];
        alpha_.boundaryFieldRef()[patchi] +=
            alphauQGD.boundaryField()[patchi];
   }
}

bool Foam::rhoQGDThermo::read()
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

const Foam::volScalarField& Foam::rhoQGDThermo::c() const
{
    return this->c_;
}

const Foam::volScalarField& Foam::rhoQGDThermo::hQGD() const
{
    return qgdCoeffsPtr_->hQGD();
}

const Foam::volScalarField& Foam::rhoQGDThermo::tauQGD() const
{
    return qgdCoeffsPtr_->tauQGD();
}

const Foam::volScalarField& Foam::rhoQGDThermo::muQGD() const
{
    return qgdCoeffsPtr_->muQGD();
}

const Foam::volScalarField& Foam::rhoQGDThermo::alphauQGD() const
{
    return qgdCoeffsPtr_->alphauQGD();
}


Foam::Switch Foam::rhoQGDThermo::implicitDiffusion() const
{
    return implicitDiffusion_;
}

// ************************************************************************* //
