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
    muQGD_
    (
        IOobject
        (
            phasePropertyName("thermo:muQGD"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),
    alphauQGD_
    (
        IOobject
        (
            phasePropertyName("thermo:alphauQGD"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(1, -1, -1, 0, 0)
    ),
    hQGD_
    (
        IOobject
        (
            phasePropertyName("thermo:hQGD"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        qgdLength(mesh)()
    ),
    aQGD_
    (
        IOobject
        (
            "alphaQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    gamma_
    (
        "thermo:gamma",
        aQGD_*1.0
        //this->Cp() / this->Cv()
    ),
    c_
    (
        "thermo:c",
        hQGD_/mesh.time().deltaT()
    ),
    tauQGD_
    (
        IOobject
        (
            phasePropertyName("thermo:tauQGD"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 1, 0, 0)
    ),
    PrQGD_(1.0),
    ScQGD_(0.0),
    implicitDiffusion_(false)//,
    //zeroWallQGDFlux_(false)
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

Foam::tmp<Foam::volScalarField> Foam::psiQGDThermo::qgdLength(const fvMesh& mesh)
{
    tmp<volScalarField> tqgdLen
    (
        new volScalarField
        (
            IOobject
            (
                "thermo:hQGD",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionSet(0, 1, 0, 0, 0)
        )
    );
    
    volScalarField& qgdLen = tqgdLen.ref();
    
    surfaceScalarField hf
    (
       1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs())
    );
    
    scalar cdist = 0.0;
    forAll(qgdLen, celli)
    {
        const cell& c = mesh.cells()[celli];
        scalar maxlen = 0.0;
        forAll(c, facei)
        {
            if (mesh.isInternalFace(c[facei]))
            {
                cdist = hf[c[facei]];
                if (cdist > maxlen)
                {
                    maxlen = hf[c[facei]];
                }
            }
        }
        qgdLen.primitiveFieldRef()[celli] = maxlen;
    }
    
    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& fvp = mesh.boundary()[patchi];
        if (fvp.coupled())
        {
            qgdLen.boundaryFieldRef()[patchi] = hf.boundaryField()[patchi];
        }
        else
        {
            qgdLen.boundaryFieldRef()[patchi] = 2.0*hf.boundaryField()[patchi];
        }
    }
    
    return tqgdLen;
}

void Foam::psiQGDThermo::correctQGD()
{
    this->gamma_ == (this->Cp() / this->Cv());
    c_ = sqrt(gamma_ / this->psi());
    this->tauQGD_ = this->aQGD_ * this->hQGD_  / this->c_;
    
    forAll(p_.primitiveField(), celli)
    {
        muQGD_.primitiveFieldRef()[celli] = 
            p_.primitiveField()[celli] * 
            ScQGD_ *
            tauQGD_.primitiveField()[celli];
            //aQGD_.primitiveField()[celli] *
            //hQGD_.primitiveField()[celli] /
            //c_.primitiveField()[celli];
        
        alphauQGD_.primitiveFieldRef()[celli] = muQGD_.primitiveField()[celli] / PrQGD_;
        
        mu_.primitiveFieldRef()[celli] +=
            muQGD_.primitiveField()[celli];
        alpha_.primitiveFieldRef()[celli] +=
            alphauQGD_.primitiveField()[celli];
    }
    
    forAll(p_.boundaryField(), patchi)
    {
        forAll(p_.boundaryField()[patchi], facei)
        {
            muQGD_.boundaryFieldRef()[patchi][facei] = 
                p_.boundaryField()[patchi][facei] * ScQGD_ *
                tauQGD_.boundaryField()[patchi][facei];
                //aQGD_.boundaryField()[patchi][facei] *
                //hQGD_.boundaryField()[patchi][facei] /
                //c_.boundaryField()[patchi][facei];

            alphauQGD_.boundaryFieldRef()[patchi][facei] = 
                muQGD_.boundaryField()[patchi][facei] / PrQGD_;
        }
        
        mu_.boundaryFieldRef()[patchi] +=
            muQGD_.boundaryField()[patchi];
        alpha_.boundaryFieldRef()[patchi] +=
            alphauQGD_.boundaryField()[patchi];
    }
    
    
    //this->tauQGD_ = this->mu_ / (this->p_ * this->ScQGD_);
    
    //remove QGD viscosity at walls
//    Info << "zeroWallQGDFlux_ = " << zeroWallQGDFlux_ << endl;
//    if (zeroWallQGDFlux_)
//    {
//        Info << "removing qgd viscosity" << endl;
//        const fvMesh& mesh = p_.mesh();
//        forAll(mesh.boundary(), iPatch)
//        {
//            if (isA<wallFvPatch>(mesh.boundary()[iPatch]))
//            {
//                Info << "removing qgd viscosity" << endl;
//                mu_.boundaryFieldRef()[iPatch] -=
//                    muQGD_.boundaryField()[iPatch];
//                alpha_.boundaryFieldRef()[iPatch] -=
//                    alphauQGD_.boundaryField()[iPatch];
//            }
//        }
//    }
}

bool Foam::psiQGDThermo::read()
{
    if (!basicThermo::read())
    {
        return false;
    }
    
    this->subDict("QGD").lookup("ScQGD") >> ScQGD_;
    this->subDict("QGD").lookup("PrQGD") >> PrQGD_;
    this->subDict("QGD").lookup("implicitDiffusion") >> implicitDiffusion_;
    //if (this->subDict("QGD").found("zeroWallQGDFlux"))
    //{
    //    this->subDict("QGD").lookup("zeroWallQGDFlux") >> zeroWallQGDFlux_;
    //}
    
    return true;
}

const Foam::volScalarField& Foam::psiQGDThermo::c() const
{
    return c_;
}

const Foam::volScalarField& Foam::psiQGDThermo::hQGD() const
{
    return hQGD_;
}

const Foam::volScalarField& Foam::psiQGDThermo::tauQGD() const
{
    return tauQGD_;
}

const Foam::volScalarField& Foam::psiQGDThermo::muQGD() const
{
    return muQGD_;
}

const Foam::volScalarField& Foam::psiQGDThermo::alphauQGD() const
{
    return alphauQGD_;
}

Foam::dimensionedScalar Foam::psiQGDThermo::ScQGD() const
{
    return dimensionedScalar("ScQGD", dimless, ScQGD_);
}

Foam::Switch Foam::psiQGDThermo::implicitDiffusion() const
{
    return implicitDiffusion_;
}

// ************************************************************************* //
