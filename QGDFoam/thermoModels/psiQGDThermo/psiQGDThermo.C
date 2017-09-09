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
    alphahQGD_
    (
        IOobject
        (
            phasePropertyName("thermo:alphahQGD"),
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
        "thermo:hQGD",
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
    ScQGD_(0.0)
{
    this->read();
//    this->gamma_ = this->Cp() / this->Cv();
//    this->c_     = sqrt(gamma_ / this->psi());
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
    
    forAll(qgdLen, celli)
    {
        const cell& c = mesh.cells()[celli];
        scalar maxlen = 0.0;
        forAll(c, facei)
        {
            if (mesh.isInternalFace(c[facei]))
            {
                if (hf[c[facei]] > maxlen)
                {
                    maxlen = hf[c[facei]];
                }
            }
        }
        qgdLen[celli] = maxlen;
    }
    
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (!isA<coupledFvPatch>(mesh.boundaryMesh()[patchi]))
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

    forAll(p_.primitiveField(), celli)
    {
        muQGD_.primitiveFieldRef()[celli] = 
            p_.primitiveField()[celli] * 
            ScQGD_*
            aQGD_.primitiveField()[celli] *
            hQGD_.primitiveField()[celli] /
            c_.primitiveField()[celli];
        
        alphahQGD_.primitiveFieldRef()[celli] = muQGD_.primitiveField()[celli] / PrQGD_;
        
        mu_.primitiveFieldRef()[celli] +=
            muQGD_.primitiveField()[celli];
        alpha_.primitiveFieldRef()[celli] +=
            alphahQGD_.primitiveField()[celli];
    }
    
    forAll(p_.boundaryField(), patchi)
    {
        forAll(p_.boundaryField()[patchi], facei)
        {
            muQGD_.boundaryFieldRef()[patchi][facei] = 
                p_.boundaryField()[patchi][facei] * ScQGD_ *
                aQGD_.boundaryField()[patchi][facei] *
                hQGD_.boundaryField()[patchi][facei] /
                c_.boundaryField()[patchi][facei];

            alphahQGD_.boundaryFieldRef()[patchi][facei] = 
                muQGD_.boundaryField()[patchi][facei] / PrQGD_;
        }

        mu_.boundaryFieldRef()[patchi] +=
            muQGD_.boundaryField()[patchi];
        alpha_.boundaryFieldRef()[patchi] +=
            alphahQGD_.boundaryField()[patchi];
    }
    
    this->tauQGD_ = this->mu_ / (this->p_ * this->ScQGD_);
}

bool Foam::psiQGDThermo::read()
{
    if (!basicThermo::read())
    {
        return false;
    }
    
    this->subDict("QGD").lookup("ScQGD") >> ScQGD_;
    this->subDict("QGD").lookup("PrQGD") >> PrQGD_;
    
    return true;
}

Foam::tmp<Foam::volScalarField> Foam::psiQGDThermo::c() const
{
    return c_;
}

Foam::tmp<Foam::volScalarField> Foam::psiQGDThermo::hQGD() const
{
    return hQGD_;
}

Foam::tmp<Foam::volScalarField> Foam::psiQGDThermo::tauQGD() const
{
    return tauQGD_;
}

Foam::dimensionedScalar Foam::psiQGDThermo::ScQGD() const
{
    return dimensionedScalar("ScQGD", dimless, ScQGD_);
}

// ************************************************************************* //
