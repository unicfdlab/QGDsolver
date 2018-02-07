/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "QGDCoeffs.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledFvsPatchFields.H"
#include "rhoQGDThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(QGDCoeffs, 0);
    defineRunTimeSelectionTable(QGDCoeffs, dictionary);
}
}

namespace Foam
{
namespace qgd
{

autoPtr<QGDCoeffs> QGDCoeffs::New
(
    const word& qgdCoeffsType,
    const fvMesh& mesh,
    const dictionary& dict
)
{
    Info<< "Selecting QGD coeffs evaluation approach type " << qgdCoeffsType << endl;
    
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(qgdCoeffsType);
    
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "QGDCoeffs::New(const word&, const fvMesh&)"
        )   << "Unknown QGD coeffs evaluation approach type " << qgdCoeffsType << nl << nl
        << "Valid model types are:" << nl
        << dictionaryConstructorTablePtr_->sortedToc()
        << exit(FatalError);
    }
    
    return autoPtr<QGDCoeffs>
    (
        cstrIter()
        (
            IOobject
            (
                qgdCoeffsType,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dict.subDict(qgdCoeffsType + "Dict")
        )
    );
}

//
QGDCoeffs::QGDCoeffs(const IOobject& io, const fvMesh& mesh, const dictionary& dict)
:
    regIOobject(io, false),
    refCount(),
    mesh_(mesh),
    runTime_(mesh_.time()),
    muQGD_
    (
        IOobject
        (
            "muQGD",
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
            "alphauQGD",
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
            "hQGD",
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
    tauQGD_
    (
        IOobject
        (
            "tauQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 1, 0, 0)
    ),
    PrQGD_
    (
        IOobject
        (
            "PrQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    ),
    ScQGD_
    (
        IOobject
        (
            "ScQGD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionSet(0, 0, 0, 0, 0)
    )
{
}

QGDCoeffs::~QGDCoeffs()
{
}

void Foam::qgd::QGDCoeffs::correct(const rhoQGDThermo& qgdThermo)
{
    forAll(tauQGD_, celli)
    {
        tauQGD_.primitiveFieldRef()[celli] = 0.0;
        muQGD_.primitiveFieldRef()[celli] = 0.0;
        alphauQGD_.primitiveFieldRef()[celli] = 0.0;
        ScQGD_.primitiveFieldRef()[celli] = 1.0;
        PrQGD_.primitiveFieldRef()[celli] = 1.0;
    }
    forAll(tauQGD_.boundaryField(), patchi)
    {
        forAll(tauQGD_.boundaryField()[patchi], facei)
        {
            tauQGD_.boundaryFieldRef()[patchi][facei] = 
                0.0;
            muQGD_.boundaryFieldRef()[patchi][facei] = 
                0.0;
            alphauQGD_.boundaryFieldRef()[patchi][facei] = 
                0.0;
            PrQGD_.boundaryFieldRef()[patchi][facei] = 
                1.0;
            ScQGD_.boundaryFieldRef()[patchi][facei] = 
                1.0;
        }
    }
}

Foam::tmp<Foam::volScalarField> Foam::qgd::QGDCoeffs::qgdLength(const fvMesh& mesh)
{
    tmp<volScalarField> tqgdLen
    (
        new volScalarField
        (
            IOobject
            (
                "hQGD",
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

const Foam::volScalarField& Foam::qgd::QGDCoeffs::hQGD() const
{
    return hQGD_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::alphauQGD() const
{
    return alphauQGD_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::muQGD() const
{
    return muQGD_;
}

const Foam::volScalarField& Foam::qgd::QGDCoeffs::tauQGD() const
{
    return tauQGD_;
}

}; //namespace qgd

}; //namespace Foam


//END-OF-FILE

