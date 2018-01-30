#include "constScPrModel1QGDCoeffs.H"
#include "psiQGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(constScPrModel1QGDCoeffs,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        constScPrModel1QGDCoeffs,
        dictionary
    );
}
}

Foam::qgd::
constScPrModel1QGDCoeffs::constScPrModel1QGDCoeffs
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar PrQGD = 1.0;
    
    if (dict.found("PrQGD"))
    {
        dict.lookup("PrQGD") >> PrQGD;
    }
    
    PrQGD_.primitiveFieldRef() = PrQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    
    IOobject ScHeader
    (
        "ScQGD",
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    
    if (ScHeader.headerOk())
    {
        //do nothing
    }
    else
    {
        scalar ScQGD = 1.0;
        dict.lookup("ScQGD") >> ScQGD;
        ScQGD_.primitiveFieldRef() = ScQGD;
        ScQGD_.boundaryFieldRef() = ScQGD;
    }
}

Foam::qgd::
constScPrModel1QGDCoeffs::~constScPrModel1QGDCoeffs()
{
}

void Foam::qgd::
constScPrModel1QGDCoeffs::correct(const Foam::psiQGDThermo& qgdThermo)
{
    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    
    this->tauQGDf_= linearInterpolate(this->aQGD_ / cSound) * hQGDf_;
    this->tauQGD_ = this->aQGD_ * this->hQGD_  / cSound;
    
    forAll(p.primitiveField(), celli)
    {
        muQGD_.primitiveFieldRef()[celli] = 
            p.primitiveField()[celli] * 
            ScQGD_.primitiveField()[celli] *
            tauQGD_.primitiveField()[celli];
        
        alphauQGD_.primitiveFieldRef()[celli] = muQGD_.primitiveField()[celli] / 
            PrQGD_.primitiveField()[celli];
    }
    
    forAll(p.boundaryField(), patchi)
    {
        forAll(p.boundaryField()[patchi], facei)
        {
            muQGD_.boundaryFieldRef()[patchi][facei] = 
                p.boundaryField()[patchi][facei] * 
                ScQGD_.boundaryField()[patchi][facei] *
                tauQGD_.boundaryField()[patchi][facei];

            alphauQGD_.boundaryFieldRef()[patchi][facei] = 
                muQGD_.boundaryFieldRef()[patchi][facei] /
                PrQGD_.boundaryField()[patchi][facei];
        }
    }
    
    if (runTime_.outputTime())
    {
        ScQGD_.write();
    }
}

//
//END-OF-FILE
//


