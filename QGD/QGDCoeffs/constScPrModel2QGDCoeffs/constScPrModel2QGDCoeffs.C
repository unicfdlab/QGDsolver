#include "constScPrModel2QGDCoeffs.H"
#include "psiQGDThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(constScPrModel2QGDCoeffs,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        constScPrModel2QGDCoeffs,
        dictionary
    );
}
}

Foam::qgd::
constScPrModel2QGDCoeffs::constScPrModel2QGDCoeffs
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar ScQGD = 1.0, PrQGD = 1.0;
    
    dict.lookup("ScQGD") >> ScQGD;
    dict.lookup("PrQGD") >> PrQGD;
    
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    
    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
}

Foam::qgd::
constScPrModel2QGDCoeffs::~constScPrModel2QGDCoeffs()
{
}

void Foam::qgd::
constScPrModel2QGDCoeffs::correct(const Foam::psiQGDThermo& qgdThermo)
{
    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    const volScalarField& mu     = qgdThermo.mu();
    
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
    
    //update tau with molecular viscosity
    this->tauQGD_ += mu / (p * this->ScQGD_);
}

//
//END-OF-FILE
//


