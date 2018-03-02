#include "constScPrModel1QGDCoeffs.H"
#include "rhoQGDThermo.H"
#include "addToRunTimeSelectionTable.H"

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
    Info << "Done" << endl;
    scalar ScQGD = 1.0, PrQGD = 1.0;
//    
    dict.lookup("ScQGD") >> ScQGD;
    dict.lookup("PrQGD") >> PrQGD;
    
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    
    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
}

Foam::qgd::
constScPrModel1QGDCoeffs::~constScPrModel1QGDCoeffs()
{
}

void Foam::qgd::
constScPrModel1QGDCoeffs::correct(const Foam::rhoQGDThermo& qgdThermo)
{
    Info << "Done" << endl;

    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    
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

}

//
//END-OF-FILE
//


