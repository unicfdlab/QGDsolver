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
    scalar ScQGD = 0.0, PrQGD = 1.0, uQGD;
//    
//  dict.lookup("ScQGD") >> ScQGD;
    dict.lookup("uQGD") >> uQGD;
    
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    uQGD_.primitiveFieldRef()  = uQGD;
    muQGD_.primitiveFieldRef() = 0.0;
    alphauQGD_.primitiveFieldRef() = 0.0;
    
    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    uQGD_.boundaryFieldRef()  = uQGD;
    muQGD_.boundaryFieldRef() = 0.0;
    alphauQGD_.boundaryFieldRef() = 0.0;
}

Foam::qgd::
constScPrModel1QGDCoeffs::~constScPrModel1QGDCoeffs()
{
}

void Foam::qgd::
constScPrModel1QGDCoeffs::correct(const Foam::rhoQGDThermo& qgdThermo)
{
    const volScalarField  nu     = qgdThermo.mu()/qgdThermo.rho();

    this->tauQGD_ = this->aQGD_ * this-> hQGD_/uQGD_;
}

//
//END-OF-FILE
//


