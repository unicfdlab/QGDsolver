#include "alphaHbyUQHD.H"
#include "rhoQGDThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(alphaHbyUQHD,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        alphaHbyUQHD,
        dictionary
    );
}
}

Foam::qgd::
alphaHbyUQHD::alphaHbyUQHD
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict),
    uQGD_(0.0)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;
    
//    
    dict.lookup("uQGD") >> uQGD_;
    
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    muQGD_.primitiveFieldRef() = 0.0;
    alphauQGD_.primitiveFieldRef() = 0.0;
    
    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    muQGD_.boundaryFieldRef() = 0.0;
    alphauQGD_.boundaryFieldRef() = 0.0;
}

Foam::qgd::
alphaHbyUQHD::~alphaHbyUQHD()
{
}

void Foam::qgd::
alphaHbyUQHD::correct(const Foam::rhoQGDThermo& qgdThermo)
{
    this->tauQGD_ = this->aQGD_ * this->hQGD_ /dimensionedScalar("uQGD", dimLength/dimTime, uQGD_);
}

//
//END-OF-FILE
//


