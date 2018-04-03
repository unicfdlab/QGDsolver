#include "alphaH2bynuQHD.H"
#include "rhoQGDThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(alphaH2bynuQHD,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        alphaH2bynuQHD,
        dictionary
    );
}
}

Foam::qgd::
alphaH2bynuQHD::alphaH2bynuQHD
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;
//    
    
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
alphaH2bynuQHD::~alphaH2bynuQHD()
{
}

void Foam::qgd::
alphaH2bynuQHD::correct(const Foam::rhoQGDThermo& qgdThermo)
{
    const volScalarField  nu     = qgdThermo.mu()/qgdThermo.rho();
    this->tauQGD_ = this->aQGD_ * sqr(this->hQGD_) / nu;
}

//
//END-OF-FILE
//


