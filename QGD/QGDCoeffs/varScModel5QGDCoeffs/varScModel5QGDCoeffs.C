#include "varScModel5QGDCoeffs.H"
#include "psiQGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSmooth.H"
#include "fvcGrad.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(varScModel5QGDCoeffs,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        varScModel5QGDCoeffs,
        dictionary
    );
}
}

Foam::qgd::
varScModel5QGDCoeffs::varScModel5QGDCoeffs
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict),
    smoothCoeff_(0.02)
{
    scalar ScQGD = 1.0, PrQGD = 1.0;
    
    dict.lookup("ScQGD") >> ScQGD;
    dict.lookup("PrQGD") >> PrQGD;
    
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    
    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    
    if (dict.found("smoothCoeff"))
    {
        dict.lookup("smoothCoeff") >> smoothCoeff_;
    }
    
    fvc::smooth(this->hQGD_, smoothCoeff_);
}

Foam::qgd::
varScModel5QGDCoeffs::~varScModel5QGDCoeffs()
{
}

void Foam::qgd::
varScModel5QGDCoeffs::correct(const Foam::psiQGDThermo& qgdThermo)
{
    const volScalarField& cSound = qgdThermo.c();
    const volScalarField& p      = qgdThermo.p();
    
    this->tauQGD_ = this->aQGD_ * this->hQGD_  / cSound;
    
    scalar rC = 0.1;
    
    this->ScQGD_ = 
        rC *
        (mag(fvc::grad(p)) * hQGD_ / p) +
        (1.0 - rC) * ScQGD_;
    
    this->ScQGD_ = 
        max(this->ScQGD_, 0.05);
    
    if (runTime_.outputTime())
    {
        this->ScQGD_.write();
    }
    
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


