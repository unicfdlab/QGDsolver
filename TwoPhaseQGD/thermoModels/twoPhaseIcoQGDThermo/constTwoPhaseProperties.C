#include "constTwoPhaseProperties.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constTwoPhaseProperties::constTwoPhaseProperties
(
    const dictionary& dict,
    const word& phase1name,
    const word& phase2name
)
:
    nu1_("nu"+phase1name,dict),
    nu2_("nu"+phase2name,dict),
    rho1_("rho"+phase1name,dict),
    rho2_("rho"+phase2name,dict),
    Tau1_("tau"+phase1name,dict),
    Tau2_("tau"+phase2name,dict)
{

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constTwoPhaseProperties::~constTwoPhaseProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::nu1() const
{
    return nu1_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::nu2() const
{
    return nu2_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::rho1() const
{
    return rho1_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::rho2() const
{
    return rho2_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::Tau1() const
{
    return Tau1_;
}

const Foam::dimensionedScalar& Foam::constTwoPhaseProperties::Tau2() const
{
    return Tau2_;
}


// ************************************************************************* //
