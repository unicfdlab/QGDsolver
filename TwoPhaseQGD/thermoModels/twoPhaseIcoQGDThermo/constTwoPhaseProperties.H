#ifndef constTwoPhaseProperties_H
#define constTwoPhaseProperties_H

#include "dictionary.H"
#include "dimensionedType.H"
#include "scalar.H"
#include "word.H"

namespace Foam
{

class constTwoPhaseProperties
{

    //-
    dimensionedScalar nu1_;
    
    //-
    dimensionedScalar nu2_;
    
    //-
    dimensionedScalar rho1_;
    
    //-
    dimensionedScalar rho2_;
    
    //-
    dimensionedScalar Tau1_;
    
    //-
    dimensionedScalar Tau2_;

public:

    //-
    constTwoPhaseProperties(const dictionary& dict, const word& phase1name, const word& phase2name);
    
    //-
    ~constTwoPhaseProperties();

    //-
    const dimensionedScalar& nu1() const;
    
    //-
    const dimensionedScalar& nu2() const;
    
    //-
    const dimensionedScalar& rho1() const;
    
    //-
    const dimensionedScalar& rho2() const;
    
    //-
    const dimensionedScalar& Tau1() const;
    
    //-
    const dimensionedScalar& Tau2() const;

};


}

#endif

//
//END-OF-FILE
//


