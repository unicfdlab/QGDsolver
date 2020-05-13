/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of QGDsolver, based on OpenFOAM library.
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
    
Group
    grpTwoPhaseIcoQGDThermo
    
Description
    Constant properties for two phase package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).
SourceFiles
    constTwoPhaseProperties.C

\*---------------------------------------------------------------------------*/

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
