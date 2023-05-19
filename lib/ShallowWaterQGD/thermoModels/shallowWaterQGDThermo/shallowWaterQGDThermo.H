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

Class
    Foam::shallowWaterQGDThermo

Group
    grpShallowWaterQGDThermo

Description
    Thermodynamics and mechanics class for incompressible
    two-phase mixture of immiscible components

SourceFiles
    shallowWaterQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef shallowWaterQGDThermo_H
#define shallowWaterQGDThermo_H

#include "IOdictionary.H"
#include "QGDThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class shallowWaterQGDThermo Declaration
\*---------------------------------------------------------------------------*/

class shallowWaterQGDThermo
:
    public IOdictionary,
    public QGDThermo
{

protected:

    // Protected Member Functions
    
    //- Tuning parameter alpha
    scalar alpha_;   
    
    //- using well-balanced scheme flag
    Switch wellBalancedScheme_;
    
    //- using dry zone condition flag
    Switch dryZoneCondition_; 
    
    //- Additional regularizer of the type of Navier–Stokes viscosity coefficient
    scalar NS_;
    
    //- cut-off parameter for dry zone condition
    scalar eps0_;
    
    //- Additional term for calculating tau
    scalar tauU_;
    
    //- Manning coefficient
    scalar n_; 
     
    //- Construct as copy (not implemented)
    shallowWaterQGDThermo(const shallowWaterQGDThermo&);
    
    //-
    void validateQGDCoeffs();


public:

    //- Runtime type information
    TypeName("shallowWaterQGDThermo");

    // Constructors

        //- Construct from mesh
        // and velocity field
        shallowWaterQGDThermo
        (
            const fvMesh&
        );

    //- Destructor
    virtual ~shallowWaterQGDThermo();


    // Member functions

    //-
    bool read();
    
    //-
    const scalar& alpha() const;
    
    //-
    const scalar& NS() const;
    
    //-
    const scalar& eps0() const;
    
    //-
    const scalar& tauU() const;
    
    //-
    const scalar& n() const;
    
    //-
    const Switch& wellBalancedScheme() const;
    
    //-
    const Switch& dryZoneCondition() const;
    
    //-
    virtual const volScalarField& p() const;

    //-
    volScalarField& p();

    //-
    virtual tmp<volScalarField> rho() const;

    //-
    virtual tmp<volScalarField> mu() const;
    
    //- Deprecated for incompressible models
    virtual const volScalarField& c() const;
    
    //-
    virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
