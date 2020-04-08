/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
                QGDsolver   | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
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
    Foam::twoPhaseIcoQGDThermo

Description
    Thermodynamics and mechanics class for incompressible
    two-phase mixture of immiscible components

SourceFiles
    twoPhaseIcoQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseIcoQGDThermo_H
#define twoPhaseIcoQGDThermo_H

#include "IOdictionary.H"
#include "twoPhaseMixture.H"
#include "qInterfaceProperties.H"

#include "QGDThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class twoPhaseIcoQGDThermo Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseIcoQGDThermo
:
    public IOdictionary,
    public twoPhaseMixture,
    public qInterfaceProperties,
    public QGDThermo
{

protected:

    // Protected Member Functions
    
    //- mixture pressure field
    volScalarField p_;
    
    //- density of phase 1
    dimensionedScalar rho1_;
    
    //- density of phase 2
    dimensionedScalar rho2_;
    
    //- kinematic viscosity of phase 1
    dimensionedScalar nu1_;
    
    //- kinematic viscosity of phase 2
    dimensionedScalar nu2_;

    //- regularization time of phase 1
    dimensionedScalar Tau1_;
    
    //- regularization time of phase 2
    dimensionedScalar Tau2_;


    //- Construct as copy (not implemented)
    twoPhaseIcoQGDThermo(const twoPhaseIcoQGDThermo&);


public:

    //- Runtime type information
    TypeName("twoPhaseIcoQGDThermo");

    // Constructors

        //- Construct from mesh
        // and velocity field
        twoPhaseIcoQGDThermo
        (
            const fvMesh&,
            const volVectorField& U
        );

    //- Destructor
    virtual ~twoPhaseIcoQGDThermo();


    // Member functions

    //-
    bool read();
    
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
