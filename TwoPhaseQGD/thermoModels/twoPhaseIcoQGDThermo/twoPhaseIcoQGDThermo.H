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
    Foam::twoPhaseIcoQGDThermo

Group
    grpTwoPhaseIcoQGDThermo

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
#include "constTwoPhaseProperties.H"
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
    public constTwoPhaseProperties,
    public qInterfaceProperties,
    public QGDThermo
{

protected:

    // Protected Member Functions
    
    //- mixture pressure field
    volScalarField p_;
    
    //- Construct as copy (not implemented)
    twoPhaseIcoQGDThermo(const twoPhaseIcoQGDThermo&);
    
    //-
    void validateQGDCoeffs();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
