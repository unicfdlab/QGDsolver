/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::psiQGDThermo

Description
    Basic thermodynamic properties based on compressibility

SourceFiles
    psiQGDThermo.C

\*---------------------------------------------------------------------------*/

#ifndef psiQGDThermo_H
#define psiQGDThermo_H

#include "psiThermo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class psiQGDThermo Declaration
\*---------------------------------------------------------------------------*/

class psiQGDThermo
:
    public psiThermo
{

protected:

    // Protected Member Functions

        //- Construct as copy (not implemented)
        psiQGDThermo(const psiQGDThermo&);
        
        //-
        volScalarField muQGD_;
        
        //-
        volScalarField alphahQGD_;
        
        //-
        volScalarField hQGD_;
        
        //-
        volScalarField aQGD_;
        
        //-
        volScalarField gamma_;
        
        //-
        volScalarField c_;
        
        //-
        volScalarField tauQGD_;
        
        //-
        scalar PrQGD_;
        
        //-
        scalar ScQGD_;
        
        //-
        void correctQGD();

        static tmp<volScalarField> qgdLength(const fvMesh& mesh);

public:

    //- Runtime type information
    TypeName("psiQGDThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        psiQGDThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName)
    );


    // Constructors

        //- Construct from mesh and phase name
        psiQGDThermo
        (
            const fvMesh&,
            const word& phaseName
        );


    //- Selector
    static autoPtr<psiQGDThermo> New
    (
        const fvMesh& mesh,
        const word& phaseName=word::null
    );


    //- Destructor
    virtual ~psiQGDThermo();


    // Member functions
    
    //-
    bool read();
    
    //-
    tmp<volScalarField> c() const;
    
    //-
    tmp<volScalarField> hQGD() const;
    
    //-
    dimensionedScalar ScQGD() const;
    
    //-
    tmp<volScalarField> tauQGD() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
