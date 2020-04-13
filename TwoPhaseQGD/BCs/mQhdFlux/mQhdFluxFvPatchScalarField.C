/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "mQhdFluxFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mQhdFluxFvPatchScalarField::mQhdFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::mQhdFluxFvPatchScalarField::mQhdFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    patchType() = dict.lookupOrDefault<word>("patchType", word::null);

    if (dict.found("value") && dict.found("gradient"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::mQhdFluxFvPatchScalarField::mQhdFluxFvPatchScalarField
(
    const mQhdFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(p, iF)
{
    patchType() = ptf.patchType();

    // Map gradient. Set unmapped values and overwrite with mapped ptf
    gradient() = 0.0;
    gradient().map(ptf.gradient(), mapper);

    // Evaluate the value field from the gradient if the internal field is valid
    if (notNull(iF))
    {
        if (iF.size())
        {
            // Note: cannot ask for nf() if zero faces

            scalarField::operator=
            (
                //patchInternalField() + gradient()/patch().deltaCoeffs()
                // ***HGW Hack to avoid the construction of mesh.deltaCoeffs
                // which fails for AMI patches for some mapping operations
                patchInternalField()
              + gradient()*(patch().nf() & patch().delta())
            );
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value. This
        // constructor is used when reconstructing fields
        this->map(ptf, mapper);
        }
    }
    else
    {
        // Enforce mapping of values so we have a valid starting value
        this->map(ptf, mapper);
    }
}


Foam::mQhdFluxFvPatchScalarField::mQhdFluxFvPatchScalarField
(
    const mQhdFluxFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf)
{}


Foam::mQhdFluxFvPatchScalarField::mQhdFluxFvPatchScalarField
(
    const mQhdFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//void Foam::mQhdFluxFvPatchScalarField::updateSnGrad
//(
//    const scalarField& snGradp
//)
//{
//     if (updated())
//     {
//         return;
//     }
//
//     curTimeIndex_ = this->db().time().timeIndex();
//
//     gradient() = snGradp;
//     fixedGradientFvPatchScalarField::updateCoeffs();
// }


void Foam::mQhdFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if (this->patch().boundaryMesh().mesh().thisDb().
        foundObject<surfaceScalarField>("phiwm")
       )
    {
        const surfaceScalarField& phiwm =
            this->patch().boundaryMesh().mesh().
            thisDb().lookupObject<surfaceScalarField>
            ("phiwm");

        if (this->patch().boundaryMesh().mesh().thisDb().
            foundObject<surfaceScalarField>("coeffp")
           )
        {
            const surfaceScalarField& coeffp = 
            this->patch().boundaryMesh().mesh().
            thisDb().lookupObject<surfaceScalarField>
            ("coeffp");

            scalarField fluxSnGrad
            (
                phiwm.boundaryField()[patch().index()]
                /
                coeffp.boundaryField()[patch().index()]
                /
                patch().magSf()
            );
            this->gradient() = fluxSnGrad;
        }
    }
    
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::mQhdFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mQhdFluxFvPatchScalarField
    );
}


// ************************************************************************* //
