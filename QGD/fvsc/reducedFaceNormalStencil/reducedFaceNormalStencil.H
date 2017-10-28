#ifndef reducedFaceNormalStencil_H
#define reducedFaceNormalStencil_H

#include "fvscStencil.H"
#include "volFields.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class extendedFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class reducedFaceNormalStencil
:
    public fvscStencil
{

public:

        TypeName ("reducedFaceNormalStencil");

    // Constructors

        //- Construct from IOobject. 
        reducedFaceNormalStencil(const IOobject&);

        //- Destructor
        ~reducedFaceNormalStencil();
        
        //-
        tmp<surfaceVectorField> Grad(const volScalarField& iF);
        
        //-
        tmp<surfaceTensorField> Grad(const volVectorField& iVF);
        
        //-
        tmp<surfaceScalarField> Div(const volVectorField& iVF);
        
        //-
        tmp<surfaceVectorField> Div(const volTensorField& iTF);
};

} //namespace fvsc

} //namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

//END-OF-FILE

