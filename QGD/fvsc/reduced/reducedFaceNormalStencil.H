#ifndef reducedFaceNormalStencil_H
#define reducedFaceNormalStencil_H

#include "fvscStencil.H"
#include "volFields.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class reduced Declaration
\*---------------------------------------------------------------------------*/

class reduced
:
    public fvscStencil
{

public:

        TypeName ("reduced");

    // Constructors

        //- Construct from IOobject. 
        reduced(const IOobject&);

        //- Destructor
        ~reduced();
        
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

