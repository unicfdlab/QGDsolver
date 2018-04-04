#ifndef leastSquaresStencil_H
#define leastSquaresStencil_H

#include "fvscStencil.H"
#include "regIOobject.H"
#include "labelList.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"
#include "vector.H"
#include "List.H"

#include "processorFvPatch.H"
#include "leastSquaresBase.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class leastSquares Declaration
\*---------------------------------------------------------------------------*/

class leastSquares
:
    public fvscStencil, public leastSquaresBase
{

public:

        TypeName ("leastSquares");

// Constructors

        //- Construct from IOobject. Optional flag for if IOobject is the
        //  top level regIOobject.
        leastSquares(const IOobject&);
        
        //- Destructor
        ~leastSquares();

        tmp<surfaceVectorField> Grad(const volScalarField& iF);

        tmp<surfaceTensorField> Grad(const volVectorField& iVF);

        tmp<surfaceScalarField> Div(const volVectorField& iVF);

        tmp<surfaceVectorField> Div(const volTensorField& iTF);
};

} //fvsc

} //Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

