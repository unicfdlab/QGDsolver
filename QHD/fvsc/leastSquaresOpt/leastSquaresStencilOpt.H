#ifndef leastSquaresOptStencil_H
#define leastSquaresOptStencil_H

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
#include "leastSquaresBase.H"

namespace Foam
{

namespace fvsc
{

/*---------------------------------------------------------------------------*\
                    Class leastSquaresOpt Declaration
\*---------------------------------------------------------------------------*/

class leastSquaresOpt
:
    public fvscStencil, public leastSquaresBase
{

protected:

        void faceScalarDer(const Field<scalar>& iF,const Field<scalar>& sF,int com, surfaceScalarField& rField);
        void faceScalarDer(const tmp<Field<scalar>>& tiF,const tmp<Field<scalar>>& tsF, int com, tmp<surfaceScalarField>& trField);
        
        void faceScalarDer
        (
            const List3<scalar>& procVfValues,
            const surfaceScalarField& sF,
            int derComp,
            surfaceScalarField& rField
        );
        
        void faceScalarDer
        (
            const tmp<List3<scalar>>& tprocVfValues,
            const tmp<surfaceScalarField>& tsF,
            int derComp,
            tmp<surfaceScalarField>& trField
        );

        #include "extendedFaceStencilFormVfValues.H"

public:

        TypeName ("leastSquaresOpt");

// Constructors

        //- Construct from IOobject. Optional flag for if IOobject is the
        //  top level regIOobject.
        leastSquaresOpt(const IOobject&);
        
        //- Destructor
        ~leastSquaresOpt();

        tmp<surfaceVectorField> Grad(const volScalarField& iF);
        tmp<surfaceVectorField> Grad(const volScalarField& iF, const surfaceScalarField&);

        tmp<surfaceTensorField> Grad(const volVectorField& iVF);
        tmp<surfaceTensorField> Grad(const volVectorField& iVF, const surfaceVectorField&);

        tmp<surfaceScalarField> Div(const volVectorField& iVF);

        tmp<surfaceVectorField> Div(const volTensorField& iTF);
        
};

} //fvsc

} //Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

