#ifndef GaussVolPointBase_H
#define GaussVolPointBase_H

#include "fvMesh.H"

#include "GaussVolPointBase1D.H"
#include "GaussVolPointBase2D.H"
#include "GaussVolPointBase3D.H"

namespace Foam
{

namespace fvsc
{

class GaussVolPointBase :
    public GaussVolPointBase1D,
    public GaussVolPointBase2D,
    public GaussVolPointBase3D
{

    public:

        //-
        GaussVolPointBase(const fvMesh& mesh);

        //-
        ~GaussVolPointBase();

        //-
        void faceGrad (const volScalarField& vF, surfaceVectorField& gradf);

        //-
        void faceGrad (const volVectorField& vF, surfaceTensorField& gradf);

        //-
        void faceDiv (const volVectorField& vVF, surfaceScalarField& divf);

        //-
        void faceDiv(const volTensorField& vTF, surfaceVectorField& divf);
};

}

}

#endif

//
//
//


