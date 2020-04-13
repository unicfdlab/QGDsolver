
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "tmp.H"

#ifndef GaussVolPointBase1D_H
#define GaussVolPointBase1D_H

namespace Foam
{

class fvMesh;

namespace fvsc
{

class GaussVolPointBase1D
{

    //-
    tmp<surfaceVectorField> nfRef_;

public:

    //-
    GaussVolPointBase1D(const fvMesh& mesh);
    
    //-
    ~GaussVolPointBase1D();
    
    //-
    void faceGrad(const volScalarField& f, surfaceVectorField& gradf);

    //-
    void faceGrad(const volVectorField& f, surfaceTensorField& gradf);

    //-
    void faceDiv(const volVectorField& f, surfaceScalarField& divf);
    
    //-
    void faceDiv(const volTensorField& f, surfaceVectorField& divf);

};

}

}

#endif

//
//END-OF-FILE
//

