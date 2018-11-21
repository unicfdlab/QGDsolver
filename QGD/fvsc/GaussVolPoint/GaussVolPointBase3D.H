
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "tmp.H"
#include "List.H"
#include "DynamicList.H"
#include "label.H"
#include "scalar.H"
#include "autoPtr.H"

#ifndef GaussVolPointBase3D_H
#define GaussVolPointBase3D_H

namespace Foam
{

class fvMesh;
class face;
class volPointInterpolation;

typedef List<face> faceList;

namespace fvsc
{

class GaussVolPointBase3D
{
    
    //-
    const volPointInterpolation& volPoint_;

    //+
    tmp<surfaceVectorField> nfRef_;

    //+
    List<List<label> > bgfid_;
    
    //+
    List<bool> processorPatch_;

    /* coefficients for quad faces */

    //+
    DynamicList<label> qf_;
    
    //+
    List<List<scalar> > aqx_;
    
    //+
    List<List<scalar> > aqy_;
    
    //+
    List<List<scalar> > aqz_;
    
    //+
    List<scalar> vq_;
    
    //+
    List<List<label> > bqf_;
    
    //+
    List<List<List<scalar> > > baqx_;

    //+
    List<List<List<scalar> > > baqy_;

    //+
    List<List<List<scalar> > > baqz_;
    
    //+
    List<List<scalar> > bvq_;
    
    //+
    List<List<scalar> > bmv65_;

    /* coefficients for other faces */
    
    //+
    DynamicList<label> of_;
    
    //+
    List<List<label> > bof_;

protected:

    //-
    void calcGradfIF
    (
        const volScalarField& sf,
        const pointScalarField& pf,
        const faceList& faces,
        surfaceVectorField& gradf,
        const surfaceVectorField& dfdn
    );
    
    //-
    void calcGradfBF
    (
        const volScalarField& sf,
        const pointScalarField& pf,
        const faceList& faces,
        surfaceVectorField& gradf,
        const surfaceVectorField& dfdn
    );

    //-
    void calcGradfIF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceTensorField& gradf,
        const surfaceTensorField& dfdn
    );
    
    //-
    void calcGradfBF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceTensorField& gradf,
        const surfaceTensorField& dfdn
    );

    //-
    void calcDivfIF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceScalarField& divf,
        const surfaceScalarField& dfdn
    );

    //-
    void calcDivfBF
    (
        const volVectorField& sf,
        const pointVectorField& pf,
        const faceList& faces,
        surfaceScalarField& divf,
        const surfaceScalarField& dfdn
    );

    //-
    void calcDivfIF
    (
        const volTensorField& tf,
        const pointTensorField& pf,
        const faceList& faces,
        surfaceVectorField& divf,
        const surfaceVectorField& dfdn
    );

    //-
    void calcDivfBF
    (
        const volTensorField& sf,
        const pointTensorField& pf,
        const faceList& faces,
        surfaceVectorField& divf,
        const surfaceVectorField& dfdn
    );

public:

    //-
    GaussVolPointBase3D(const fvMesh& mesh);
    
    //-
    ~GaussVolPointBase3D();
    
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

