#include "fvc.H"
#include "reducedFaceNormalStencil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(reducedFaceNormalStencil,0);
    addToRunTimeSelectionTable
    (
        fvscStencil,
        reducedFaceNormalStencil,
        components
    );
}
}

// constructors
Foam::fvsc::reducedFaceNormalStencil::reducedFaceNormalStencil(const IOobject& io)
:
    fvscStencil(io)
{
}

Foam::fvsc::reducedFaceNormalStencil::~reducedFaceNormalStencil()
{
}

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::reducedFaceNormalStencil::Grad(const volScalarField& vF)
{
    surfaceScalarField sF = linearInterpolate(vF);
    
    tmp<surfaceVectorField> tgradIF(fvc::snGrad(vF)  * nf_);
    
    return tgradIF;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::reducedFaceNormalStencil::Grad(const volVectorField& iVF)
{

    tmp<surfaceTensorField> tgradIVF(fvc::snGrad(iVF) * nf_);

    return tgradIVF;
};

Foam::tmp<Foam::surfaceScalarField> Foam::fvsc::reducedFaceNormalStencil::Div(const volVectorField& iVF)
{
    tmp<surfaceScalarField> tdivIVF(fvc::snGrad(iVF) & nf_);
    
    return tdivIVF;
};

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::reducedFaceNormalStencil::Div(const volTensorField& iTF)
{
    tmp<surfaceVectorField> tdivITF(fvc::snGrad(iTF) & nf_);
    
    return tdivITF;
}


//
//END-OF-FILE
//


