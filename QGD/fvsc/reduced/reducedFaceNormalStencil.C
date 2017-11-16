#include "fvc.H"
#include "reducedFaceNormalStencil.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace fvsc
{
    defineTypeNameAndDebug(reduced,0);
    addToRunTimeSelectionTable
    (
        fvscStencil,
        reduced,
        components
    );
}
}

// constructors
Foam::fvsc::reduced::reduced(const IOobject& io)
:
    fvscStencil(io)
{
}

Foam::fvsc::reduced::~reduced()
{
}

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::reduced::Grad(const volScalarField& vF)
{
    surfaceScalarField sF = linearInterpolate(vF);
    
    tmp<surfaceVectorField> tgradIF(nf_ * fvc::snGrad(vF));
    
    return tgradIF;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
Foam::tmp<Foam::surfaceTensorField> Foam::fvsc::reduced::Grad(const volVectorField& iVF)
{

    tmp<surfaceTensorField> tgradIVF(nf_* fvc::snGrad(iVF));

    return tgradIVF;
};

Foam::tmp<Foam::surfaceScalarField> Foam::fvsc::reduced::Div(const volVectorField& iVF)
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
Foam::tmp<Foam::surfaceVectorField> Foam::fvsc::reduced::Div(const volTensorField& iTF)
{
    tmp<surfaceVectorField> tdivITF(Foam::T(fvc::snGrad(iTF)) & nf_);
    
    return tdivITF;
}


//
//END-OF-FILE
//


