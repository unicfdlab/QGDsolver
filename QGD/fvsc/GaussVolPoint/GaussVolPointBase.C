
#include "GaussVolPointBase.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "coupledFvPatch.H"
#include "processorFvPatch.H"
#include "emptyFvPatch.H"
#include "wedgeFvPatch.H"
#include "fvcSnGrad.H"

Foam::fvsc::GaussVolPointBase::GaussVolPointBase(const fvMesh& mesh)
    :
    GaussVolPointBase1D(mesh),
    GaussVolPointBase2D(mesh),
    GaussVolPointBase3D(mesh)
{
};

Foam::fvsc::GaussVolPointBase::~GaussVolPointBase()
{
};

void Foam::fvsc::GaussVolPointBase::faceGrad(const volScalarField& vF,
        surfaceVectorField& gradf)
{
    if (vF.mesh().nGeometricD() == 1)
    {
        GaussVolPointBase1D::faceGrad(vF, gradf);
        return;
    }
    else if (vF.mesh().nGeometricD() == 2)
    {
        GaussVolPointBase2D::faceGrad(vF, gradf);
        return;
    }
    else
    {
        GaussVolPointBase3D::faceGrad(vF, gradf);
    }
}

void Foam::fvsc::GaussVolPointBase::faceGrad(const volVectorField& vF,
        surfaceTensorField& gradf)
{
    if (vF.mesh().nGeometricD() == 1)
    {
        GaussVolPointBase1D::faceGrad(vF, gradf);
        return;
    }
    else if (vF.mesh().nGeometricD() == 2)
    {
        surfaceVectorField gradU (vector::zero * fvc::snGrad(vF.component(0)));
        surfaceVectorField gradV (gradU * 0.0);
        surfaceVectorField gradW (gradU * 0.0);
        GaussVolPointBase2D::faceGrad(vF.component(0), gradU);
        GaussVolPointBase2D::faceGrad(vF.component(1), gradV);
        GaussVolPointBase2D::faceGrad(vF.component(2), gradW);
        //Internal field
        gradf.primitiveFieldRef().replace(0, gradU.primitiveField().component(0));
        gradf.primitiveFieldRef().replace(1, gradV.primitiveField().component(0));
        gradf.primitiveFieldRef().replace(2, gradW.primitiveField().component(0));
        gradf.primitiveFieldRef().replace(3, gradU.primitiveField().component(1));
        gradf.primitiveFieldRef().replace(4, gradV.primitiveField().component(1));
        gradf.primitiveFieldRef().replace(5, gradW.primitiveField().component(1));
        gradf.primitiveFieldRef().replace(6, gradU.primitiveField().component(2));
        gradf.primitiveFieldRef().replace(7, gradV.primitiveField().component(2));
        gradf.primitiveFieldRef().replace(8, gradW.primitiveField().component(2));
        //Boundary field
        gradf.boundaryFieldRef().replace(0, gradU.boundaryField().component(0));
        gradf.boundaryFieldRef().replace(1, gradV.boundaryField().component(0));
        gradf.boundaryFieldRef().replace(2, gradW.boundaryField().component(0));
        gradf.boundaryFieldRef().replace(3, gradU.boundaryField().component(1));
        gradf.boundaryFieldRef().replace(4, gradV.boundaryField().component(1));
        gradf.boundaryFieldRef().replace(5, gradW.boundaryField().component(1));
        gradf.boundaryFieldRef().replace(6, gradU.boundaryField().component(2));
        gradf.boundaryFieldRef().replace(7, gradV.boundaryField().component(2));
        gradf.boundaryFieldRef().replace(8, gradW.boundaryField().component(2));
        return;
    }
    else
    {
        GaussVolPointBase3D::faceGrad(vF, gradf);
    }
}

void Foam::fvsc::GaussVolPointBase::faceDiv(const volVectorField& vVF,
        surfaceScalarField& divf)
{
    if (vVF.mesh().nGeometricD() == 1)
    {
        GaussVolPointBase1D::faceDiv(vVF, divf);
        return;
    }
    else if (vVF.mesh().nGeometricD() == 2)
    {
        GaussVolPointBase2D::faceDiv(vVF, divf);
        return;
    }
    else
    {
        GaussVolPointBase3D::faceDiv(vVF, divf);
    }
}

void Foam::fvsc::GaussVolPointBase::faceDiv(const volTensorField& vTF,
        surfaceVectorField& divf)
{
    if (vTF.mesh().nGeometricD() == 1)
    {
        GaussVolPointBase1D::faceDiv(vTF, divf);
        return;
    }
    else if (vTF.mesh().nGeometricD() == 2)
    {
        GaussVolPointBase2D::faceDiv(vTF, divf);
        return;
    }
    else
    {
        GaussVolPointBase3D::faceDiv(vTF, divf);
    }
}

//
//END-OF-FILE
//


