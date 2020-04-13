
#include "GaussVolPointBase1D.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"

Foam::fvsc::GaussVolPointBase1D::GaussVolPointBase1D(const fvMesh& mesh)
:
    nfRef_(mesh.thisDb().lookupObject<surfaceVectorField>("nf"))
{
};


Foam::fvsc::GaussVolPointBase1D::~GaussVolPointBase1D()
{
}

void Foam::fvsc::GaussVolPointBase1D::faceGrad(const volScalarField& f, surfaceVectorField& gradf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        gradf = nfRef_() * fvc::snGrad(f);
    }
};

void Foam::fvsc::GaussVolPointBase1D::faceGrad(const volVectorField& f, surfaceTensorField& gradf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        gradf = nfRef_() * fvc::snGrad(f);
    }
};

void Foam::fvsc::GaussVolPointBase1D::faceDiv(const volVectorField& f, surfaceScalarField& divf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        divf = nfRef_() & fvc::snGrad(f);
    }
};

void Foam::fvsc::GaussVolPointBase1D::faceDiv(const volTensorField& f, surfaceVectorField& divf)
{
    if (f.mesh().nGeometricD() == 1)
    {
        divf = nfRef_() & fvc::snGrad(f);
    }
}

//
//END-OF-FILE
//


