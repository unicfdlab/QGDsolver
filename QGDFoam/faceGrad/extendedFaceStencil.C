#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>

// constructors
Foam::extendedFaceStencil::extendedFaceStencil(const IOobject& io, const bool isTime)
:
    regIOobject(io, isTime),
    mesh_(refCast<const fvMesh>(io.db()))
{
    findNeighbours();
    calculateWeights();
}


Foam::extendedFaceStencil::extendedFaceStencil(const regIOobject& rio)
:
    regIOobject(rio),
    mesh_(refCast<const fvMesh>(rio.db()))
{
    findNeighbours();
    calculateWeights();
}

Foam::extendedFaceStencil::extendedFaceStencil(const regIOobject& rio, bool registerCopy)
:
    regIOobject(rio, registerCopy),
    mesh_(refCast<const fvMesh>(rio.db()))
{
    findNeighbours();
    calculateWeights();
}

Foam::extendedFaceStencil::~extendedFaceStencil()
{
}

void Foam::extendedFaceStencil::rename(const word& newName)
{
};

bool Foam::extendedFaceStencil::readData(Istream&)
{
  return false;
};

bool Foam::extendedFaceStencil::read()
{
  return false;
};

bool Foam::extendedFaceStencil::modified()  const
{
  return false;
};

bool Foam::extendedFaceStencil::readIfModified()
{
  return false;
};

bool Foam::extendedFaceStencil::writeData(Ostream&) const
{
  return false;
};

bool Foam::extendedFaceStencil::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
)  const
{
  return false;
};

bool Foam::extendedFaceStencil::write()  const
{
  return false;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
tmp<surfaceTensorField> Foam::extendedFaceStencil::faceVectorGrad(const volVectorField& iVF)
{
    surfaceVectorField gradComp0col = faceScalarGrad(iVF.component(0));
    surfaceVectorField gradComp1col = faceScalarGrad(iVF.component(1));
    surfaceVectorField gradComp2col = faceScalarGrad(iVF.component(2));

    tmp<surfaceTensorField> tgradIVF(0*fvc::snGrad(iVF) * mesh_.Sf() / mesh_.magSf());
    surfaceTensorField& gradIVF = tgradIVF.ref();

    forAll(mesh_.faces(), facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            gradIVF[facei].component(0) = gradComp0col[facei].component(0);
            gradIVF[facei].component(1) = gradComp1col[facei].component(0);
            gradIVF[facei].component(2) = gradComp2col[facei].component(0);

            gradIVF[facei].component(3) = gradComp0col[facei].component(1);
            gradIVF[facei].component(4) = gradComp1col[facei].component(1);
            gradIVF[facei].component(5) = gradComp2col[facei].component(1);

            gradIVF[facei].component(6) = gradComp0col[facei].component(2);
            gradIVF[facei].component(7) = gradComp1col[facei].component(2);
            gradIVF[facei].component(8) = gradComp2col[facei].component(2);
        }
    }
    
    //processor & coupled boundaries
    forAll(mesh_.boundary(), patchIndex)
    {
        if (isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
        {
            forAll(mesh_.boundary()[patchIndex], facei)
            {
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(0) = 
                    gradComp0col.boundaryField()[patchIndex][facei].component(0);
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(1) = 
                    gradComp1col.boundaryField()[patchIndex][facei].component(0);
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(2) = 
                    gradComp2col.boundaryField()[patchIndex][facei].component(0);

                gradIVF.boundaryFieldRef()[patchIndex][facei].component(3) = 
                    gradComp0col.boundaryField()[patchIndex][facei].component(1);
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(4) = 
                    gradComp1col.boundaryField()[patchIndex][facei].component(1);
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(5) = 
                    gradComp2col.boundaryField()[patchIndex][facei].component(1);

                gradIVF.boundaryFieldRef()[patchIndex][facei].component(6) = 
                    gradComp0col.boundaryField()[patchIndex][facei].component(2);
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(7) = 
                    gradComp1col.boundaryField()[patchIndex][facei].component(2);
                gradIVF.boundaryFieldRef()[patchIndex][facei].component(8) = 
                    gradComp2col.boundaryField()[patchIndex][facei].component(2);
            }
        }
    }

    return tgradIVF;
};

tmp<surfaceTensorField> Foam::extendedFaceStencil::faceVectorGrad(const tmp<volVectorField>& tiVF)
{
    return faceVectorGrad(tiVF());
}

//- Calculate divergence of volume vector field on the faces.
//
// \param iVF        Internal vector field.
//                   Allowable values: constant reference to the volVectorField.
//
// \return           Divergence of iVF (scalar field) which was computed on the faces of mesh.
tmp<surfaceScalarField> Foam::extendedFaceStencil::faceVectorDiv(const volVectorField& iVF)
{
    surfaceVectorField gradComp0 = faceScalarGrad(iVF.component(0));
    surfaceVectorField gradComp1 = faceScalarGrad(iVF.component(1));
    surfaceVectorField gradComp2 = faceScalarGrad(iVF.component(2));

    tmp<surfaceScalarField> tdivIVF(0*fvc::snGrad(iVF) & mesh_.Sf() / mesh_.magSf());
    surfaceScalarField& divIVF = tdivIVF.ref();

    forAll(mesh_.faces(), facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            divIVF[facei] = gradComp0[facei].component(0) + gradComp1[facei].component(1) + gradComp2[facei].component(2);
        }
    }

    forAll(mesh_.boundary(), patchIndex)
    {
        if (isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
        {
            divIVF.boundaryFieldRef()[patchIndex] = 
                gradComp0.boundaryField()[patchIndex].component(0)
                +
                gradComp1.boundaryField()[patchIndex].component(1)
                +
                gradComp2.boundaryField()[patchIndex].component(2);
        }
    }
    
    return tdivIVF;
};

tmp<surfaceScalarField> Foam::extendedFaceStencil::faceVectorDiv(const tmp<volVectorField>& tiVF)
{
    return faceVectorDiv(tiVF());
}

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
tmp<surfaceVectorField> Foam::extendedFaceStencil::faceTensorDiv(const volTensorField& iTF)
{
    tmp<surfaceVectorField> gradComp0 (faceScalarGrad(iTF.component(0)));
    tmp<surfaceVectorField> gradComp1 (faceScalarGrad(iTF.component(1)));
    tmp<surfaceVectorField> gradComp2 (faceScalarGrad(iTF.component(2)));
    
    tmp<surfaceVectorField> gradComp3 (faceScalarGrad(iTF.component(3)));
    tmp<surfaceVectorField> gradComp4 (faceScalarGrad(iTF.component(4)));
    tmp<surfaceVectorField> gradComp5 (faceScalarGrad(iTF.component(5)));

    tmp<surfaceVectorField> gradComp6 (faceScalarGrad(iTF.component(6)));
    tmp<surfaceVectorField> gradComp7 (faceScalarGrad(iTF.component(7)));
    tmp<surfaceVectorField> gradComp8 (faceScalarGrad(iTF.component(8)));

    tmp<surfaceScalarField> divComp0 (gradComp0().component(0) + gradComp3().component(1) + gradComp6().component(2));
    tmp<surfaceScalarField> divComp1 (gradComp1().component(0) + gradComp4().component(1) + gradComp7().component(2));
    tmp<surfaceScalarField> divComp2 (gradComp2().component(0) + gradComp5().component(1) + gradComp8().component(2));

    tmp<surfaceVectorField> tdivITF(0*fvc::snGrad(iTF.component(0)) * mesh_.Sf() / mesh_.magSf());
    surfaceVectorField& divITF = tdivITF.ref();
    
    forAll(mesh_.faces(), facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            divITF[facei].component(0) = divComp0()[facei];
            divITF[facei].component(1) = divComp1()[facei];
            divITF[facei].component(2) = divComp2()[facei];
        }
    }
    
    forAll(mesh_.boundary(), patchIndex)
    {
        if (isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
        {
            forAll(mesh_.boundary()[patchIndex], facei)
            {
                divITF.boundaryFieldRef()[patchIndex][facei].component(0) = 
                    divComp0().boundaryField()[patchIndex][facei];
                divITF.boundaryFieldRef()[patchIndex][facei].component(1) = 
                    divComp1().boundaryField()[patchIndex][facei];
                divITF.boundaryFieldRef()[patchIndex][facei].component(2) = 
                    divComp2().boundaryField()[patchIndex][facei];
            }
        }
    }
    
    return tdivITF;
}

tmp<surfaceVectorField> Foam::extendedFaceStencil::faceTensorDiv(const tmp<volTensorField>& tiTF)
{
    return faceTensorDiv(tiTF());
}

//
//END-OF-FILE
//


