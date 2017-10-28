
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "tmp.H"

#ifndef fvsc_H
#define fvsc_H

namespace Foam
{

namespace fvsc
{

    //returns gradient of volScalarField at the surface mesh
    tmp<surfaceVectorField> grad(const volScalarField& vF);
    
    //-
    tmp<surfaceVectorField> grad(const tmp<volScalarField>& vF);

    //returns gradient of volVectorField at the surface mesh
    tmp<surfaceTensorField> grad(const volVectorField& vF);
    
    //-
    tmp<surfaceTensorField> grad(const tmp<volVectorField>& vF);

    //returns divergence of volVectorField at the surface mesh
    tmp<surfaceScalarField> div(const volVectorField& vF);
    
    //-
    tmp<surfaceScalarField> div(const tmp<volVectorField>& vF);

    //returns divergence of volTensorField at the surface mesh
    tmp<surfaceVectorField> div(const volTensorField& vF);
    
    //-
    tmp<surfaceVectorField> div(const tmp<volTensorField>& vF);

} //namespace fvsc

} //namespace Foam

#endif

//END-OF-FILE

