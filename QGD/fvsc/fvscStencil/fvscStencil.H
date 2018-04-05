
#include "fvMesh.H"
#include "regIOobject.H"
#include "runTimeSelectionTables.H"
#include "surfaceFields.H"

#ifndef fvscStencil_H
#define fvscStencil_H

namespace Foam
{

namespace fvsc
{

class fvscStencil : public regIOobject, public refCount
{

protected:
    
    //-
    const fvMesh& mesh_;
    
    //-
    const Time& runTime_;
    
    //-
    surfaceVectorField nf_;
    
    //-
    static PtrList<fvscStencil> stencils_;

public:
 
    //-
    TypeName("fvscStencil");
    
    //-
    declareRunTimeSelectionTable
    (
        autoPtr,
        fvscStencil,
        components,
        (
            const IOobject& io
        ),
        (io)
    );
    
    //- Construct from components
    fvscStencil
    (
        const IOobject& io
    );
    
    //- Return a reference to the selected fvscStencil model
    static autoPtr<fvscStencil> New
    (
        const word& name,
        const fvMesh& mesh
    );
    
    //-
    //static tmp<fvscStencil> lookupOrNew
    static fvscStencil& lookupOrNew
    (
        const word& nname,
        const fvMesh& mesh
    );
    
    //-
    virtual ~fvscStencil();
    
    //-
    virtual tmp<surfaceVectorField> Grad(const volScalarField& vF)
    {
        notImplemented("tmp<surfaceVectorField> Grad(const volScalarField& vF)");
        return tmp<surfaceVectorField>(NULL);
    }

    //-
    virtual tmp<surfaceTensorField> Grad(const volVectorField& iVF)
    {
        notImplemented("tmp<surfaceTensorField> Grad(const volVectorField& vF)");
        return tmp<surfaceTensorField>(NULL);
    }

    //-
    virtual tmp<surfaceScalarField> Div(const volVectorField& iVF)
    {
        notImplemented("tmp<surfaceScalarField> Grad(const volVectorField& vF)");
        return tmp<surfaceScalarField>(NULL);
    }

    //-
    virtual tmp<surfaceVectorField> Div(const volTensorField& iTF)
    {
        notImplemented("tmp<surfaceVectorField> Grad(const volTensorField& vF)");
        return tmp<surfaceVectorField>(NULL);
    }


    virtual bool writeData(Ostream&) const
    {
        return true;
    }
};

}

}

#endif

//END-OF-FILE

