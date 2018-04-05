#include "fvsc.H"
#include "fvscStencil.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
namespace fvsc
{
    word fvscOpName(const Foam::fvMesh& mesh, Foam::word termName);
}
}

Foam::word
Foam::fvsc::fvscOpName(const Foam::fvMesh& mesh, Foam::word termName)
{
    word opname = "none";
    if (mesh.schemesDict().subDict("fvsc").found(termName))
    {
        mesh.schemesDict().subDict("fvsc").lookup(termName) >> opname;
    }
    else
    {
        mesh.schemesDict().subDict("fvsc").lookup("default") >> opname;
    }
    
    return opname;
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::grad(const volScalarField& vf)
{
    word tname = "grad(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    tmp<surfaceVectorField> tGrad(Stencil.Grad(vf));
    
    return tGrad;
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::grad(const tmp<volScalarField>& tvf)
{
    return Foam::fvsc::grad(tvf());
}

Foam::tmp<Foam::surfaceTensorField>
Foam::fvsc::grad(const volVectorField& vf)
{
    word tname = "grad(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    return Stencil.Grad(vf);
}

Foam::tmp<Foam::surfaceTensorField>
Foam::fvsc::grad(const tmp<volVectorField>& tvf)
{
    return Foam::fvsc::grad(tvf());
}

Foam::tmp<Foam::surfaceScalarField>
Foam::fvsc::div(const volVectorField& vf)
{
    word tname = "div(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    return Stencil.Div(vf);
}

Foam::tmp<Foam::surfaceScalarField>
Foam::fvsc::div(const tmp<volVectorField>& tvf)
{
    return Foam::fvsc::div(tvf());
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::div(const volTensorField& vf)
{
    word tname = "div(" + vf.name() + ")";
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        fvscOpName(vf.mesh(),tname),
        vf.mesh()
    );
    
    return Stencil.Div(vf);
}

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::div(const tmp<volTensorField>& tvf)
{
    return Foam::fvsc::div(tvf());
}


//END-OF-FILE


