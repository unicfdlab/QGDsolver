#include "fvsc.H"
#include "fvscStencil.H"
#include "volFields.H"
#include "surfaceFields.H"

Foam::tmp<Foam::surfaceVectorField>
Foam::fvsc::grad(const volScalarField& vf)
{
    word opname = "none";
    
    vf.mesh().schemesDict().subDict
    (
        "fvsc"
    ).lookup("grad(" + vf.name() + ")") >> opname;
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        opname,
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
    word opname = "none";
    
    vf.mesh().schemesDict().subDict
    (
        "fvsc"
    ).lookup("grad(" + vf.name() + ")") >> opname;
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        opname,
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
    word opname = "none";
    
    vf.mesh().schemesDict().subDict
    (
        "fvsc"
    ).lookup("div(" + vf.name() + ")") >> opname;
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        opname,
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
    word opname = "none";
    
    vf.mesh().schemesDict().subDict
    (
        "fvsc"
    ).lookup("div(" + vf.name() + ")") >> opname;
    
    fvscStencil& Stencil = fvscStencil::lookupOrNew
    (
        opname,
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


