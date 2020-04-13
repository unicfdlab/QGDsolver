
namespace Foam
{

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdInterpolate(const GeometricField<T, Foam::fvPatchField, Foam::volMesh > &psi)
{
    if
    (
        psi.mesh().schemesDict().subDict("interpolationSchemes").found
        ("interpolate("+psi.name()+")")
    )
    {
        return fvc::interpolate(psi);
    }
    else if 
    (
        psi.mesh().schemesDict().subDict("interpolationSchemes").found
        ("default")
    )
    {
        if
        (
            psi.mesh().schemesDict().subDict("interpolationSchemes").template get<word>
            ("default") == "none"
        )
        {
            return linearInterpolate(psi);
        }
        return fvc::interpolate(psi);
    }
    return linearInterpolate(psi);
}

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdInterpolate(const tmp<GeometricField<T, Foam::fvPatchField, Foam::volMesh >> &tpsi)
{
    return qgdInterpolate(tpsi());
}

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdFlux
(
    const GeometricField<scalar, Foam::fvsPatchField, Foam::surfaceMesh>& flux,
    const GeometricField<T, Foam::fvPatchField, Foam::volMesh > &psi,
    const GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh>& psif,
    const word fluxName
)
{
    if (psi.mesh().schemesDict().subDict("divSchemes").found(fluxName))
    {
        //if
        //(
        //    psi.mesh().schemesDict().subDict("divSchemes").template get<word>(fluxName)
        //    ==
        //    "linear"
        //)
        //{
        //    return flux*psif;
        //}
        return fvc::flux
        (
            flux,
            psi,
            fluxName
        );
    }
    return flux*psif;
}

template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdFlux
(
    const GeometricField<scalar, Foam::fvsPatchField, Foam::surfaceMesh>& flux,
    const GeometricField<T, Foam::fvPatchField, Foam::volMesh > &psi,
    const GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh>& psif
)
{
    word fluxName = "div("+flux.name()+","+psi.name()+")";
    return qgdFlux(flux,psi,psif,fluxName);
}


template<class T>
tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh> >
qgdFlux
(
    const tmp<GeometricField<scalar, Foam::fvsPatchField, Foam::surfaceMesh>>& flux,
    const tmp<GeometricField<T, Foam::fvPatchField, Foam::volMesh >> &psi,
    const tmp<GeometricField<T, Foam::fvsPatchField, Foam::surfaceMesh>>& psif
)
{
    return qgdFlux(flux(),psi(),psif());
}

}

//
//END-OF-FILE
//

