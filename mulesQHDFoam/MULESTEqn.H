{
    Field<scalar> maxT (T.size(), gMax(T));
    Field<scalar> minT (T.size(), gMin(T));
    
    MULES::explicitSolve
    (
        geometricOneField(),
        T,
        phi,
        phiTf,
        zeroField(),
        zeroField(),
        maxT,
        minT
    );

    if (implicitDiffusion)
    {
        solve(fvm::ddt(T) - fvc::ddt(T) - fvm::laplacian(Hif, T) == TSu);
    }
    else
    {
        solve(fvm::ddt(T) - fvc::ddt(T) - fvc::laplacian(Hif, T) == TSu);
    }
}

