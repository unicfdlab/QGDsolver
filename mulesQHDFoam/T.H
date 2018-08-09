{
    if (MULESCorr)
    {
        solve
        (
              fvm::ddt(T)
            + fvc::div(phiTf)
            - fvc::laplacian(Hif,T)
        );
    }
}
