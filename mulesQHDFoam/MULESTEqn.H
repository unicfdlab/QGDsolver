{
    word alphaScheme("div(phiTf)");
    word alpharScheme("laplacian(Hif,T)");

    tmp<fv::ddtScheme<scalar>> ddtAlpha
    (
        fv::ddtScheme<scalar>::New
        (
            mesh,
            mesh.ddtScheme("ddt(T)")
        )
    );

    MULES::explicitSolve(T, phi, phiTf, gMax(T), gMin(T));

    solve(fvm::ddt(T)-fvc::ddt(T)-fvc::laplacian(Hif,T));

    Info <<"exit mules" << endl;
}

