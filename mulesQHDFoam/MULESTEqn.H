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

    if (MULESCorr)
    {
        fvScalarMatrix TEqn
        (
              fvm::ddt(T)
            + fvc::div(phiTf)
            - fvc::laplacian(Hif,T)
               
        );

        TEqn.solve();
        
        alpha1 = T;


        if (alphaApplyPrevCorr && talphaPhiCorr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct(alpha1, phiTf, talphaPhiCorr0.ref(), 1, 0);

            phiTf += talphaPhiCorr0();
        }

        // Cache the upwind-flux
        talphaPhiCorr0 = phiTf;

    }
    else
    {
        MULES::explicitSolve(alpha1, phi, phiTf, 1, 0);
    }
}

