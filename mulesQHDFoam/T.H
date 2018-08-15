for
(
    subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
    !(++alphaSubCycle).end();
)
{
    tmp<surfaceScalarField> tTPhiCorr0;
    
    word TFScheme("div(phiTf)")
    word TLScheme("laplacian(Hif,T)")

    tmp<fv::ddtScheme<scalar>> ddtT
    (
        fv::ddtScheme<scalar>::New
        (
            mesh,
            mesh.ddtScheme("ddt(T)")
        )
    );

    // Set the off-centering coefficient according to ddt scheme
    scalar ocCoeff = 0;
    if
    (
        isType<fv::EulerDdtScheme<scalar>>(ddtT())
     || isType<fv::localEulerDdtScheme<scalar>>(ddtT())
    )
    {
        ocCoeff = 0;
    }
    else if (isType<fv::CrankNicolsonDdtScheme<scalar>>(ddtT()))
    {
        if (nAlphaSubCycles > 1)
        {
            FatalErrorInFunction
                << "Sub-cycling is not supported "
                   "with the CrankNicolson ddt scheme"
                << exit(FatalError);
        }

        ocCoeff =
            refCast<const fv::CrankNicolsonDdtScheme<scalar>>(ddtT())
           .ocCoeff();
    }
    else
    {
        FatalErrorInFunction
            << "Only Euler and CrankNicolson ddt schemes are supported"
            << exit(FatalError);
    }

    scalar cnCoeff = 1.0/(1.0 + ocCoeff);

    // Standard face-flux compression coefficient
    surfaceScalarField phic(mixture.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (mixture.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    surfaceScalarField::Boundary& phicBf =
        phic.boundaryFieldRef();


    // Do not compress interface at non-coupled boundary faces
    // (inlets, outlets etc.)
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phicBf[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }

    tmp<surfaceScalarField> phiCN(phi);

   // Calculate the Crank-Nicolson off-centred volumetric flux
    if (ocCoeff > 0)
    {
        phiCN = cnCoeff*phi + (1.0 - cnCoeff)*phi.oldTime();
    }

    if (MULESCorr)
    {
        fvScalarMatrix TEqn
        (
              fvm::ddt(T)
            + fvc::div(phiTf)
            - fvc::laplacian(Hif,T)
        );

	TEqn.solve();

	tmp<surfaceScalarField> tPhiUD(TEqn.flux());
	tPhi = tPhiUD();

        if (alphaApplyPrevCorr && tTPhiCorr0.valid())
        {
            Info<< "Applying the previous iteration compression flux" << endl;
            MULES::correct(T, phiTf, tTPhiCorr0.ref(), 1, 0);

            PhiTf += tTPhiCorr0();
        }

        // Cache the upwind-flux
        tTPhiCorr0 = tTPhiUD;

	mixture.correct();
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(phic*mixture.nHatf());

        tmp<surfaceScalarField> tTPhiUn
        (
            fvc::flux
            (
                phi,
                T,
                TScheme
            )
          + fvc::flux
            (
               -fvc::flux(-phir, T, TLScheme),
                T,
                TLScheme
            )
        );
     
     


    if (ocCoeff > 0)
        {
            tTPhiUn =
                cnCoeff*tTPhiUn + (1.0 - cnCoeff)*phiTf.oldTime();
        }

    if (MULESCorr)
        {
            tmp<surfaceScalarField> tTPhiCorr(tTPhiUn() - phiTf);
            volScalarField T10("T10", T);

            MULES::correct(T, tTPhiUn(), tTPhiCorr.ref(), 1, 0);

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                phiTf += tTPhiCorr();
            }
            else
            {
                T = 0.5*T + 0.5*T10;
                phiTf += 0.5*tTPhiCorr();
            }
        }
        else
        {
            phiTf = tTPhiUn;

            MULES::explicitSolve(T1, phiCN, phiTf, 1, 0);
        }

	mixture.correct();
    }

}
