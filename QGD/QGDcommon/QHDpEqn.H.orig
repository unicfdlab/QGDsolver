{
    //Continuity equation
    p.correctBoundaryConditions();
    fvScalarMatrix pEqn
    (
         fvc::div(phiu)
        -fvc::div(phiwo)
        -fvm::laplacian(taubyrhof,p)
    );

    pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

    pEqn.solve();

    phi = phiu - phiwo + pEqn.flux();
}
