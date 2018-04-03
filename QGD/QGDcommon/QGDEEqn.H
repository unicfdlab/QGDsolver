    solve
    (
        fvm::ddt(rhoE)
      + fvc::div(phiJmH)
      + fvc::div(phiQ)
      - fvc::div(phiPiU)
      - fvc::div(phiSigmaDotU)
    );

    // Correct energy
    e = rhoE/rho - 0.5*magSqr(U);
    e.correctBoundaryConditions();
    
    // Solve diffusive QGD & NS part
    if (implicitDiffusion)
    {
        solve
        (
            fvm::ddt(rho, e) - fvc::ddt(rho,e)
          - fvm::laplacian(alphauf, e)
        );
        
        rhoE = rho*(e + 0.5*magSqr(U));
    }
    
    rhoE.boundaryFieldRef() == rho.boundaryField()*
        (e.boundaryField() + 0.5*magSqr(U.boundaryField()));
//
//END-OF-FILE
//


