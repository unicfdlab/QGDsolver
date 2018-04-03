    solve
    (
        fvm::ddt(rhoU)
        +
        fvc::div(phiJmU)
        +
        fvc::div(phiP)
        -
        fvc::div(phiPi)
    );
    
    // Correct velocity
    U.ref() =
        rhoU()
       /rho();
    U.correctBoundaryConditions();
    
   // Solve diffusive QGD & NS part
    if (implicitDiffusion)
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(rho, U) - fvc::ddt(rho,U)
            - fvm::laplacian(muf, U)
            - fvc::div(phiTauMC)
        );
        
        solve
        (
            UEqn
        );
        
        rhoU = rho*U;
        
        sigmaDotUPtr() = (muf*linearInterpolate(fvc::grad(U)) + tauMCPtr()) & Uf;
        
        phiSigmaDotU = mesh.Sf() & sigmaDotUPtr(); //or eqn.flux()?
    }
    rhoU.boundaryFieldRef() == rho.boundaryField()*
    U.boundaryField();
    
//
//END-OF-FILE
//


