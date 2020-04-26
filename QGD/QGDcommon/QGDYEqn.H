{
    combustion->correct();
    Qdot = combustion->Qdot();
    volScalarField Yt(0.0*Y[0]);
    
    forAll(Y, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Yi = Y[i];
            const surfaceScalarField& phiJmYi = phiJmY[i];
            
            if (implicitDiffusion)
            {
                fvScalarMatrix YEqn
                (
                    fvm::ddt(rho,Yi)
                    +
                    fvc::div(phiJmYi)
                    -
                    fvm::laplacian(muf/ScNumbers[i],Yi)
                    ==
                    combustion->R(Yi)
                    +
                    parcels.SYi(i, Yi)
                );
                
                YEqn.solve();
                
                diffusiveFlux[i] += YEqn.flux();
                diffusiveFlux[inertIndex] -= diffusiveFlux[i];
            }
            else
            {
                solve
                (
                    fvm::ddt(rho,Yi)
                    +
                    fvc::div(phiJmYi)
                    -
                    fvc::laplacian(muf/ScNumbers[i],Yi)
                    ==
                    combustion->R(Yi)
                    +
                    parcels.SYi(i, Yi)
                );
                
                diffusiveFlux[i] += (muf / ScNumbers[i]) * fvc::snGrad(Yi.oldTime()) * mesh.magSf();
                diffusiveFlux[inertIndex] -= diffusiveFlux[i];
            }
            
            Yi.max(0.0);
            Yt += Yi;
        }
        
        Y[inertIndex] = scalar(1) - Yt;
        Y[inertIndex].max(0.0);
    }
}

//
//END-OF-FILE
//


