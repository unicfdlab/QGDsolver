{

    volScalarField  Hi
    (
        IOobject
        (
            "Hi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass
    );
    
    surfaceScalarField Hif
    (
        "Hif",
        Hf*0.0
    );
    
    forAll(Y, iCmpt)
    {
        scalarField& hiIF = Hi.ref();
        const scalarField& pIF  = p.internalField();
        const scalarField& TIF  = T.internalField();
        
        forAll(hiIF, iCell)
        {
            hiIF[iCell] = thermo.composition().Es(iCmpt, pIF[iCell], TIF[iCell]);
        }
        
        forAll(Hi.boundaryField(), iPatch)
        {
            fvPatchScalarField& hip = Hi.boundaryFieldRef()[iPatch];
            const fvPatchScalarField& pp = p.boundaryField()[iPatch];
            const fvPatchScalarField& Tp = T.boundaryField()[iPatch];
            forAll(hip, iFace)
            {
                hip[iFace] = thermo.composition().Hs(iCmpt, pp[iFace], Tp[iFace]);
            }
        }
        
        Hif = qgdInterpolate(Hi);
        
        phiQ +=
            //contribution from relative motion of components (diffusion in Yi transport)
            qgdFlux
            (
                diffusiveFlux[iCmpt],
                Hi,
                Hif
            )
            +
            // contribution due to casting diffusion term from temperature to system enthalpy
            qgdFlux
            (
                fvc::snGrad(Y[iCmpt]) * mesh.magSf(),
                Hi * turbulence->alphaEff(),
                Hif * alphauf
            );
    }
}

//
//END-OF-FILE
//

