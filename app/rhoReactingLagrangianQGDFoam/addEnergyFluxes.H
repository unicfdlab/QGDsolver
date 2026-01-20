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

    volScalarField  Ei
    (
        IOobject
        (
            "Ei",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass
    );
    
    surfaceScalarField Eif
    (
        "Eif",
        Hf*0.0
    );
    
    forAll(Y, iCmpt)
    {
        scalarField& hiIF = Hi.ref();
        scalarField& eiIF = Ei.ref();
        const scalarField& pIF  = p.internalField();
        const scalarField& TIF  = T.internalField();
        
        forAll(hiIF, iCell)
        {
            hiIF[iCell] = thermo.composition().Hs(iCmpt, pIF[iCell], TIF[iCell]);
            eiIF[iCell] = thermo.composition().Es(iCmpt, pIF[iCell], TIF[iCell]);
        }
        
        forAll(Hi.boundaryField(), iPatch)
        {
            fvPatchScalarField& hip = Hi.boundaryFieldRef()[iPatch];
            fvPatchScalarField& eip = Ei.boundaryFieldRef()[iPatch];
            const fvPatchScalarField& pp = p.boundaryField()[iPatch];
            const fvPatchScalarField& Tp = T.boundaryField()[iPatch];
            forAll(hip, iFace)
            {
                hip[iFace] = thermo.composition().Hs(iCmpt, pp[iFace], Tp[iFace]);
                eip[iFace] = thermo.composition().Es(iCmpt, pp[iFace], Tp[iFace]);
            }
        }
        
        Hif = qgdInterpolate(Hi);
        Eif = qgdInterpolate(Ei);
        
        phiQ +=
            //contribution from relative motion of components (diffusion in Yi transport)
            //which acts actually as convection of energy
            qgdFlux
            (
                diffusiveFlux[iCmpt],
                Hi,
                Hif
            );
        phiQ +=
            // contribution due to casting diffusion term from temperature to system enthalpy
            qgdFlux
            (
                fvc::snGrad(Y[iCmpt]) * mesh.magSf(),
                Ei * turbulence->alphaEff(),
                Eif * alphauf
            );
    }
}

//
//END-OF-FILE
//

