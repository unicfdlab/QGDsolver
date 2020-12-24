
#ifdef QHD_ALGORITHM
if (runTime.controlDict().lookupOrDefault<bool>("checkFo",false))
{
    dimensionedScalar smallu ("smallu", U.dimensions(), SMALL);
    scalar maxDiffusionDt =
        min
        (
            0.5*
            magSqr(1.0/mesh.deltaCoeffs())/
            (thermo.tauQGDf()*sqr(mag(phiu/mesh.magSf())+smallu))
        ).value();
    Info << "Foruier limited max Dt: " << maxDiffusionDt << endl;
    
    maxDeltaT1 = min(maxDeltaT1,maxDiffusionDt);
}
#endif

//
//END-OF-FILE
//


