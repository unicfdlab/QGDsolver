Info<< "\nReading shallowWaterProperties" << endl;

shallowWaterQGDThermo thermo(mesh);
thermo.correct();

scalar alpha = thermo.alpha();

Switch wellBalancedScheme = thermo.wellBalancedScheme();
Switch dryZoneCondition = thermo.dryZoneCondition();
scalar NS = thermo.NS();
scalar eps0 = thermo.eps0();
scalar tauU = thermo.tauU();

scalar n0 = thermo.n();

scalar kappa = 1.0;

dimensionedScalar n = dimensionedScalar("n",  dimensionSet(0, -1.0/3.0, 1.0, 0, 0, 0, 0), n0);