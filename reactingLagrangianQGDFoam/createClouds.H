Info<< "\nConstructing reacting cloud" << endl;
basicReactingCloud parcels
(
    "reactingCloud",
    rho,
    U,
    g,
    slgThermo
);