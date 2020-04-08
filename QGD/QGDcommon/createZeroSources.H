
fvScalarMatrix rhoSu
(
    rho,
    rho.dimensions()*dimVolume/dimTime
);

fvVectorMatrix rhoUSu
(
    U,
    rho.dimensions()*U.dimensions()*dimVolume/dimTime
);

fvScalarMatrix rhoESu
(
    e,
    rho.dimensions()*e.dimensions()*dimVolume/dimTime
);

