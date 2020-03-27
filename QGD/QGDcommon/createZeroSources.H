
fvScalarMatrix rhoSu
(
    rho,
    rho.dimensions()/dimTime
);

fvVectorMatrix rhoUSu
(
    U,
    rho.dimensions()*U.dimensions()/dimTime
);

fvScalarMatrix rhoESu
(
    e,
    rho.dimensions()*e.dimensions()/dimTime
);

