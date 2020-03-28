
fvVectorMatrix USu
(
    U,
    U.dimensions()*dimVolume/dimTime
);

fvScalarMatrix TSu
(
    T,
    T.dimensions()*dimVolume/dimTime
);

