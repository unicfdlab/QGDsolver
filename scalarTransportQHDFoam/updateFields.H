gradTf = fvsc::grad(T);

alphaf = qgdInterpolate(thermo.alpha());

rhof   = qgdInterpolate(rho);
    
Tf     = qgdInterpolate(T);

Hif    = alphaf/rhof;
