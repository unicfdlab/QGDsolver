//Gradients and divergence
//---------Start---------
gradPf = faceStencil.faceScalarGrad(p);

gradUf = faceStencil.faceVectorGrad(U);

divUf = faceStencil.faceVectorDiv(U);
//---------End---------

//Continuity equation fluxes
//---------Start---------
phivf = Uf & mesh.Sf();

rhoW = tauQGDf * (faceStencil.faceTensorDiv(rhoUU) + gradPf);

jm = rhofUf - rhoW;

phiJm = jm & mesh.Sf();
phi = phiJm;
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------

phiJmU = (jm * Uf) & mesh.Sf();

phiP = pf*mesh.Sf();

Pif = 
    //QGD diffusive fluxes
    tauQGDf * 
    (
        Uf * (rhof * (Uf & gradUf) + gradPf)
        +
        I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
    );

phiPi = Pif & mesh.Sf();

//viscous stress explicit terms interpolation
if (!inviscid)
{
    tauMCPtr() = muf*(Foam::T(gradUf) - (2.0/3.0)*I*divUf);
//    tauMCPtr() = muf*linearInterpolate(Foam::T(fvc::grad(U)) - (2.0/3.0)*I*fvc::div(U));
    phiTauMC = tauMCPtr() & mesh.Sf();
}
//---------End---------

// Fluxes for energy balance equation
//---------Start---------
phiJmH = phiJm * Hf;

qf =
    -
    tauQGDf*rhofUf*
    (
        (Uf & faceStencil.faceScalarGrad(e))
        +
        (pf * (Uf & faceStencil.faceScalarGrad(1/rho)))
    );

phiQ = qf & mesh.Sf();

phiPiU = (Pif & Uf) & mesh.Sf();

if (!inviscid)
{
    sigmaDotUPtr() = (muf*gradUf + tauMCPtr()) & Uf;
//    sigmaDotUPtr() = (muf*linearInterpolate(fvc::grad(U)) + tauMCPtr()) & Uf;
    phiSigmaDotU = sigmaDotUPtr() & mesh.Sf();
}
//---------End---------
