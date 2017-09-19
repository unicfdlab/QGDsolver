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

//force rhoW to be zero at wall patches
forAll(mesh.boundaryMesh(), patchi)
{
    if (isA<wallPolyPatch>(mesh.boundaryMesh()[patchi]))
    {
        rhoW.boundaryFieldRef()[patchi] = vector::zero;
    }
}

//divergent scheme
//rhoW = tauQGDf * 
//        (
//            (faceStencil.faceVectorGrad(rhoU) &  Uf)
//            + 
//            (rhofUf * divUf)
//            + gradPf 
//        );


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
//    + //viscous fluxes (QGD + NS) are now implicit
//    muf*
//    (
//        gradUf
//        +
//        Foam::T(gradUf)
//        -
//        (2.0/3.0)*I*divUf
//    );

phiPi = Pif & mesh.Sf();

//viscous stress explicit terms interpolation
if (!inviscid)
{
//    tauMCPtr() = muf*(Foam::T(gradUf) - (2.0/3.0)*I*divUf);
    tauMCPtr() = muf*linearInterpolate(Foam::T(fvc::grad(U)) - (2.0/3.0)*I*fvc::div(U));
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
//    //NS and QGD heat flux are now implicit
//    -
//    muf*faceStencil.faceScalarGrad(e);

phiQ = qf & mesh.Sf();

phiPiU = (Pif & Uf) & mesh.Sf();

//viscous heat
if (!inviscid)
{
//    sigmaDotUPtr() = (muf*gradUf + tauMCPtr()) & Uf;
    sigmaDotUPtr() = (muf*linearInterpolate(fvc::grad(U)) + tauMCPtr()) & Uf;
    phiSigmaDotU = sigmaDotUPtr() & mesh.Sf();
}
//---------End---------
