
/*
  Field for debugging parallel execution
  and tau-terms evaluation
*/
/*
volScalarField cellNo
(
    IOobject
    (
        "cellNo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*0.0
);

forAll(cellNo, celli)
{
    cellNo[celli] = mesh.C()[celli].z();
}
forAll(cellNo.boundaryField(), patchi)
{
    forAll(cellNo.boundaryField()[patchi], facei)
    {
        cellNo.boundaryFieldRef()[patchi][facei] = mesh.Cf().boundaryField()[patchi][facei].z();
    }
}
surfaceVectorField cellNoGradf
(
    "cellNoGradf", fvsc::grad(cellNo)
);

cellNo.write();
cellNoGradf.write();
return 0;
*/
//Gradients and divergence
//---------Start---------


surfaceVectorField gradAlphalf
(
    "gradTf", fvsc::grad(alphal)
);

surfaceTensorField gradUf
(
    "gradUf",
    fvsc::grad(U)
);

surfaceVectorField gradpf
(
    "gradpf",
    fvsc::grad(p)
);

surfaceScalarField divUf
(
    "divUf",
    tr(gradUf)
);

//---------End---------
//P equation fluxes
//---------Start---------
//---------End---------

// Fluxes for momentum balance equation
//---------Start---------
surfaceVectorField phiUfRhof
(
    "phiUf",
    phi * Uf * rhof
);

surfaceScalarField phiAlphalf
(
    "phiTf",
    phi*alphalf
);

surfaceVectorField phiPi
(
    "phiPi",
    muf/rhof * mesh.Sf() & (gradUf+gradUf.T())
);

//---------End---------

// Fluxes for T balance equation

surfaceScalarField Hif
(
    "Hif",
    alphaf/rhof
);
Info << "Creating Face Fluxes is done"<< endl;

//---------End---------
