#include "readGravitationalAcceleration.H"
#include "readShallowWaterProperties.H"

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field h\n" << endl;

volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field hU\n" << endl;
volVectorField hU
(
    IOobject
    (
        "hU",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h*U
);

//bathymetry
volScalarField b
(
    IOobject
    (
        "b",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    h*0
);

volVectorField h_
(
    IOobject
    (
        "h_",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    h*fvc::grad(b)
);

volScalarField ksi
(
    IOobject
    (
        "ksi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h+b
);
volScalarField epsilon
(
    IOobject
    (
        "epsilon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    h*0
);

volVectorField divPiNS
(
    IOobject
    (
        "divPiNS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U*mag(U)
);

volScalarField hQGD_ = thermo.hQGD();

surfaceScalarField hQGDf_("hQGDf", 1.0 / mag(mesh.surfaceInterpolation::deltaCoeffs()));

volScalarField tau
(
    IOobject
    (
        "tau",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE //NO_WRITE
    ),
    alpha*hQGD_/sqrt(magg*max(h,dimensionedScalar("h0", dimLength, 1e-6)))
);
volScalarField hSumX
(
    IOobject
    (
        "hSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hSumX", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 )
);

volScalarField hhSumX
(
    IOobject
    (
        "hhSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hhSumX", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 )
);

volScalarField hSumY
(
    IOobject
    (
        "hSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hSumY", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 )
);

volScalarField hhSumY
(
    IOobject
    (
        "hhSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hhSumY", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 )
);

volScalarField hSumZ
(
    IOobject
    (
        "hSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hSumZ", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0 )
);

volScalarField hhSumZ
(
    IOobject
    (
        "hhSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("hhSumZ", dimensionSet(0, 1, 0, 0, 0, 0, 0), 0.0 )
);

volVectorField ghGradB
(
    IOobject
    (
        "ghgradB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    magg*h*fvc::grad(b)
);

const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();


if (dryZoneCondition)
{
    Info << "Creating cutoff parameter epsilon" << nl;
    
    forAll(epsilon,cellI)
	{
        labelList neighbours = mesh.cellCells()[cellI];
        scalar deltaB = eps0;
        
        forAll(neighbours,cellJ)
        {
            deltaB = max(deltaB, mag(b[neighbours[cellJ]] - b[cellI]));
        }
        
        epsilon[cellI] = kappa * deltaB;
	}
    
    Info << "==========" << nl;
    
    // Loop over all boundary cells
    forAll(mesh.boundary(), patch) 
    {
        scalarField bNei;
        if (isA<processorPolyPatch>(boundaryMesh[patch]))
        {
            bNei = b.boundaryField()[patch].patchNeighbourField();
        }
        
        forAll (mesh.boundary()[patch],facei) 
        {                 
            const label& bCell = boundaryMesh[patch].faceCells()[facei]; 
            
            labelList neighbours = mesh.cellCells()[bCell];
            scalar deltaB = eps0;
            
            
            forAll(neighbours,cellJ)
            {
                deltaB = max(deltaB,mag(b[neighbours[cellJ]] - b[bCell]));
                
                // If boundary is between processors - use b value from neighbour
                if (isA<processorPolyPatch>(boundaryMesh[patch]))
                {
                    deltaB = max(deltaB, mag(bNei[facei] - b[bCell]));
                }
            }
            epsilon[bCell] = kappa * max(epsilon[bCell], deltaB);    
        }
    }   
    epsilon.max(eps0);
}