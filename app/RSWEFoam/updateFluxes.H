/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd .
    Copyright (C) 2016-2019 ISP RAS (www.ispras.ru) UniCFD Group (www.unicfd.ru)
-------------------------------------------------------------------------------
License
    This file is part of QGDsolver library, based on OpenFOAM+.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
Global
    updateFluxes
Description
    Updates fluxes for continuity equation, momentum balance equation,
    energy balance equation.
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//Gradients and divergence
//---------Start---------
gradUf = fvsc::grad(U);

divUf  = tr(gradUf);

divhUf = fvsc::div(hU);

gradksif = fvsc::grad(ksi);

gradhf = fvsc::grad(h);

surfaceVectorField gradH2 = fvsc::grad(h*h);

//---------End---------

//Continuity equation fluxes
hW = tauf *
        (
            ((Uf * gradhf) & Uf)
            +
            (hUf * divUf)
            +
            (hUf & gradUf)
            +
            magg * hf * gradbf
            +
            magg * 0.5 * gradH2
        );

jm = hUf - hW;

phiJm = mesh.Sf() & jm;


//---------End---------
		
// Fluxes for momentum balance equation
//---------Start---------


if (dryZoneCondition)
{
    forAll(faces, faceI)
    { 
        if (mesh.isInternalFace(faceI))
        {
            if (h[faceOwner[faceI]] <= eps0 && phiJm[faceI] > 0.0 || h[faceNeighbour[faceI]] <= eps0 && phiJm[faceI] < 0.0)
            {
                phiJm[faceI] = 0.0;
            }
            if (h[faceOwner[faceI]] <= eps0 || h[faceNeighbour[faceI]] <= eps0)
            {
                tauf[faceI] = 0.0;
            }
        }
        else
        {
            const label& patch = boundaryMesh.whichPatch(faceI);       // Boundary patch index
            const label& facei = boundaryMesh[patch].whichFace(faceI); // Local boundary face index      
            
            if (isA<processorPolyPatch>(boundaryMesh[patch]))
            {
                scalarField hOwn = h.boundaryField()[patch].patchInternalField();
                scalarField hNei = h.boundaryField()[patch].patchNeighbourField();
                
                scalarField epsOwn = epsilon.boundaryField()[patch].patchInternalField();
                scalarField epsNei = epsilon.boundaryField()[patch].patchNeighbourField();

                if (hOwn[facei] <= eps0 && phiJm.boundaryField()[patch][facei] > 0.0 || hNei[facei] <= eps0 && phiJm.boundaryField()[patch][facei] < 0.0)
                {
                    phiJm.boundaryFieldRef()[patch][facei] = 0.0;
                }
                if (hOwn[facei] <= eps0 || hNei[facei] <= eps0)
                {
                    tauf.boundaryFieldRef()[patch][facei] = 0.0;
                }
            }
        }
    }
}

Pif = 
    //QGD diffusive fluxes
	tauf*
    (
        (UhUf & gradUf)
        +
        magg * Uf *0.5 * gradH2
        +
        magg * Uf * hf * gradbf
    )
    +
    magg * tauf * 
    (
        I * ( hf * divhUf)
    );   

phiPi = mesh.Sf() & Pif;

phiJmU = qgdFlux(phiJm,U,Uf);

//---------End---------	
   		
// ************************************************************************* //