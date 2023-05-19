/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    updateFields
Description
    Updates fields.
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
 *
 * Linear interpolation of fields from volumes to face centers
 *
 */

// Height
hf = qgdInterpolate(h);

// Height^2
h2f = qgdInterpolate(h*h);

// Bathymetry
bf = qgdInterpolate(b);

// Velocity
Uf = qgdInterpolate(U);

// Momentum
hUf = hf*Uf;


// Nonlinear Tensor
UhUf = hf*qgdInterpolate(U*U); 

if (!dryZoneCondition)
{
    if (gMin(h) <= 0)
    {
        FatalErrorIn("RSWEFoam.C") << "Can't calculate cases h <=0 without dryZoneCondition = true." << nl << exit(FatalError);
    }
    Info << gMin(h) << endl;
    tau = alpha*hQGD_/(sqrt(magg*h) + tauU*mag(U));
}
else
{
    tau = alpha*hQGD_/(sqrt(magg*max(h,dimensionedScalar("h0", dimLength, eps0))) + tauU*mag(U));
    forAll(tau, celli)
    {
        if (h[celli] <= eps0) 
        {
            tau[celli] = 0;
        }
    }
    forAll(mesh.boundary(), patchID) 
    {
        forAll (mesh.boundary()[patchID],facei) 
        {         
            if (h.boundaryField()[patchID][facei] <= eps0) 
            {
                tau.boundaryFieldRef()[patchID][facei] = 0;
            }
        }
    }  
} 

tauf = qgdInterpolate(tau);

vector e1(1, 0, 0);
vector e2(0, 1, 0);

volVectorField gradB = fvc::div(bf * mesh.Sf());


phih2 = mesh.Sf() * h2f;


volVectorField gradh2 = fvc::div(phih2);
volVectorField gradh = fvc::div(mesh.Sf() * hf);
volVectorField gradksi = fvc::grad(ksi);


if (wellBalancedScheme)
{  
    
	hSumX = gradh.component(vector::X);
	hhSumX = gradh2.component(vector::X);
    
    hSumY = gradh.component(vector::Y);
	hhSumY = gradh2.component(vector::Y);
    
    scalar smallVal = 1e-15;

	forAll(h_,celli)
	{
		if (mag(hSumX[celli]) > smallVal && h[celli] > epsilon[celli])
		{
			h_[celli].component(vector::X) = hhSumX[celli] / hSumX[celli];           
		}
		else
		{
            h_[celli].component(vector::X) = h[celli];
		}
        if (mag(hSumY[celli]) > smallVal && h[celli] > epsilon[celli])
		{
            h_[celli].component(vector::Y ) = hhSumY[celli] / hSumY[celli];
		}
		else
		{
            h_[celli].component(vector::Y) = h[celli];
		}
	}
    
    ghGradB =  e1 * (0.5 * magg * h_.component(vector::X) * gradksi.component(vector::X) ) + e2 * (0.5 * magg * h_.component(vector::Y) * gradksi.component(vector::Y));
}
else
{
    ghGradB = magg * h * gradB + 0.5 * magg * gradh2;
} 

nuf = tauf * magg * hf;

divPiNS = fvc::laplacian(nuf * hf,U)
        +
		fvc::div(nuf * hf * mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))));
               

// ************************************************************************* //