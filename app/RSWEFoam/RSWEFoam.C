/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    
Application
    RSWEFoam
    
Description
    Solver for inviscid  shallow water flows governed by 
    regularized shallow water equations (RSWE). 
    The regularization method is based on quasi-gasdynamic approach (QGD).
    
    RSWE has been developed by scientific group from
    Keldysh Institute of Applied Mathematics,
    see http://elizarova.imamod.ru/selection-of-papers.html
    
    A comprehensive description of RSWE equations and their applications
    can be found here:
    \verbatim
    Elizarova, T.G. Bulatov O.V. "Regularized Shallow Water Equations and 
    an Efficient Method for Numerical Simulation of Shallow Water Flows."
    J. Computational Mathematics and Mathematical Physics,  
    2011, Vol. 51, No 1, pp. 160-173.
    \endverbatim
    
    A brief of theory on QGD and QHD system of equations:
    \verbatim
    Elizarova, T.G. and Sheretov, Y.V.
    "Theoretical and numerical analysis of quasi-gasdynamic and quasi-hydrodynamic
    equations"
    J. Computational Mathematics and Mathematical Physics, vol. 41, no. 2, pp 219-234,
    2001
    \endverbatim
    
    Developed by UniCFD group (www.unicfd.ru) of ISP RAS (www.ispras.ru).
    
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "QGD.H"
#include "fvOptions.H"
#include "shallowWaterQGDThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFaceFields.H"
    #include "createFaceFluxes.H"
    #include "createTimeControls.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {

		/*
         *
         * Update fields
         *
         */
        #include "updateFields.H"
        /*
         *
         * Update fluxes
         *
         */
        #include "updateFluxes.H"

        /*
         *
         * Update time step
         *
         */	
		
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"		
        #include "setDeltaT.H"
  
		solve
		(
			fvm::ddt(h)
			+
			fvc::div(phiJm)
		);		  
        
        solve
		(
			fvm::ddt(hU)
			+
			fvc::div(phiJmU)
			+ 
            ghGradB
            -
            magg * tau * fvc::div(hU) * fvc::grad(b)
            +
            magg * n * n * hU * mag(U)/pow(max(h,dimensionedScalar("h0", dimLength, eps0)),4.0/3.0)
			-
			fvc::div(phiPi)
			- 
			NS * divPiNS
		); 
        	
        if (!dryZoneCondition)
        {
            if (gMin(h) <= 0)
            {
                FatalErrorIn("RSWEFoam.C") << "Can't calculate cases h <=0 without dryZoneCondition = true." << nl << exit(FatalError);
            }
            Info << gMin(h) << endl;
            U.ref() =
                hU()
               /h();
            U.correctBoundaryConditions();
        }
        else
        {
            U == hU/max(h,dimensionedScalar("h0", dimLength, eps0));
            forAll(U,celli)
            {
                if (h[celli] <= epsilon[celli]) 
                {
                    U[celli] = vector::zero;
                    hU[celli] = vector::zero;
                }
            }  
            
            forAll(mesh.boundary(), patchID) 
            {
                forAll (mesh.boundary()[patchID],facei) 
                {                 
                    if (h.boundaryField()[patchID][facei] <= epsilon.boundaryField()[patchID][facei]) 
                    {
                        U.boundaryFieldRef()[patchID][facei] = vector::zero;
                        hU.boundaryFieldRef()[patchID][facei] = vector::zero;
                    }
                }
            }        
        }    

		ksi == b + h;       
		
        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
