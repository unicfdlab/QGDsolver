/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    QGDFoam

Description
    Solver for unsteady 2D (3D is under development) turbulent flow of perfect gas governed by
    quasi-gas dynamic (QGD) equations at high Mach numbers (from 2 to
    infinity).
    
    QGD system of equations has been developed by scientific group from
    Keldysh Institute of Applied Mathematics, 
    see http://elizarova.imamod.ru/selection-of-papers.html
    
    A comprehensive description of QGD equations and their applications can be found here:
    \verbatim
    Elizarova, T.G.
    "Quasi-Gas Dynamic equations"
    Springer, 2009
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
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFaceFields.H"
    #include "createFaceFluxes.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
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
        #include "readTimeControls.H"
        #include "QGDCourantNo.H"
        #include "setDeltaT.H"
        
        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- Store old time values
        U.oldTime();
        T.oldTime();
         	
        //Continuity equation
        fvScalarMatrix pEqn
        (
             fvc::div(phiu)
	    +fvc::div(phiwo)
            -fvm::laplacian(taubyrhof,p)
        );

	pEqn.setReference(pRefCell, pRefValue);
        
        pEqn.solve();

	Info << "Solve of continuity finished" << endl;
        
        phi = phiu - phiwo + pEqn.flux();
        
	gradPf = fvsc::grad(p);
        
	Wf = tauQGDf*((Uf & gradUf) + gradPf/rhof + beta*g*Tf);
   
	phiUf = (phi * Uf) - (mesh.Sf() & (Uf * Wf));

      	// --- Solve U
        solve
        (
            fvm::ddt(U)
            + 
            fvc::div(phiUf)
            +
            fvc::grad(p)/rho
            -
            //fvc::div(phiPi)
	    fvm::laplacian(muf/rhof,U)
	    -
	    fvc::div(muf/rhof * mesh.Sf() & linearInterpolate(Foam::T(fvc::grad(U))))
            +
            BdFrc
        );
       
	phiTf = phi * Tf;
        
        // --- Solve T
        solve
        (
            fvm::ddt(T)
          + fvc::div(phiTf)
          - fvc::laplacian(Hif,T)
        );      

        thermo.correct();
      
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
        if (runTime.outputTime())
        {
            tauQGD.write();
        }
        
        Info<< "max/min T:    "<< max(T).value()  << "/" << min(T).value()   << endl;
        Info<< "max/min p:    "<< max(p).value()  << "/" << min(p).value()   << endl;
        Info<< "max/min rho:  "<< max(rho).value()<< "/" << min(rho).value() << endl;
        Info<< "max/min U:    "<< max(U).value()  << "/" << min(U).value()   << endl;
    }
    
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
