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
    sQGDFoam

Description
    Solver for unsteady 3D turbulent flow of perfect gas governed by
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
#include "fvOptions.H"
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
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        /*
         *
         * Update QGD viscosity
         *
         */
        turbulence->correct();
        
        
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
        rho.oldTime();
        rhoU.oldTime();
        U.oldTime();
        rhoE.oldTime();
        e.oldTime();
        
        // --- Solve density
        solve
        (
            fvm::ddt(rho)
            +
            fvc::div(phiJm)
        );
        
        // --- Solve momentum
        solve
        (
            fvm::ddt(rhoU)
            + 
            fvc::div(phiJmU)
            +
            fvc::div(phiP)
            -
            fvc::div(phiPi)
        );
        
        // Correct velocity
        U.ref() =
            rhoU()
           /rho();
        U.correctBoundaryConditions();
        
        // Solve diffusive QGD & NS part
        if (implicitDiffusion)
        {
            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U) - fvc::ddt(rho,U)
              - fvm::laplacian(muf, U)
              - fvc::div(phiTauMC)
            );
            
            solve
            (
                UEqn
            );
            
            rhoU = rho*U;
            
            sigmaDotUPtr() = (muf*linearInterpolate(fvc::grad(U)) + tauMCPtr()) & Uf;
            
            phiSigmaDotU = mesh.Sf() & sigmaDotUPtr(); //or eqn.flux()?
        }
        rhoU.boundaryFieldRef() == rho.boundaryField()*
            U.boundaryField();
        
        //--- Solve energy
        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiJmH)
          + fvc::div(phiQ)
          - fvc::div(phiPiU)
          - fvc::div(phiSigmaDotU)
        );
        
        // Correct energy
        e = rhoE/rho - 0.5*magSqr(U);
        fvOptions.correct(e);
        e.correctBoundaryConditions();
        
        // Solve diffusive QGD & NS part
        if (implicitDiffusion)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho,e)
              - fvm::laplacian(alphauf, e)
            );
            
            rhoE = rho*(e + 0.5*magSqr(U));
        }
        rhoE.boundaryFieldRef() == rho.boundaryField()*
            (e.boundaryField() + 0.5*magSqr(U.boundaryField()));
        
        if ( (min(e).value() <= 0.0) || (min(rho).value() <= 0.0) )
        {
            U.write();
            e.write();
            rho.write();
        }
        
        thermo.correct();
        
        // Correct pressure
        p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField();
        
        runTime.write();
        
        if (runTime.outputTime())
        {
            e.write();
        }
        
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
