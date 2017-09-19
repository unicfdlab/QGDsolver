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
    Quasi-Gasdynamic solver.
    Now only for orthogonal mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiQGDThermo.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "faceGrad/extendedFaceStencil.H"
#include "limitedSurfaceInterpolationScheme.H"
#include "wallFvPatch.H"

#warning "insert QGD includes in separate file"

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
        rhoU.boundaryFieldRef() == rho.boundaryField()*
            U.boundaryField();
        
        // Solve diffusive QGD & NS part
        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muf, U)
              - fvc::div(phiTauMC)
            );
            rhoU = rho*U;
        }

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
        e.correctBoundaryConditions();
        thermo.correct();
        rhoE.boundaryFieldRef() == rho.boundaryField()*
            (e.boundaryField() + 0.5*magSqr(U.boundaryField()));
        
        // Solve diffusive QGD & NS part
        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), e)
            );
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }
        
        // Correct pressure
        p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField();
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
        if (runTime.outputTime())
        {
//            phiJm.write();
//            rhoW.write();
//            phiJmU.write();
//            phiP.write();
//            phiPi.write();
//            gradPf.write();
//            tauQGDf.write();
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
