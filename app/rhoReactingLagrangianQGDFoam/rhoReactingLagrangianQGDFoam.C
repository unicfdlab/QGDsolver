
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
    reactinLagrangianQGDFoam

Description
    Solver for unsteady 3D turbulent flow of perfect gases mixture governed by
    quasi-gas dynamic (QGD) equations at all Mach numbers (from 0 to
    infinity) coupled with chemistry reactions.

    QGD system of equations has been developed by scientific group from
    Keldysh Institute of Applied Mathematics,
    see http://elizarova.imamod.ru/selection-of-papers.html

    Developed by UniCFD group (www.unicfd.ru) of ISP RAS (www.ispras.ru).


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mcQGD.H"
#include "fvOptions.H"
#include "turbulentFluidThermoModel.H"
#include "basicReactingCloud.H"
#include "CombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    Info << "Basic fields were created" << endl;
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
        #include "setDeltaT-QGDQHD.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        parcels.evolve();

        // --- Store old time values
        rho.oldTime();
        rhoU.oldTime();
        U.oldTime();
        rhoE.oldTime();
        e.oldTime();
        forAll(Y,i)
        {
            Y[i].oldTime();
        }

        // --- Solve density
        rhoSu = parcels.Srho(rho);
        #include "QGDRhoEqn.H"

        // --- Solve for mass fractions
        #include "QGDYEqn.H"

        // --- Solve momentum
        rhoUSu = parcels.SU(U);
        #include "QGDUEqn.H"

        //--- Solve energy
        rhoESu = parcels.Sh(e) + Qdot;
        #include "addEnergyFluxes.H"
        #include "QGDEEqn.H"

        if ( (min(T).value() <= 0.0) || (min(rho).value() <= 0.0) )
        {
            U.write();
            T.write();
            rho.write();
            forAll(Y, i)
                Y[i].write();
            p.write();
        }

        thermo.correct();

        // Correct pressure
        //p.ref() =
        //    rho()
        //   /psi();
        //p.correctBoundaryConditions();
        //rho.boundaryFieldRef() = psi.boundaryField()*p.boundaryField();
        thermo.findPressure(rho);
        rho.boundaryFieldRef() = thermo.rho()().boundaryField();

        runTime.write();
        if (runTime.outputTime())
        {
            turbulence->alphaEff()().write();
            thermo.c().write();
            thermo.psi().write();
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
