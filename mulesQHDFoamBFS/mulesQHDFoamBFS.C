/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
                QGDsolver   | Copyright (C) 2016-2018 ISP RAS (www.unicfd.ru)
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
    QHDFoam

Description
    Solver for unsteady 3D turbulent flow of incompressible fluid governed by
    quasi-hydrodynamic dynamic (QHD) equations.

    QHD system of equations has been developed by scientific group from
    Keldysh Institute of Applied Mathematics,
    see http://elizarova.imamod.ru/selection-of-papers.html

    A comprehensive description of QGD equations and their applications
    can be found here:
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
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "QHD.H"
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

    pimpleControl pimple(mesh);


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
        #include "QHDCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Store old time values
        U.oldTime();
        T.oldTime();
        turbulence->correct();
        //Continuity equation
	
        while (pimple.correctNonOrthogonal())
        {
                // Pressure corrector

            fvScalarMatrix pEqn
            (
                 fvc::div(phiu)
	        -fvc::div(phiwo)
                -fvm::laplacian(taubyrhof,p_rgh)
            );
	    
            pEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

            pEqn.solve();

                if (pimple.finalNonOrthogonalIter())
                {
                    phi = phiu - phiwo + pEqn.flux();
                }
        }
        
        gradPf = fvsc::grad(p_rgh);
        
        Wf = tauQGDf*((Uf & gradUf) + gradPf/rhof - BdFrcf);
        
        phiUf = (phi * Uf) - (mesh.Sf() & (Uf * Wf));

        // --- Solve U
        solve
        (
            fvm::ddt(U)
            +
            fvc::div(phiUf)
            -
            fvc::laplacian(muf/rhof,U)
            -
            fvc::div(muf/rhof * mesh.Sf() & linearInterpolate(Foam::T(fvc::grad(U))))
            ==
            -
            fvc::grad(p_rgh)/rho
            +
            BdFrc
        );
        
        phiTf = phi * Tf;

        // --- Solve T
        #include "alphaControls.H"
        #include "MULESTEqn.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        if (runTime.outputTime())
        {
            thermo.tauQGD().write();
        }

//        p = p_rgh + rho *(1-beta*T) * gh;
//        
        if (p.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pRefCell)
            );
            p_rgh = p - rho * (1-beta*T) * gh;
        }

        Info<< "max/min T:    "<< max(T).value()  << "/" << min(T).value()   << endl;
        Info<< "max/min p:    "<< max(p).value()  << "/" << min(p).value()   << endl;
        Info<< "max/min rho:  "<< max(rho).value()<< "/" << min(rho).value() << endl;
        Info<< "max/min U:    "<< max(U).value()  << "/" << min(U).value()   << endl;
        
        if(runTime.write())
        {
            phi.write();
            Wf.write();
            BdFrc.write();
            BdFrcf.write();
            phiwo.write();
            phiu.write();
            phi.write();
        }

    }


    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
