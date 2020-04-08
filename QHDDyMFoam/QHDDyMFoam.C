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
    QHDDyMFoam

Description
    Solver for unsteady 3D turbulent flow of incompressible fluid governed by
    quasi-hydrodynamic dynamic (QHD) equations with support for deforming
    meshes.

    QHD system of equations has been developed by scientific group from
    Keldysh Institute of Applied Mathematics,
    see http://elizarova.imamod.ru/selection-of-papers.html

    Developed by UniCFD group (www.unicfd.ru) of ISP RAS (www.ispras.ru).


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "QHD.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    bool checkMeshCourantNo
    (
        thermo.subDict("QGD").lookupOrDefault("checkMeshCourantNo", false)
    );
    #include "createFaceFields.H"
    #include "createFaceFluxes.H"
    #include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        checkMeshCourantNo = thermo.subDict("QGD").lookupOrDefault
        (
            "checkMeshCourantNo",
            checkMeshCourantNo
        );

      #include "readTimeControls.H"
      #include "QHDCourantNo.H"
      #include "setDeltaT-QGDQHD.H"

        /*
         *
         * Update time step
         *
         */
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Store old time values
        U.oldTime();
        T.oldTime();

        /*
         *
         * Update the mesh
         *
         */
        mesh.update();

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
        
        #include "QHDpEqn.H"
        
        if (mesh.changing())
        {
            fvc::makeRelative(phi, U);

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }
        
        turbulence->correct();
        
        #include "QHDUEqn.H"
        
        #include "QHDTEqn.H"
        
        if (p.needReference())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pRefValue - getRefCellValue(p, pRefCell)
            );
        }
        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }



    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
