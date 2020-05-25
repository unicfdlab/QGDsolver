

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
    scalarTransportQHDFoam

Description
    Evolves a passive scalar transport equation with regularized terms.

    Developed by UniCFD group (www.unicfd.ru) of ISP RAS (www.ispras.ru).


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "QHD.H"

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

    //solve for pressure to get zero divergence flux field
    #include "updateFluxes.H"
    
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
         * Update time step
         *
         */
        #include "readTimeControls.H"
        if(runTime.controlDict().lookupOrDefault<bool>("adjustTimeStep", false))
        {
            surfaceScalarField Cof
            (
                "Cof",
                runTime.deltaT()*
                (
                    mag(Uf)/hQGDf
                )
            );
            CoNum = max(Cof).value();
            Info << "Courant Number = " << CoNum << endl;
        }
        
        #include "setDeltaT-QGDQHD.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Store old time values
        T.oldTime();
        
        {
            phiTf = qgdFlux(phiu,T,Tf);
            surfaceScalarField phiTauTReg = tauQGDf*phiu*(Uf & gradTf);
            
            // --- Solve T
            if (implicitDiffusion)
            {
                solve
                (
                    fvm::ddt(T)
                    + fvc::div(phiTf) - fvc::Sp(fvc::div(phiu),T)
                    - fvm::laplacian(Hif,T)
                    - fvc::div(phiTauTReg)
                    ==
                    TSu
                );
            }
            
            Info << "max/min of T: " << max(T).value() << "/" << min(T).value() << endl;
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
