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
    Solver for unsteady 3D turbulent flow of 2 incompressible immisicible 
    fluids governed by quasi-hydrodynamic dynamic (QHD) equations.

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
#include "wallFvPatch.H"
#include "fvsc.H"
#include "twoPhaseIcoQGDThermo.H"
#include "QGDInterpolate.H"
#include "MULES.H"

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
        //#include "QHDCourantNo.H"
        if(runTime.controlDict().lookupOrDefault<bool>("adjustTimeStep", false))
        {
            surfaceScalarField Unf
            (
                "Unf",
                Uf & mesh.Sf() / mesh.magSf()
            );
            
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
        //#include "setDeltaT-QGDQHD.H"
        if (adjustTimeStep)
        {
            scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
            scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);
            scalar maxDeltaT1 = min(alpha1f*Tau1+alpha2f*Tau2).value();
            maxDeltaT1 = min(maxDeltaT,maxDeltaT1)*0.75;
            
            runTime.setDeltaT
            (
                min
                (
                    deltaTFact*runTime.deltaTValue(),
                    maxDeltaT1
                )
            );
            
            Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
        }

        runTime++;
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- Store old time values
        U.oldTime();
        
        //Continuity equation
        phiwm = -phiwo1*alpha1f-phiwo2*alpha2f;
        surfaceScalarField tphi = phiu*da1dtf*(Tau1-Tau2);
        tphi.setOriented(true);
        phiwm += tphi;
        coeffp = alpha1f*Tau1/rho1 + alpha2f*Tau2/rho2;
        p.correctBoundaryConditions();
        
        //solve for pressure
        {
            fvScalarMatrix pEqn1
            (
                -fvm::laplacian(Tau1/rho1,p)
            );
            fvScalarMatrix pEqn2
            (
                -fvm::laplacian(Tau2/rho2,p)
            );
            fvScalarMatrix pEqn
            (
                 fvc::div(phiu)
                +fvc::div(phiwm)
                -fvm::laplacian(coeffp,p)
            );
            pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));
            pEqn.solve();
            
            phiw1 = phiwo1 - pEqn1.flux();
            phiw2 = phiwo2 - pEqn2.flux();
            phi1 = phiu - phiw1;
            phi2 = phiu - phiw2;

            surfaceScalarField phiq1 = phiu*da1dtf*Tau1;
            phiq1.setOriented(true);
            surfaceScalarField phiq2 = phiu*da1dtf*Tau2;
            phiq2.setOriented(true);
            phi = phiu+phiwm+pEqn.flux();
        }
        
        gradpf = fvsc::grad(p);
        
        W1 = ((Uf & gradUf) + (1./rho1)*gradpf - g - cFrcf/rho1)*Tau1;
        W2 = ((Uf & gradUf) + (1./rho2)*gradpf - g - cFrcf/rho2)*Tau2;
        phiWr =
            fvc::flux(-phiw2+phiw1,(1.0 - alpha1),"div(phi,alpha1)");
        
        surfaceScalarField phic(thermo.cAlpha()*mag(phi/mesh.magSf()));
        // Do not compress interface at non-coupled boundary faces
        // (inlets, outlets etc.)
        {
            surfaceScalarField::Boundary& phicBf =
                phic.boundaryFieldRef();
            
            forAll(phic.boundaryField(), patchi)
            {
                fvsPatchScalarField& phicp = phicBf[patchi];
            
                if (!phicp.coupled())
                {
                    phicp == 0;
                }
            }
        }
        surfaceScalarField phir(phic*thermo.nHatf());
        
        phiAlpha1f = 
            fvc::flux
            (
                phi,
                alpha1,
                "div(phi,alpha1)"
            )
            +
            fvc::flux
            (
                -fvc::flux(-phir, alpha2, "div(phir,alphar)"),
                alpha1,
                "div(phir,alphar)"
            );
        
        // --- Solve for volume fraction
        {
            //limit flux with MULES limiter
            MULES::limit
            (
                1.0 / runTime.deltaTValue(),
                geometricOneField(),
                alpha1,
                phi,
                phiAlpha1f,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField(),
                false //return total flux
            );
            
            //add qgd terms from Wr
            phiAlpha1f+=
                fvc::flux
                (
                    -phiWr,
                    alpha1,
                    "div(phi,alpha1)"
                );
                //+
                //phiq1;
                //-
                //fvc::flux
                //(
                //    tphi,
                //    alpha1,
                //    "div(phi,alpha1)"
                //);
            solve
            (
                fvm::ddt(alpha1)
              + fvc::div(phiAlpha1f)
            );
            //
            Info << "max/min alpha1: " << max(alpha1).value() << "/" << min(alpha1).value() << endl;
            alpha1 = max(min(alpha1,1.0),0.0);
            alpha2 = 1.0 - alpha1;
        }
        
        phiAlpha2f = phi - phiAlpha1f;
        rhoPhi = phiAlpha1f*rho1 + phiAlpha2f*rho2;
        
        phiRhofWf = phiu*(alpha1f*rho1*W1 + alpha2f*rho2*W2);
        phiRhofWf.setOriented(true);
        phiUfRhof = qgdFlux
        (
            rhoPhi,
            U,
            Uf
        )
        -
        phiRhofWf;
        
        // --- Solve U
        solve
        (
            fvm::ddt(rhoU)
            +
            fvc::div(phiUfRhof)
            +
            fvc::grad(p)
            *(1.0 + da1dt*(Tau1-Tau2))
            -
            fvc::laplacian(muf,U)
            -
            fvc::div(muf*mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))))
            -
            BdFrc
            -
            cFrc*(1.0 + da1dt*(Tau1-Tau2))
        );
        
        rho = thermo.rho();
        U = rhoU / rho;
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() = rho.boundaryField()*U.boundaryField();
	
        /*
         *
         * Update physical properties
         *
         */
        thermo.correct(); //curvature, mechanics and so on
        //turbulence.correct();

        if(runTime.write())
        {
            phi.write();
        }
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
	
    }
    
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
