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
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("smoothAlpha");
    argList::addOption("nSmoothIters");
    argList::addOption("smoothCoeff");
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFaceFields.H"
    #include "createFaceFluxes.H"
    #include "createTimeControls.H"
    
    #include "smoothSolution.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        /*
         *
         * Update physical properties
         *
         */
        thermo.correct(); 
        //turbulence.correct();
        
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
        #include "setDeltaT-QGDQHD.H"
        
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- Store old time values
        U.oldTime();
        alpha1.oldTime();
        rho.oldTime();
        
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

            phi = phiu+phiwm+pEqn.flux();
        }
        
        gradpf = fvsc::grad(p);
        //W1 = ((Uf & gradUf) + (1./rho1)*gradpf - g)*Tau1;
        //W2 = ((Uf & gradUf) + (1./rho2)*gradpf - g)*Tau2;
        W1 = ((Uf & gradUf) + (1./rho1)*gradpf - g - linearInterpolate(cFrc)/rho1/*cFrcf/rho1*/)*Tau1;
        W2 = ((Uf & gradUf) + (1./rho2)*gradpf - g - linearInterpolate(cFrc)/rho2/*cFrcf/rho2*/)*Tau2;
        
        phiWr =
            qgdFlux
            (
                -phiw2+phiw1,
                alpha2,
                alpha2f,
                "div(phi,alpha1)"
            );
        phiAlpha1f = 
            qgdFlux
            (
                phi,
                alpha1,
                alpha1f,
                "div(phi,alpha1)"
            )
            +
            //add qgd terms from Wr
            qgdFlux
            (
                -phiWr,
                alpha1,
                alpha1f,
                "div(phi,alpha1)"
            );
        
        //add terms from da1dt
        {
            surfaceScalarField DeltaTauFlux =
                phiu*da1dtf*(Tau1 - alpha1f*(Tau1-Tau2));
                //phiu*da1dtf*Tau1;
            DeltaTauFlux.setOriented(true);
            phiAlpha1f += DeltaTauFlux;
        }
        
        //apply compressive fluxes and MULES limiter
        if (thermo.cAlpha() > SMALL)
        {
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
            
            phiAlpha1f +=
                fvc::flux
                (
                    -fvc::flux(-phir, alpha2, "div(phir,alphar)"),
                    alpha1,
                    "div(phir,alphar)"
                );
            
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
        }
        
        // --- Solve for volume fraction
        {

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
        
        // --- Update density field 
        rho = thermo.rho();
        
        // --- Update mass fluxes
        {
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
        }
        
        // --- Solve U
        if (implicitDiffusion)
        {
            solve
            (
                fvm::ddt(rho,U)
                +
                fvc::div(phiUfRhof)
                -
                fvm::laplacian(muf,U)
                -
                fvc::div(muf*mesh.Sf() & qgdInterpolate(Foam::T(fvc::grad(U))))
                -
                BdFrc
                //+
                //(
                //    fvc::grad(p)
                //    -
                //    cFrc
                //)*(1.0 + da1dt*(Tau1-Tau2))
                +
                (
                    fvc::reconstruct
                    (
                        fvc::snGrad(p)*mesh.magSf()
                    )
                    -
                    cFrc
                )*(1.0 + da1dt*(Tau1-Tau2))
            );
        }
        else
        {
            solve
            (
                fvm::ddt(rho,U)
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
