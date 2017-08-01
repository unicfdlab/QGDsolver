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
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "faceGrad/extendedFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    extendedFaceStencil test
    (
          IOobject(
            "test",
            runTime.timeName(),
            mesh,
            regIOobject::NO_READ,
            regIOobject::NO_WRITE),
          false
    );

    IOdictionary paramsQGDDict
    (
        IOobject
        (
            "paramsQGD",            // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    
    scalar ScQGD(readScalar(paramsQGDDict.lookup("ScQGD")));
    scalar PrQGD(readScalar(paramsQGDDict.lookup("PrQGD")));

    surfaceScalarField alphaQGD
    (
        IOobject
        (
            "alphaQGD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        /*
         *
         * Auxiliary fields
         *
         */
        
        // Inversed compressibility
        volScalarField rPsi
        (
            "rPsi", 
            1.0/psi
        );
        
        // Heat capacities ratio
        volScalarField gamma
        (
            "gamma",
            thermo.Cp()/thermo.Cv()
        );
        
        // Speed of sound
        volScalarField c
        (
            "c",
            sqrt(gamma * rPsi)
        );
        
        /*
         *
         * Linear interpolation of fields from volumes to face centers
         *
         */
        
        // Density
        surfaceScalarField rhof
        (
            "rhof",
            linearInterpolate(rho)
        );

        // Velocity
        surfaceVectorField Uf
        (
            "Uf",
            linearInterpolate(U)
        );
        
        // Pressure
        surfaceScalarField pf
        (
            "pf",
            linearInterpolate(p)
        );

        // Heat capacities ratio
        surfaceScalarField gammaf 
        (
            "gammaf",
            linearInterpolate(gamma)
        );
        
        surfaceScalarField gammam1f
        (
            "gammam1",
            gammaf - 1
        );

        // Speed of sound
        surfaceScalarField cf
        (
            "cf",linearInterpolate(c)
        );

        // Heat capacity at constant pressure
        surfaceScalarField Cpf
        (
            "Cpf",
            linearInterpolate(thermo.Cp())
        );
        
        surfaceScalarField Hf
        (
            "Hf",
            linearInterpolate((rhoE + p)/rho)
            //or (rhoEf + pf)/rhof?
        );

        /*
         *
         * QGD coefficients
         *
         */
        surfaceScalarField hQGD
        (
            "hQGD", 
            1.0 / mesh.surfaceInterpolation::deltaCoeffs()
        );
        
        surfaceScalarField tauQGD
        (
            "tauQGD",
            alphaQGD * hQGD / cf
        );
        
        surfaceScalarField muQGD
        (
            "muQGD",
            tauQGD*pf*ScQGD
        );

        surfaceScalarField kappaQGD
        (
            "kappaQGDf",
            muQGD*Cpf / PrQGD
        );

        /*
         *
         * Fluxes
         *
         */
        
        //Gradients and divergence
        //---------Start---------
        surfaceVectorField gradPf 
        (
            "gradPf", test.faceScalarGrad(p)
        );
        
        surfaceTensorField gradUf
        (
            "gradUf",
            test.faceVectorGrad(U)
        );
        
        surfaceScalarField divUf
        (
            "divUf",
            test.faceVectorDiv(U)
        );
        //---------End---------
        
        //Continuity equation fluxes
        //---------Start---------
        surfaceScalarField phivf
        (
            "phivf",
            Uf & mesh.Sf()
        );
        
        phi = phivf*rhof;
        
        surfaceVectorField rhoW1
        (
            "rhoW1",
            tauQGD * test.faceTensorDiv(rho * (U * U))
        );
        
        surfaceScalarField phiRhoW1
        (
            "phiRhoW1", rhoW1 & mesh.Sf()
        );
        
        surfaceVectorField rhoW2
        (
            "rhoW2",
            tauQGD * gradPf
        );
        
        surfaceScalarField phiRhoW2
        (
            "phiRhoW2",
            rhoW2 & mesh.Sf()
        );
        
        surfaceVectorField jm
        (
            "jm",
            Uf*rhof - rhoW1 - rhoW2
        );
        
        surfaceScalarField phiJm
        (
            "phiJm",
            jm & mesh.Sf()
        );
        //---------End---------
        
        // Fluxes for momentum balance equation
        //---------Start---------
        surfaceVectorField phiJmU
        (
            "phiJmU",
            (jm * Uf) & mesh.Sf()
        );
        
        surfaceVectorField phiP
        (
            "phiP",
            pf*mesh.Sf()
        );
        
        surfaceTensorField Pif
        (
            "Pif",
            //QGD diffusive fluxes
            tauQGD * 
            (
                Uf * (rhof * (Uf & gradUf) + gradPf)
                +
                I * ( (Uf & gradPf) + (gammaf * pf * divUf) )
            )
            +
            //NS diffusive fluxes
            (
                muQGD*(gradUf + gradUf.T())
                -
                muQGD*I*divUf
            )
        );
        
        surfaceVectorField phiPi
        (
            "phiPi",
            Pif & mesh.Sf()
        );
        
        //---------End---------

        // Fluxes for energy balance equation
        //---------Start---------
        surfaceScalarField phiJmH
        (
            "phiJmH",
            phiJm * Hf
        );
        
        surfaceVectorField qf
        (
            "qf",
            -
            kappaQGD*test.faceScalarGrad(T)
            -
            tauQGD* 
            ( 
                ((fvc::interpolate(rho*(U*U))) & test.faceScalarGrad(e))
                +
                (pf * rhof * Uf * (Uf & test.faceScalarGrad(1/rho)))
            )
        );
        
        surfaceScalarField phiQ
        (
            "phiQ",
            qf & mesh.Sf()
        );
        
        surfaceScalarField phiPiU
        (
            "phiPiU",
            (Pif & Uf) & mesh.Sf()
        );
        
        
        // End for third equation
        // ******************************************************************* //

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
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        //--- Solve energy
        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiJmH)
          + fvc::div(phiQ)
          - fvc::div(phiPiU)
        );
        
        // Correct energy
        e = rhoE/rho - 0.5*magSqr(U);
        e.correctBoundaryConditions();
        rhoE.boundaryFieldRef() == rho.boundaryField()*
        (e.boundaryField() + 0.5*magSqr(U.boundaryField()));
        
        thermo.correct();
        
        // Correct pressure
        p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
        
        Info<< "max/min T:    "<< max(T).value()  << "/" << min(T).value()   << endl;
        Info<< "max/min p:    "<< max(p).value()  << "/" << min(p).value()   << endl;
        Info<< "max/min rho:  "<< max(rho).value()<< "/" << min(rho).value() << endl;
        Info<< "max/min U:    "<< max(U).value()  << "/" << min(U).value()   << endl;
    }



    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
