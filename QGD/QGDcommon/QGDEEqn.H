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
    This file is part of QGDsolver, based on OpenFOAM library.

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

description:
    Equation of energy for QGD solver
source file:
    QGDFoam.C

\*----------------------------------------------------------------------------*/

    fvScalarMatrix EEqn
    (
        fvm::ddt(rhoE)
      + fvc::div(phiJmH)
      + fvc::div(phiQ)
      - fvc::div(phiPiU)
      - fvc::div(phiSigmaDotU)
    );
    
    EEqn.solve();
    
    // Correct energy
    e = rhoE/rho - 0.5*magSqr(U);
    e.correctBoundaryConditions();
    
    // Solve diffusive QGD & NS part
    if (implicitDiffusion)
    {
        solve
        (
            fvm::ddt(rho, e) - fvc::ddt(rho,e)
          - fvm::laplacian(alphauf, e)
          ==
          rhoESu
        );
        
        rhoE = rho*(e + 0.5*magSqr(U));
    }
    else
    {
        solve
        (
            fvm::ddt(rho,e) - fvc::ddt(rhoE)
            ==
            rhoESu
        );
    }
    
    rhoE.boundaryFieldRef() == rho.boundaryField()*
        (e.boundaryField() + 0.5*magSqr(U.boundaryField()));
//
//END-OF-FILE
//


