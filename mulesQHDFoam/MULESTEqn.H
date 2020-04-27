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

Global
    MULESTEqn

Description
    Solves temperature scalar transport equation using MULES (Multidimensional 
    universal limiter for explicit solution) technology.    

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

{
    Field<scalar> maxT (T.size(), gMax(T));
    Field<scalar> minT (T.size(), gMin(T));
    
    MULES::explicitSolve
    (
        geometricOneField(),
        T,
        phi,
        phiTf,
        zeroField(),
        zeroField(),
        maxT,
        minT
    );

    if (implicitDiffusion)
    {
        solve(fvm::ddt(T) - fvc::ddt(T) - fvm::laplacian(Hif, T) == TSu);
    }
    else
    {
        solve(fvm::ddt(T) - fvc::ddt(T) - fvc::laplacian(Hif, T) == TSu);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //