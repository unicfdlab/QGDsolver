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

Class
    Foam::RSWETau

Description
    Class for calculating tau for Regularized Shallow Water Equations model.

    $$    \tau_{RSWE} = \alpha \cdot \frac{\Delta x}{\sqrt{g h}}  $$

    where $\tau_{RSWE}$ - regularization parameter, $\alpha$ - tuning parameter
    which depends by user in constant foulder, $\Delta x$ - size of mesh,
    $g$ - gravity factor, $h$ -  thickness of the water layer measured from the 
    bottom.


SourceFiles
    RSWETau.C

\*---------------------------------------------------------------------------*/


#include "QGDCoeffs.H"

#ifndef RSWETau_H
#define RSWETau_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace qgd
{

class RSWETau : public QGDCoeffs
{

public:

    //-
    TypeName ("RSWETau");

    //-
    RSWETau
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~RSWETau();

    //-
    void correct(const QGDThermo&);
};


}//namespace qgd

}//namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// ************************************************************************* //