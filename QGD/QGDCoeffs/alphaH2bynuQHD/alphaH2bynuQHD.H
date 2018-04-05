/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                QGDsolver   | Copyright (C) 2016-2017 ISP RAS (www.unicfd.ru)
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

Class
    alphaH2bynuQHD

Description
    Class for one of possible ways of tau calculating.

$$    \tau_QGD = \alpha \cdot \frac{h^2}{\nu}  $$

    where $\tau_QGD$ - regulatisation parameter, $\alpha$ - parameter
    which depends by user in constant foulder, $h$ - size of mesh,
    $\nu$ - kinetic viscosity.


SourceFiles
    alphaH2bynuQHD.C

\*---------------------------------------------------------------------------*/




#include "QGDCoeffs.H"

#ifndef alphaH2bynuQHD_H
#define alphaH2bynuQHD_H

namespace Foam
{

namespace qgd
{

class alphaH2bynuQHD : public QGDCoeffs
{

public:

    //-
    TypeName ("alphaH2bynuQHD");

    //-
    alphaH2bynuQHD
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~alphaH2bynuQHD();

    //-
    void correct(const QGDThermo&);
};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
