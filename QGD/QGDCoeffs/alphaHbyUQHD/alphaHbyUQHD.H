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
    alphaHbyUQHD

Description
    Class for one of possible ways of $\tau$ calculating.
    in this case $\tau$ parameter is calculating by

    $$ \tau = \alpha \cdot \frac{h}{U_0} $$

    Where:

    $\tau$ is regularisation parameter, $\alpha$ is parameter
    which depends by user in constant foulder, $h$ is size of the mesh
    and $U_0$ is characteristic velocity

SourceFiles
    alphaHbyUQHD.C

\*---------------------------------------------------------------------------*/




#include "QGDCoeffs.H"

#ifndef alphaHbyUQHD_H
#define alphaHbyUQHD_H

namespace Foam
{

namespace qgd
{

class alphaHbyUQHD : public QGDCoeffs
{

public:

    //-
    TypeName ("alphaHbyUQHD");

    //-
    alphaHbyUQHD
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~alphaHbyUQHD();

    //-
    void correct(const QGDThermo&);

private:

    //-
    scalar UQHD_ ;

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
