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

Class
    Foam::qgd::constScPrModel2

Description
    Class describing QGD model coefficients using dynamic values for
    Sc^{QGD}, Pr^{QGD} and specified \alpha field:

    Sc^{QGD} = min(\h^QGD abs(\nabla \rho)/\rho, 1)

    \tau^{QGD} = \alpha h / c + \mu / (p Sc)

    \mu^{QGD} = (\alpha h / c) p Sc^{QGD}

    \lambda^{QGD} = \mu^{QGD} Cp / Pr^{QGD}

SourceFiles
    constScPrModel2.C

\*---------------------------------------------------------------------------*/
#include "QGDCoeffs.H"

#ifndef varScModel5_H
#define varScModel5_H

namespace Foam
{

class cellSet;

namespace qgd
{

class varScModel5 : public QGDCoeffs
{

protected:

    scalar smoothCoeff_;
    scalar rC_;
    scalar minSc_;
    scalar maxSc_;
    scalar badQualitySc_;
    scalar qgdAspectRatioThreshold_;
    scalar constSc_;

    scalarField cqSc_;

    autoPtr<cellSet> constScCellSetPtr_;


public:

    //-
    TypeName ("varScModel5");

    //-
    varScModel5
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~varScModel5();

    //-
    void correct(const QGDThermo&);

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
