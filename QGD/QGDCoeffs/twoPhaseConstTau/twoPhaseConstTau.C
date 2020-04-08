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

\*---------------------------------------------------------------------------*/

#include "twoPhaseConstTau.H"
#include "QGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

namespace Foam
{
namespace qgd
{
    defineTypeNameAndDebug(twoPhaseConstTau,0);
    addToRunTimeSelectionTable
    (
        QGDCoeffs,
        twoPhaseConstTau,
        dictionary
    );
}
}

Foam::qgd::
twoPhaseConstTau::twoPhaseConstTau
(
    const IOobject& io,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    QGDCoeffs(io, mesh, dict)
{
    scalar ScQGD = 0.0, PrQGD = 1.0;
//
    ScQGD_.primitiveFieldRef() = ScQGD;
    PrQGD_.primitiveFieldRef() = PrQGD;
    muQGD_.primitiveFieldRef() = 0.0;
    alphauQGD_.primitiveFieldRef() = 0.0;

    ScQGD_.boundaryFieldRef() = ScQGD;
    PrQGD_.boundaryFieldRef() = PrQGD;
    muQGD_.boundaryFieldRef() = 0.0;
    alphauQGD_.boundaryFieldRef() = 0.0;
    
    dict.lookup("tau1") >> tau1_;
    dict.lookup("tau2") >> tau2_;
}

Foam::qgd::
twoPhaseConstTau::~twoPhaseConstTau()
{
}

void Foam::qgd::
twoPhaseConstTau::correct(const QGDThermo& qgdThermo)
{
    const volScalarField  nu     = qgdThermo.mu()/qgdThermo.rho();
                            //alphal
    if (mesh_.thisDb().found("alphal"))
    {
        const volScalarField& alpha1 = 
            mesh_.thisDb().lookupObject<volScalarField>("alphal");
        dimensionedScalar Tau1
        (
            "Tau1",
            dimensionSet(0,0,1,0,0),
            tau1_
        );
        dimensionedScalar Tau2
        (
            "Tau2",
            dimensionSet(0,0,1,0,0),
            tau2_
        );
        //this->tauQGD_ = alpha1*Tau1 + (1.0 - alpha1)*Tau2;
        this->tauQGD_ = pow(alpha1/Tau1 + (1.0 - alpha1)/Tau2,-1.0);
    }
    else
    {
        dimensionedScalar Tau
        (
            "Tau",
            dimensionSet(0,0,1,0,0),
            max(tau1_, tau2_)
        );
        this->tauQGD_ = Tau;
    }
    
    this->tauQGDf_= linearInterpolate(this->tauQGD_);
}

//
//END-OF-FILE
//
