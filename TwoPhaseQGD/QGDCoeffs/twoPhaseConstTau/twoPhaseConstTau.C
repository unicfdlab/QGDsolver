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

\*---------------------------------------------------------------------------*/

#include "twoPhaseConstTau.H"
#include "QGDThermo.H"
#include "twoPhaseIcoQGDThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
}

Foam::qgd::
twoPhaseConstTau::~twoPhaseConstTau()
{
}

void Foam::qgd::
twoPhaseConstTau::correct(const QGDThermo& qgdThermo)
{
    const volScalarField  nu     = qgdThermo.mu()/qgdThermo.rho();
    if (isA<twoPhaseIcoQGDThermo>(qgdThermo))
    {
        const twoPhaseIcoQGDThermo& refThermo = 
            refCast<const twoPhaseIcoQGDThermo>(qgdThermo);
        
        const dimensionedScalar& Tau1 = refThermo.Tau1();
        const dimensionedScalar& Tau2 = refThermo.Tau2();
        const volScalarField& alpha1  = refThermo.alpha1();
        this->tauQGD_ = alpha1*Tau1 + (1.0 - alpha1)*Tau2;
        Info << "max/min tauQGD: " << max(tauQGD_).value() << "/" << min(tauQGD_).value() << endl;
    }
    else
    {
        FatalErrorInFunction
        << "twoPhaseConstTau::correct(const QGDThermo& qgdThermo)"
        << "MUST be only with two-phase solver interQHDFoam"
        << exit(FatalError);
    }

    this->tauQGDf_= linearInterpolate(this->tauQGD_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //