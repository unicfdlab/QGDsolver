/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    selectbadCells

Description
    Picks up cells which can produce oscillations in QGD solution and puts
    them into the badCells cell set. Works only in serial mode.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "primitiveMeshTools.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "faceSet.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char* argv[])
{
    //    argList::addNote
    //    (
    //        "Create a cellSet for cells with their centres inside the defined "
    //        "surface.\n"
    //        "Surface must be closed and singly connected."
    //    );
    argList::noParallel();
    //argList::validArgs.append("surfaceFile");
    //argList::validArgs.append("cellSet");
#include "setRootCase.H"
#include "createTime.H"
#include "createPolyMesh.H"
    IOdictionary QGDCellQuality
    (
        IOobject
        (
            "QGDCellQuality",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    scalar faceCosine (readScalar(QGDCellQuality.lookup("faceCosine")));
    scalar maxAspectRatio(readScalar(QGDCellQuality.lookup("maxAspectRatio")));
    // Destination cellSet.
    cellSet badFaceAngle (mesh, "badFaceAngle", IOobject::NO_READ);
    const labelListList& pointFaces =
        mesh.pointFaces();
    const label nInternalFaces = mesh.nInternalFaces();
    forAll(pointFaces, ipoint)
    {
        forAll(pointFaces[ipoint], iface)
        {
            label ifaceId = pointFaces[ipoint][iface];
            label inei = -1;
            label iown = mesh.faceOwner()[ifaceId];

            if (ifaceId < nInternalFaces)
            {
                inei = mesh.faceNeighbour()[ifaceId];
            }

            forAll(pointFaces[ipoint], kface)
            {
                label kfaceId = pointFaces[ipoint][kface];

                if (ifaceId == kfaceId)
                {
                    continue;
                }

                label knei = -1;
                label kown = mesh.faceOwner()[kfaceId];

                if (kfaceId < nInternalFaces)
                {
                    knei = mesh.faceNeighbour()[kfaceId];
                }

                vector ni = mesh.faceAreas()[ifaceId] / mag(mesh.faceAreas()[ifaceId]);
                vector nk = mesh.faceAreas()[kfaceId] / mag(mesh.faceAreas()[kfaceId]);
                scalar dotnf = mag( ni & nk );

                if ((inei == knei) && (inei >= 0))
                {
                    if (dotnf >= faceCosine) // 5 degrees
                    {
                        badFaceAngle.insert(inei);
                    }
                }

                if ((iown == kown) || (iown == knei))
                {
                    if (dotnf >= faceCosine) // 5 degrees
                    {
                        badFaceAngle.insert(iown);
                    }
                }

                if (inei == kown)
                {
                    if (dotnf >= faceCosine) // 5 degrees
                    {
                        badFaceAngle.insert(inei);
                    }
                }
            }
        }
    }
    badFaceAngle.write();
    //
    // Select cells with high aspect ratio
    //
    scalarField openness(mesh.cellVolumes().size(), 0);
    scalarField aspectRatio(mesh.cellVolumes().size(), 1);
    primitiveMeshTools::cellClosedness
    (
        mesh,
        mesh.geometricD(),
        mesh.faceAreas(),
        mesh.cellVolumes(),
        openness,
        aspectRatio
    );
    cellSet highAspectRatio (mesh, "highAspectRatio", IOobject::NO_READ);
    forAll(aspectRatio, iCell)
    {
        if (aspectRatio[iCell] > maxAspectRatio)
        {
            highAspectRatio.insert(iCell);
        }
    }
    //    faceSet badFacesSet (mesh, "badFaces", IOobject::NO_READ);
    //
    //    mesh.checkFaceFlatness(true, 0.95, &badFacesSet);
    //
    //    badFacesSet.write();
    Info << "End\n" << endl;
    return 0;
}


// ************************************************************************* //
