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
#include "polyMesh.H"
#include "cellSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
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
    
    // Destination cellSet.
    cellSet badCellsSet (mesh, "badCells", IOobject::NO_READ);
    
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
                    
                    if (dotnf >= 0.996194698) // 5 degrees
                    {
                        badCellsSet.insert(inei);
                    }
                }
                
                if ((iown == kown) || (iown == knei))
                {

                    if (dotnf >= 0.996194698) // 5 degrees
                    {
                        badCellsSet.insert(iown);
                    }
                }
                
                if (inei == kown)
                {
                    if (dotnf >= 0.996194698) // 5 degrees
                    {
                        badCellsSet.insert(inei);
                    }
                }
                
            }
        }
    }
    
//    Info<< "Selected " << insideCells.size() << " of " << mesh.nCells()
//        << " cells" << nl << nl
//        << "Writing selected cells to cellSet " << insideCells.name()
//        << nl << nl
//        << "Use this cellSet e.g. with subsetMesh : " << nl << nl
//        << "    subsetMesh " << insideCells.name()
//        << nl << endl;
    
    badCellsSet.write();
    
    Info<< "End\n" << endl;
    
    return 0;
}


// ************************************************************************* //
