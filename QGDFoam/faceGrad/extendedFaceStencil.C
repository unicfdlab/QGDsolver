#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <iostream>   // std::cout
#include <string>
#include <HashTable.H>

// constructors
Foam::extendedFaceStencil::extendedFaceStencil(const IOobject& io, const bool isTime)
:
    regIOobject(io, isTime),
    mesh_(refCast<const fvMesh>(io.db()))
{
    findNeighbours();
    calculateWeights();
}


Foam::extendedFaceStencil::extendedFaceStencil(const regIOobject& rio)
:
    regIOobject(rio),
    mesh_(refCast<const fvMesh>(rio.db()))
{
    findNeighbours();
    calculateWeights();
}

Foam::extendedFaceStencil::extendedFaceStencil(const regIOobject& rio, bool registerCopy)
:
    regIOobject(rio, registerCopy),
    mesh_(refCast<const fvMesh>(rio.db()))
{
    findNeighbours();
    calculateWeights();
}

Foam::extendedFaceStencil::~extendedFaceStencil()
{
}

void Foam::extendedFaceStencil::rename(const word& newName)
{
};

bool Foam::extendedFaceStencil::readData(Istream&)
{
  return false;
};

bool Foam::extendedFaceStencil::read()
{
  return false;
};

bool Foam::extendedFaceStencil::modified()  const
{
  return false;
};

bool Foam::extendedFaceStencil::readIfModified()
{
  return false;
};

bool Foam::extendedFaceStencil::writeData(Ostream&) const
{
  return false;
};

bool Foam::extendedFaceStencil::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber,
    IOstream::compressionType
)  const
{
  return false;
};

bool Foam::extendedFaceStencil::write()  const
{
  return false;
};

//- Find neighbour cells for each face (throught face points).
void Foam::extendedFaceStencil::findNeighbours()
{
    Pout << "findNeighbours()" << endl;
    // List of faces
    const faceList& faces = mesh_.faces();

    // Or for all faces?
    labelListList neighbourCellsForFace(mesh_.nInternalFaces());

    forAll(faces, facei)
    {
        if(mesh_.isInternalFace(facei))
        {
            labelList neighbourCells;
            labelList pointsFacei = faces[facei];

            forAll(pointsFacei, k)
            {
                label pointi = pointsFacei[k];
                // cells should be added to list if it wasn't contain in the list yet
                labelList neighbourCellsForPointI = mesh_.pointCells()[pointi];

                forAll(neighbourCellsForPointI, j)
                {
                    label celli = neighbourCellsForPointI[j];

                    bool contained = false;

                    forAll(neighbourCells, f)
                    {
                        label ncell=neighbourCells[f];
                        if(ncell==celli)
                        {
                            contained = true;
                        }
                    }

                    if(!contained)
                    {
                        neighbourCells.append(celli);
                    }
                }
            }
            neighbourCellsForFace[facei].append(neighbourCells);
        }
    }

    neighbourCells_ = neighbourCellsForFace;

    if(Pstream::parRun())
    {
        List<List<label> > processCellToGlobalAddr_;
        //List<List<label> > globalCellToProcessAddr_;
        List<List<label> > processPointToGlobalAddr_;
        //List<List< List<label> > > globalPointToProcessAddr_;

        const fvMesh& localMesh = mesh_;
        const Time& localTime   = mesh_.time();

        //
        // Step 1:
        // Initialize globalCellToProcessAddr_
        //
        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                const word globalConstant = localTime.rootPath() + "/" + localTime.globalCaseName() + "/constant";
                const fileName gRootPath = localTime.rootPath();
                const fileName gCaseName = localTime.globalCaseName();

                // suspend MPI
                Pstream::parRun() = false;
                //label comm        = Pstream::worldComm;
                List<label> oldProcIDs_       = Pstream::procID(Pstream::worldComm);//comm);
                List<label> newProcIDs_       = List<label> (1, 0);
                Pstream::procID(Pstream::worldComm) = newProcIDs_;//comm) = newProcIDs_;
                // end suspend MPI

                autoPtr<Time> globalTimePtr_
                (
                    new Time
                    (
                        gRootPath,
                        gCaseName
                    )
                );

                Time& globalTime = globalTimePtr_();

                autoPtr<fvMesh> globalMeshPtr_
                (
                    new fvMesh
                    (
                        IOobject
                        (
                            fvMesh::defaultRegion,
                            globalTime.timeName(),
                            globalTime,
                            IOobject::MUST_READ
                        )
                    )
                );

                // resume MPI
                //label comm            = Pstream::worldComm;
                Pstream::procID(Pstream::worldComm) = oldProcIDs_;//comm) = oldProcIDs_;
                Pstream::parRun()     = true;
                // end resume MPI

                globalCellToProcessAddr_.resize(globalMeshPtr_().V().size());
                globalPointToProcessAddr_.resize(globalMeshPtr_().points().size());
                List<vector> globalC_ = globalMeshPtr_().C().internalField();
                List<scalar> globalV_ = globalMeshPtr_().V();

                //send number of cells in global mesh
                for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
                {
                    OPstream toSlave(Pstream::scheduled, jSlave);
                    toSlave << globalCellToProcessAddr_.size();
                }
            }
            else
            {
                //receive number of cells in global mesh from master process
                IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
                label nGlobalCells = 0;
                fromMaster >> nGlobalCells;
                globalCellToProcessAddr_.resize(nGlobalCells);
            }
        }

        // for globalPointToProcessAddr_
        if (Pstream::parRun())
        {
            if (Pstream::master())
            {
                //send number of points in global mesh
                for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
                {
                    OPstream toSlave(Pstream::scheduled, jSlave);
                    toSlave << globalPointToProcessAddr_.size();
                }
            }
            else
            {
                //receive number of points in global mesh from master process
                IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
                label nGlobalPoints = 0;
                fromMaster >> nGlobalPoints;
                globalPointToProcessAddr_.resize(nGlobalPoints);
            }
        }

        //
        // Step 2:
        // Filling globalCellToProcessAddr_ by numbers
        //
        processCellToGlobalAddr_.resize
        (
            Pstream::nProcs()
        );

        processPointToGlobalAddr_.resize
        (
            Pstream::nProcs()
        );

        //read local cell addressing
        labelIOList localCellProcAddr
        (
           IOobject
            (
                "cellProcAddressing",
                localMesh.facesInstance(),
                localMesh.meshSubDir,
                localMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        labelIOList localPointProcAddr
        (
           IOobject
            (
                "pointProcAddressing",
                localMesh.facesInstance(),
                localMesh.meshSubDir,
                localMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        processCellToGlobalAddr_[Pstream::myProcNo()] = localCellProcAddr;
        processPointToGlobalAddr_[Pstream::myProcNo()] = localPointProcAddr;

        //send local cell addressing master process
        if (Pstream::master())
        {
            for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
            {
                IPstream fromSlave(Pstream::scheduled, jSlave);
                label nSlaveCells = 0;
                fromSlave >> nSlaveCells;
                processCellToGlobalAddr_[jSlave].resize(nSlaveCells);
                labelList& slaveCellProcAddr = processCellToGlobalAddr_[jSlave];

                forAll(slaveCellProcAddr, iCell)
                {
                    fromSlave >> slaveCellProcAddr[iCell];
                }
            }
        }
        else
        {
            OPstream toMaster (Pstream::scheduled, Pstream::masterNo());
            toMaster << localCellProcAddr.size();

            forAll(localCellProcAddr, iCell)
            {
                toMaster << localCellProcAddr[iCell];
            }
        }

        //send local point addressing master process
        if (Pstream::master())
        {
            for (label jSlave=Pstream::firstSlave(); jSlave<=Pstream::lastSlave(); jSlave++)
            {
                IPstream fromSlave(Pstream::scheduled, jSlave);
                label nSlavePoints = 0;
                fromSlave >> nSlavePoints;
                processPointToGlobalAddr_[jSlave].resize(nSlavePoints);
                labelList& slavePointProcAddr = processPointToGlobalAddr_[jSlave];

                forAll(slavePointProcAddr, iPoint)
                {
                    fromSlave >> slavePointProcAddr[iPoint];
                }
            }
        }
        else
        {
            OPstream toMaster (Pstream::scheduled, Pstream::masterNo());
            toMaster << localPointProcAddr.size();

            forAll(localPointProcAddr, iPoint)
            {
                toMaster << localPointProcAddr[iPoint];
            }
        }

        //redistribute cell addressing to slave processes
        if (Pstream::master())
        {
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                OPstream toSlave (Pstream::scheduled, jSlave);

                forAll(processCellToGlobalAddr_, iProcess)
                {
                    const labelList& thisProcessAddr = processCellToGlobalAddr_[iProcess];
                    const label nCells = thisProcessAddr.size();
                    toSlave << nCells;

                    forAll(thisProcessAddr, jCell)
                    {
                        toSlave << thisProcessAddr[jCell];
                    }
                }
            }
        }
        else
        {
            IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());

            forAll(processCellToGlobalAddr_, iProcess)
            {
                labelList& thisProcessAddr = processCellToGlobalAddr_[iProcess];
                label nCells = 0;
                fromMaster >> nCells;
                thisProcessAddr.resize(nCells);

                forAll(thisProcessAddr, jCell)
                {
                    fromMaster >> thisProcessAddr[jCell];
                }
            }
        }

        //redistribute point addressing to slave processes
        if (Pstream::master())
        {
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                OPstream toSlave (Pstream::scheduled, jSlave);

                forAll(processPointToGlobalAddr_, iProcess)
                {
                    const labelList& thisProcessAddr = processPointToGlobalAddr_[iProcess];
                    const label nPoints = thisProcessAddr.size();
                    toSlave << nPoints;

                    forAll(thisProcessAddr, jPoint)
                    {
                        toSlave << thisProcessAddr[jPoint];
                    }
                }
            }
        }
        else
        {
            IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());

            forAll(processPointToGlobalAddr_, iProcess)
            {
                labelList& thisProcessAddr = processPointToGlobalAddr_[iProcess];
                label nPoints = 0;
                fromMaster >> nPoints;
                thisProcessAddr.resize(nPoints);

                forAll(thisProcessAddr, jPoint)
                {
                    fromMaster >> thisProcessAddr[jPoint];
                }
            }
        }

        // for cells
        forAll(processCellToGlobalAddr_, jProc)
        {
            const labelList& jProcessAddr = processCellToGlobalAddr_[jProc];

            forAll(jProcessAddr, iCell)
            {
                label iGlobalCell = jProcessAddr[iCell];
                labelList procAndLocalNo(2,0);
                procAndLocalNo[0] = jProc;
                procAndLocalNo[1] = iCell;
                globalCellToProcessAddr_[iGlobalCell].append(procAndLocalNo);
            }
        }
        // for points
        forAll(processPointToGlobalAddr_, jProc)
        {
            const labelList& jProcessAddr = processPointToGlobalAddr_[jProc];

            forAll(jProcessAddr, iPoint)
            {
                label iGlobalPoint = jProcessAddr[iPoint];
                if(globalPointToProcessAddr_[iGlobalPoint].size()==0)
                {
                    labelList procNo(1, jProc);
                    labelList localNo(1, iPoint);
                    globalPointToProcessAddr_[iGlobalPoint].resize(2);
                    globalPointToProcessAddr_[iGlobalPoint][0] = procNo;
                    globalPointToProcessAddr_[iGlobalPoint][1] = localNo;
                }
                else
                {
                    labelList procNo = globalPointToProcessAddr_[iGlobalPoint][0];
                    procNo.append(jProc);
                    labelList localNo = globalPointToProcessAddr_[iGlobalPoint][1];
                    localNo.append(iPoint);
                    globalPointToProcessAddr_[iGlobalPoint][0] = procNo;
                    globalPointToProcessAddr_[iGlobalPoint][1] = localNo;
                }
            }
        }

        // for patches
        List< List< labelList > > neibCellsForPatch(mesh_.boundary().size());          //free
        List< labelList >         pointsForEachPatch(mesh_.boundary().size());
        List< List< labelList > > neibCellsForPoints(mesh_.boundary().size());
        List< List< labelList > > neibGlobalCellsForPoints(mesh_.boundary().size());

        forAll(mesh_.boundary(), patchIndex)
        {
            if (isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
            {
                const fvPatch& p = mesh_.boundary()[patchIndex];
                const processorFvPatch& pp = refCast<const processorFvPatch>(p);

                const label startFaceNo = pp.start();
                const label sizeFaces = pp.size();

                labelList facesForPatch(sizeFaces, 0);
                forAll(facesForPatch, faceI)
                {
                    facesForPatch[faceI] = startFaceNo + faceI;
                }

                labelList pointsForPatch;

                forAll(facesForPatch, faceI)
                {
                    labelList pointsForFace = mesh_.faces()[facesForPatch[faceI]];

                    forAll(pointsForFace, pointI)
                    {
                        bool ex = false;
                        forAll(pointsForPatch, pointJ)
                        {
                            if(pointsForFace[pointI]==pointsForPatch[pointJ])
                            {
                                ex = true;
                            }
                        }
                        if(!ex)
                        {
                            pointsForPatch.append(pointsForFace[pointI]);
                        }
                    }
                }

                pointsForEachPatch[patchIndex] = pointsForPatch;

                neibCellsForPoints[patchIndex].resize(pointsForPatch.size());
                neibGlobalCellsForPoints[patchIndex].resize(pointsForPatch.size());

                forAll(neibCellsForPoints[patchIndex], pointI)
                {
                    neibCellsForPoints[patchIndex][pointI] = mesh_.pointCells()[pointsForEachPatch[patchIndex][pointI]];
                    forAll(neibCellsForPoints[patchIndex][pointI], cellI)
                    {
                        label localCellNo = neibCellsForPoints[patchIndex][pointI][cellI];
                        neibGlobalCellsForPoints[patchIndex][pointI].append(processCellToGlobalAddr_[pp.myProcNo()][localCellNo]);
                    }
                }
            }

        }

        labelList boundPointsGlobal;
        labelList boundPointsLocal;
        List < labelList > neibCellsForBoundPoints;

        forAll(globalPointToProcessAddr_, pointIndex)
        {
            if(globalPointToProcessAddr_[pointIndex][0].size() > 1)
            {

                labelList procList = globalPointToProcessAddr_[pointIndex][0];
                labelList pointList = globalPointToProcessAddr_[pointIndex][1];
                
                while(procList.size() > 1)
                {
                    label indexOfMainProc = findMin(procList);
                    label mainProc = procList[indexOfMainProc];
                    labelList newProcList;
                    labelList newPointList;

                    forAll(procList, procIndex)
                    {
                        if(procList[procIndex] != mainProc)
                        {
                            newProcList.append(procList[procIndex]);
                            newPointList.append(pointList[procIndex]);
                        }
                    }

                    if(Pstream::myProcNo() == mainProc)
                    {
                        boundPointsGlobal.append(pointIndex);
                        boundPointsLocal.append(pointList[indexOfMainProc]);

                        labelList forAdd;

                        forAll(procList, procIndex)
                        {
                            if(procList[procIndex]!=mainProc)
                            {
                                label procI = procList[procIndex];
                                label lengthOfCells = 0;

                                IPstream fromSlave(Pstream::scheduled, procI);

                                fromSlave >> lengthOfCells;

                                labelList addToNeibCells(lengthOfCells);
                                forAll(addToNeibCells, cellI)
                                {
                                    fromSlave >> addToNeibCells[cellI];
                                }

                                forAdd.append(addToNeibCells);
                            }
                        }

                        neibCellsForBoundPoints.append(forAdd);
                    }
                    else
                    {
                        forAll(procList, index)
                        {
                            if(Pstream::myProcNo() == procList[index])
                            {
                                label patchNo = 0;
                                label pointOnPatchNo = 0;

                                forAll(mesh_.boundary(), patchIndex)
                                {
                                    if(isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
                                    {
                                        forAll(pointsForEachPatch[patchIndex], pointI)
                                        {
                                            if(pointsForEachPatch[patchIndex][pointI] == pointList[index])
                                            {
                                                patchNo = patchIndex;
                                                pointOnPatchNo = pointI;
                                            }
                                        }
                                    }
                                }

                                OPstream toSlave(Pstream::scheduled, mainProc);
                                toSlave << neibGlobalCellsForPoints[patchNo][pointOnPatchNo].size();

                                forAll(neibGlobalCellsForPoints[patchNo][pointOnPatchNo], cellI)
                                {
                                    toSlave << neibGlobalCellsForPoints[patchNo][pointOnPatchNo][cellI];
                                }
                            }
                        }
                    }


                    procList.resize(procList.size() - 1);
                    pointList.resize(pointList.size() - 1);
                    procList = newProcList;
                    pointList = newPointList;
                }
            }
        }

        forAll(globalPointToProcessAddr_, pointIndex)
        {

            if(globalPointToProcessAddr_[pointIndex][0].size() > 1)
            {
                label mainProc = globalPointToProcessAddr_[pointIndex][0][findMin(globalPointToProcessAddr_[pointIndex][0])];
                label notMainProc = globalPointToProcessAddr_[pointIndex][0][findMax(globalPointToProcessAddr_[pointIndex][0])];
                label pointOnMain = globalPointToProcessAddr_[pointIndex][1][findMin(globalPointToProcessAddr_[pointIndex][0])];

                if(Pstream::myProcNo() == mainProc)
                {
                        forAll(boundPointsGlobal, pointI)
                        {
                            if(boundPointsGlobal[pointI]==pointIndex)
                            {
                                bool iSend = false;
                                forAll(mesh_.boundary(), patchI)
                                {
                                    if(isType<processorFvPatch>(mesh_.boundary()[patchI])&&(!iSend))
                                    {
                                        forAll(pointsForEachPatch[patchI], pointJ)//global???
                                        {
                                            if(pointsForEachPatch[patchI][pointJ]==pointOnMain)
                                            {
                                                forAll(neibGlobalCellsForPoints[patchI][pointJ], cellI)
                                                {
                                                    neibCellsForBoundPoints[pointI].append(neibGlobalCellsForPoints[patchI][pointJ][cellI]);
                                                }
                                                iSend = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
            }

            if(globalPointToProcessAddr_[pointIndex][0].size() > 2)
            {
                label mainProc = globalPointToProcessAddr_[pointIndex][0][findMin(globalPointToProcessAddr_[pointIndex][0])];
                label notMainProc = globalPointToProcessAddr_[pointIndex][0][findMax(globalPointToProcessAddr_[pointIndex][0])];
                label pointOnMain = globalPointToProcessAddr_[pointIndex][1][findMin(globalPointToProcessAddr_[pointIndex][0])];

                if(Pstream::myProcNo() == mainProc)
                {
                    forAll(globalPointToProcessAddr_[pointIndex][0], procI)
                    {

                        if(globalPointToProcessAddr_[pointIndex][0][procI]!=mainProc && globalPointToProcessAddr_[pointIndex][0][procI]!=notMainProc)
                        {
                            forAll(boundPointsGlobal, pointI)
                            {
                                if(boundPointsGlobal[pointI]==pointIndex)
                                {
                                    OPstream toSlave(Pstream::scheduled, globalPointToProcessAddr_[pointIndex][0][procI]);
                                    toSlave << neibCellsForBoundPoints[pointI].size();

                                    forAll(neibCellsForBoundPoints[pointI], cellI)
                                    {
                                        toSlave << neibCellsForBoundPoints[pointI][cellI];
                                    }
                                }
                            }
                        }
                    }
                }
                else
                {
                    forAll(globalPointToProcessAddr_[pointIndex][0], procI)
                    {
                        if( Pstream::myProcNo()==globalPointToProcessAddr_[pointIndex][0][procI] && Pstream::myProcNo()!=notMainProc)
                        {

                            label lengthOfCells;
                            IPstream fromSlave(Pstream::scheduled, mainProc);
                            fromSlave >> lengthOfCells;

                            labelList neibCells(lengthOfCells);
                            forAll(neibCells, cellI)
                            {
                                fromSlave >> neibCells[cellI];
                            }

                            forAll(boundPointsGlobal, pointI)
                            {
                                if(boundPointsGlobal[pointI]==pointIndex)
                                {
                                    neibCellsForBoundPoints[pointI].resize(lengthOfCells);
                                    neibCellsForBoundPoints[pointI] = neibCells;
                                }
                            }

                        }
                    }
                }
            }
        }

        boundPointsGlobal_.resize(Pstream::nProcs());
        boundPointsLocal_.resize(Pstream::nProcs());
        neibCellsForBoundPoints_.resize(Pstream::nProcs());

        // collect information from all processors
        if (Pstream::master())
        {
            boundPointsGlobal_[Pstream::myProcNo()] = boundPointsGlobal;
            boundPointsLocal_[Pstream::myProcNo()] = boundPointsLocal;
            neibCellsForBoundPoints_[Pstream::myProcNo()] = neibCellsForBoundPoints;

            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                scalar lenBoundPoints = 0;
                IPstream fromSlave (Pstream::scheduled, jSlave);
                fromSlave >> lenBoundPoints;

                labelList addBoundPointsGlobal(lenBoundPoints);
                forAll(addBoundPointsGlobal, pointI)
                {
                    fromSlave >> addBoundPointsGlobal[pointI];
                }
                boundPointsGlobal_[jSlave] = addBoundPointsGlobal;

                labelList addBoundPointsLocal(lenBoundPoints);
                forAll(addBoundPointsLocal, pointI)
                {
                    fromSlave >> addBoundPointsLocal[pointI];
                }
                boundPointsLocal_[jSlave] = addBoundPointsLocal;

                List < labelList > addNeibCellsForBoundPoints(lenBoundPoints);
                scalar lenCells = 0;
                forAll(addNeibCellsForBoundPoints, pointI)
                {
                    fromSlave >> lenCells;
                    addNeibCellsForBoundPoints[pointI].resize(lenCells);
                    forAll(addNeibCellsForBoundPoints[pointI], cellI)
                    {
                        fromSlave >> addNeibCellsForBoundPoints[pointI][cellI];
                    }
                }
                neibCellsForBoundPoints_[jSlave] = addNeibCellsForBoundPoints;
            }
        }
        else
        {
            OPstream toMaster (Pstream::scheduled, Pstream::masterNo());

            toMaster << boundPointsGlobal.size();

            forAll(boundPointsGlobal, pointI)
            {
                toMaster << boundPointsGlobal[pointI];
            }

            forAll(boundPointsLocal, pointI)
            {
                toMaster << boundPointsLocal[pointI];
            }

            forAll(neibCellsForBoundPoints, pointI)
            {
                toMaster << neibCellsForBoundPoints[pointI].size();
                forAll(neibCellsForBoundPoints[pointI], cellI)
                {
                    toMaster << neibCellsForBoundPoints[pointI][cellI];
                }
            }

        }

        // send information for all processors
        if (Pstream::master())
        {
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                forAll (boundPointsGlobal_, procI)
                {
                    OPstream toSlave (Pstream::scheduled, jSlave);
                    toSlave << boundPointsGlobal_[procI].size();

                    forAll(boundPointsGlobal_[procI], pointI)
                    {
                        toSlave << boundPointsGlobal_[procI][pointI];
                    }

                    forAll(boundPointsLocal_[procI], pointI)
                    {
                        toSlave << boundPointsLocal_[procI][pointI];
                    }

                    forAll(neibCellsForBoundPoints_[procI], pointI)
                    {
                        toSlave << neibCellsForBoundPoints_[procI][pointI].size();
                        forAll(neibCellsForBoundPoints_[procI][pointI], cellI)
                        {
                             toSlave << neibCellsForBoundPoints_[procI][pointI][cellI];
                        }
                    }
                }
            }
        }
        else
        {
            forAll (boundPointsGlobal_, procI)
            {
                IPstream fromMaster (Pstream::scheduled, Pstream::masterNo());
                scalar lenBoundPoints = 0;
                fromMaster >> lenBoundPoints;

                labelList addBoundPointsGlobal(lenBoundPoints);
                forAll(addBoundPointsGlobal, pointI)
                {
                    fromMaster >> addBoundPointsGlobal[pointI];
                }
                boundPointsGlobal_[procI] = addBoundPointsGlobal;

                labelList addBoundPointsLocal(lenBoundPoints);
                forAll(addBoundPointsLocal, pointI)
                {
                    fromMaster >> addBoundPointsLocal[pointI];
                }
                boundPointsLocal_[procI] = addBoundPointsLocal;

                List <labelList> addNeibCellsForBoundPoints(lenBoundPoints);
                forAll(addNeibCellsForBoundPoints, pointI)
                {
                    scalar lenNeibCells = 0;
                    fromMaster >> lenNeibCells;
                    addNeibCellsForBoundPoints[pointI].resize(lenNeibCells);
                    forAll(addNeibCellsForBoundPoints[pointI], cellI)
                    {
                         fromMaster >> addNeibCellsForBoundPoints[pointI][cellI];
                    }
                }
                neibCellsForBoundPoints_[procI] = addNeibCellsForBoundPoints;
            }
        }

//        HashTable <bool, Pair<label>, Hash< Pair<label> > > existingPatches_;
        List < Pair<label> > listPatches;

        for(int iProc=0; iProc<Pstream::nProcs(); iProc++)
        {
            for(int jProc=iProc+1; jProc<Pstream::nProcs(); jProc++)
            {
                if(Pstream::myProcNo() == iProc)
                {
                    forAll(mesh_.boundary(), patchIndex)
                    {
                        if(isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
                        {
                            const fvPatch& p = mesh_.boundary()[patchIndex];
                            const processorFvPatch& pp = refCast<const processorFvPatch>(p);

                            if(pp.neighbProcNo() == jProc)
                            {
                                Pair<label> procPatchPair;
                                procPatchPair.first()=iProc;
                                procPatchPair.second()=jProc;

                                existingPatches_.insert(procPatchPair, true);
                                listPatches.append(procPatchPair);
                            }
                        }
                    }

                    forAll(neibCellsForBoundPoints, bpIndex)
                    {
                        if(neibCellsForBoundPoints[bpIndex].size() == 0)
                        {
                            label myPointGlobal = boundPointsGlobal[bpIndex];
                            labelList procList = globalPointToProcessAddr_[myPointGlobal][0];
                            label mainProc = procList[findMin(procList)];
                            label myPointLocal = globalPointToProcessAddr_[myPointGlobal][1][findMin(procList)];
                        }
                    }
                }
            }
        }

        //Pout << "existingPatches_: " << existingPatches_ << endl;
        // нужно собрать информацию о всех патчах, а потом раздать ее
        
        // собираем
        if(Pstream::master())
        {
            Pair<label> procPatchPair;
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                label pairNumber;
                IPstream fromSlave (Pstream::scheduled, jSlave);
                fromSlave >> pairNumber;
                
                for(int i=1; i<=pairNumber; i++)
                {
                    fromSlave >> procPatchPair;
                    existingPatches_.insert(procPatchPair, true);
                    listPatches.append(procPatchPair);
                }
            }
        }
        else
        {
            //send
            OPstream toMaster (Pstream::scheduled, Pstream::masterNo());
            toMaster << listPatches.size();
            
            forAll(listPatches, pairI)
            {
                toMaster << listPatches[pairI];
            }
        }
        
        // раздаем
        if(Pstream::master())
        {
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                OPstream toSlave (Pstream::scheduled, jSlave);
                toSlave << listPatches.size();

                forAll(listPatches, pairI)
                {
                    toSlave << listPatches[pairI];
                }
            }
        }
        else
        {
            existingPatches_.clear();
            listPatches.clear();

            Pair<label> procPatchPair;
            label pairNumber;

            IPstream fromMaster (Pstream::scheduled, Pstream::masterNo());
            fromMaster >> pairNumber;

            for(int i=1; i<=pairNumber; i++)
            {
                fromMaster >> procPatchPair;
                existingPatches_.insert(procPatchPair, true);
                listPatches.append(procPatchPair);
            }
        }

    }
};

//- Compute weights for least squares scheme for gradient calculation.
void Foam::extendedFaceStencil::calculateWeights()
{

    Info << "calculateWeights" << endl;
    Pout << "calculateWeights()" << endl;

    const faceList& faces = mesh_.faces();

    List< List<vector> > GdfAll(faces.size());
    List< scalarList > wf2All(faces.size());

    //for internal faces
    forAll(faces, facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            List<vector> df(neighbourCells_[facei].size());
            scalarList wf2(neighbourCells_[facei].size());
            symmTensor G(0);

            forAll(neighbourCells_[facei], i)
            {
                df[i] = mesh_.cellCentres()[neighbourCells_[facei][i]] - mesh_.faceCentres()[facei];
                wf2[i] = 1/magSqr(df[i]);
                symmTensor addToG(0);
                addToG = sqr(df[i]);
                addToG = addToG * wf2[i];
                G += addToG;
            }

            symmTensor G0(0);

            if(mesh_.nGeometricD()==1)
            {
                symmTensor G01(1, 0, 0, 1, 0, 1);
                symmTensor G02(sqr((Vector<label>::one + mesh_.geometricD())/2));
                G0 = G01 - G02;
            }
            else
            {
                if(mesh_.nGeometricD()==2)
                {
                    G0 = sqr((Vector<label>::one - mesh_.geometricD())/2);
                }
            };

            G = G + G0;
            //Info << "nSolution: " << mesh_.nGeometricD() << endl;
            //Info << "G0 = " << G0 << endl;
            G = inv(G);
            G = G - G0;

            forAll(df, i)
            {
                df[i] = G&df[i];
            }

            GdfAll[facei] = df;
            wf2All[facei] = wf2;
        }
    }
    
    Info << "End for not parallel" << endl;
    
    GdfAll_ = GdfAll;
    wf2All_ = wf2All;
    
    // если это процессор патч
    if(Pstream::parRun())
    {
        centersNeibCellsForBoundPoints_.resize(Pstream::nProcs());
    
        forAll(neibCellsForBoundPoints_, procI)
        {
            forAll(neibCellsForBoundPoints_[procI], pointI)
            {
                centersNeibCellsForBoundPoints_[procI].resize(neibCellsForBoundPoints_[procI].size());
                labelList cells = neibCellsForBoundPoints_[procI][pointI];
                forAll(cells, cellI)
                {
                    if(Pstream::master())
                    {
                        // wait
                        centersNeibCellsForBoundPoints_[procI][pointI].resize(neibCellsForBoundPoints_[procI][pointI].size());

                        label procJ = globalCellToProcessAddr_[cells[cellI]][0];
                        label cellJ = globalCellToProcessAddr_[cells[cellI]][1];
                        point cellCentre;
                        if(procJ != Pstream::myProcNo())
                        {
                            IPstream fromSlave (Pstream::scheduled, procJ);
                            fromSlave >> cellCentre;
                        }
                        else
                        {
                            cellCentre = mesh_.cellCentres()[cellJ];
                        }
                        centersNeibCellsForBoundPoints_[procI][pointI][cellI] = cellCentre;
                    }
                    else
                    {
                        //send
                        if(Pstream::myProcNo()==globalCellToProcessAddr_[cells[cellI]][0])
                        {
                            label cellJ = globalCellToProcessAddr_[cells[cellI]][1];
                            OPstream toMaster (Pstream::scheduled, Pstream::masterNo());
                            toMaster << mesh_.cellCentres()[cellJ];
                        }
                    }
                }
            }
        }

        // send information for all processors
        if (Pstream::master())
        {
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                forAll (centersNeibCellsForBoundPoints_, procI)
                {
                    OPstream toSlave (Pstream::scheduled, jSlave);
                    toSlave << centersNeibCellsForBoundPoints_[procI].size();

                    forAll(centersNeibCellsForBoundPoints_[procI], pointI)
                    {
                        toSlave << centersNeibCellsForBoundPoints_[procI][pointI].size();
                        forAll(centersNeibCellsForBoundPoints_[procI][pointI], cellI)
                        {
                            toSlave << centersNeibCellsForBoundPoints_[procI][pointI][cellI];
                        }
                    }
                }
            }
        }
        else
        {
            forAll (centersNeibCellsForBoundPoints_, procI)
            {
                IPstream fromMaster (Pstream::scheduled, Pstream::masterNo());
                scalar lenBoundPoints = 0;
                fromMaster >> lenBoundPoints;

                List < List <point> > addCentersNeibCellsForBoundPoints(lenBoundPoints);
                forAll(addCentersNeibCellsForBoundPoints, pointI)
                {
                    scalar lenNeibCells = 0;
                    fromMaster >> lenNeibCells;
                    addCentersNeibCellsForBoundPoints[pointI].resize(lenNeibCells);
                    forAll(addCentersNeibCellsForBoundPoints[pointI], cellI)
                    {
                        fromMaster >> addCentersNeibCellsForBoundPoints[pointI][cellI];
                    }
                }
                centersNeibCellsForBoundPoints_[procI] = addCentersNeibCellsForBoundPoints;
            }
        }

        // раздать centersNeibCellsForBoundPoints_ - есть!
        // посчитать веса на "своих" патчах - есть!
        // раздать веса согласно хэш-таблице - может не нужно?

        List < List < labelList >  > neibCellsForFaceForEachPatch(mesh_.boundary().size());
        //List < List < List < point > > > centersNeibCellsForFaceForEachPatch(mesh_.boundary().size());
        centersNeibCellsForFaceForEachPatch_.resize(mesh_.boundary().size());

        forAll(mesh_.boundary(), patchIndex)
        {
            if(isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
            {
                const fvPatch& p = mesh_.boundary()[patchIndex];
                const processorFvPatch& pp = refCast<const processorFvPatch>(p);
                List < labelList  > neibCellsForFaceFromPatch(pp.size());
                List < List < point > > centersNeibCellsForFaceFromPatch(pp.size());

                if(pp.myProcNo() < pp.neighbProcNo())
                {
                    labelList facesFromPatch(pp.size());
                    List < labelList > pointForFaceFromPatch(pp.size());

                    forAll(facesFromPatch, faceI)
                    {
                        facesFromPatch[faceI] = pp.start() + faceI;
                        pointForFaceFromPatch[faceI] = mesh_.faces()[facesFromPatch[faceI]];

                        forAll(pointForFaceFromPatch[faceI], pointI)
                        {
                            label currentPoint = pointForFaceFromPatch[faceI][pointI];

                            forAll(boundPointsLocal_[Pstream::myProcNo()], pointJ)
                            {
                                label boundPoint = boundPointsLocal_[Pstream::myProcNo()][pointJ];
                                if( currentPoint == boundPoint )
                                {
                                    forAll(neibCellsForBoundPoints_[Pstream::myProcNo()][pointJ], cellJ)
                                    {
                                        bool exist = false;
                                        forAll(neibCellsForFaceFromPatch[faceI], cellI)
                                        {
                                            if(neibCellsForBoundPoints_[Pstream::myProcNo()][pointJ][cellJ]==neibCellsForFaceFromPatch[faceI][cellI])
                                            {
                                                exist = true;
                                            }
                                        }
                                        if(!exist)
                                        {
                                            neibCellsForFaceFromPatch[faceI].append(neibCellsForBoundPoints_[Pstream::myProcNo()][pointJ][cellJ]);
                                            centersNeibCellsForFaceFromPatch[faceI].append(centersNeibCellsForBoundPoints_[Pstream::myProcNo()][pointJ][cellJ]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                neibCellsForFaceForEachPatch[patchIndex].append(neibCellsForFaceFromPatch);
                centersNeibCellsForFaceForEachPatch_[patchIndex].append(centersNeibCellsForFaceFromPatch);
            }
        }

        GdfForEachPatch_.resize(mesh_.boundary().size());
        wf2ForEachPatch_.resize(mesh_.boundary().size());

        forAll(mesh_.boundary(), patchI)
        {
            if(neibCellsForFaceForEachPatch[patchI].size() > 0)
            {
                if(neibCellsForFaceForEachPatch[patchI][0].size() > 0)
                {
                    const fvPatch& p = mesh_.boundary()[patchI];
                    const processorFvPatch& pp = refCast<const processorFvPatch>(p);

                    labelList facesFromPatch(pp.size());
                    forAll(facesFromPatch, faceI)
                    {
                        facesFromPatch[faceI] = pp.start() + faceI;
                    }

                    GdfForEachPatch_[patchI].resize(neibCellsForFaceForEachPatch[patchI].size());
                    wf2ForEachPatch_[patchI].resize(neibCellsForFaceForEachPatch[patchI].size());

                    forAll(mesh_.boundary()[patchI], faceI)
                    {
                        GdfForEachPatch_[patchI][faceI].resize(neibCellsForFaceForEachPatch[patchI][faceI].size());
                        wf2ForEachPatch_[patchI][faceI].resize(neibCellsForFaceForEachPatch[patchI][faceI].size());

                        vector zeroVector(0, 0, 0);
                        List<vector> df(neibCellsForFaceForEachPatch[patchI][faceI].size(), zeroVector);
                        GdfForEachPatch_[patchI][faceI] = df;

                        scalar zero(0);
                        scalarList wf2(neibCellsForFaceForEachPatch[patchI][faceI].size(), zero);
                        wf2ForEachPatch_[patchI][faceI] = wf2;

                        symmTensor G(0);

                        forAll(df, cellI)
                        {
                            df[cellI] = centersNeibCellsForFaceForEachPatch_[patchI][faceI][cellI] - mesh_.faceCentres()[facesFromPatch[faceI]];
                            wf2[cellI] = 1/magSqr(df[cellI]);
                            symmTensor addToG(0);
                            addToG = sqr(df[cellI]);
                            addToG = addToG * wf2[cellI];
                            G += addToG;
                        }

                        symmTensor G0(0);

                        if(mesh_.nSolutionD()==1)
                        {
                            symmTensor G01(1, 0, 0, 1, 0, 1);
                            symmTensor G02(sqr((Vector<label>::one + mesh_.geometricD())/2));
                            G0 = G01 - G02;
                        }
                        else
                        {
                            if(mesh_.nSolutionD()==2)
                            {
                                G0 = sqr((Vector<label>::one - mesh_.geometricD())/2);
                            }
                        };

                        G = G + G0;

                        G = inv(G);
                        G = G - G0;

                        forAll(df, cellI)
                        {
                            df[cellI] = G & df[cellI];
                        }

                        GdfForEachPatch_[patchI][faceI] = df;
                        wf2ForEachPatch_[patchI][faceI] = wf2;
                    }
                }
            }
        }

    }

};

//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
surfaceVectorField Foam::extendedFaceStencil::faceScalarGrad(const volScalarField& iF)
{

//    Pout << "faceScalarGrad" << endl;

    surfaceScalarField sF = fvc::interpolate(iF);

    surfaceVectorField gradIF("gradIF", 0*fvc::snGrad(iF)  * mesh_.Sf() / mesh_.magSf());

    // List of faces
    const faceList& faces = mesh_.faces();

    forAll(faces, facei)
    {
        vector gradF = vector::zero;

        if (mesh_.isInternalFace(facei))
        {
            forAll(GdfAll_[facei], i)
            {
                gradF = gradF + wf2All_[facei][i]*GdfAll_[facei][i]*(iF[neighbourCells_[facei][i]] - sF[facei]);
            }

            gradIF[facei] = gradF;
        }
    }

    if(!Pstream::parRun())
    {
        return gradIF;
    }

    // если это процессор патч
    if(Pstream::parRun())
    {
        valueInNeibCellsForBoundPoints_.resize(Pstream::nProcs());

        forAll(neibCellsForBoundPoints_, procI)
        {
            forAll(neibCellsForBoundPoints_[procI], pointI)
            {
                valueInNeibCellsForBoundPoints_[procI].resize(neibCellsForBoundPoints_[procI].size());
                labelList cells = neibCellsForBoundPoints_[procI][pointI];
                forAll(cells, cellI)
                {
                    if(Pstream::master())
                    {
                        // wait
                        valueInNeibCellsForBoundPoints_[procI][pointI].resize(neibCellsForBoundPoints_[procI][pointI].size());

                        label procJ = globalCellToProcessAddr_[cells[cellI]][0];
                        label cellJ = globalCellToProcessAddr_[cells[cellI]][1];
                        scalar cellValue;
                        if(procJ != Pstream::myProcNo())
                        {
                            IPstream fromSlave (Pstream::scheduled, procJ);
                            fromSlave >> cellValue;
                        }
                        else
                        {
                            cellValue = iF[cellJ];//mesh_.cellCentres()[cellJ];
                        }
                        valueInNeibCellsForBoundPoints_[procI][pointI][cellI] = cellValue;
                    }
                    else
                    {
                        //send
                        if(Pstream::myProcNo()==globalCellToProcessAddr_[cells[cellI]][0])
                        {
                            label cellJ = globalCellToProcessAddr_[cells[cellI]][1];
                            OPstream toMaster (Pstream::scheduled, Pstream::masterNo());
                            toMaster << iF[cellJ];//mesh_.cellCentres()[cellJ];
                        }
                    }
                }
            }
        }
        
//        Pout << "values2 first: " << valueInNeibCellsForBoundPoints_ << endl;

        // send information for all processors
        if (Pstream::master())
        {
            for (label jSlave = Pstream::firstSlave(); jSlave <= Pstream::lastSlave(); jSlave++)
            {
                forAll (valueInNeibCellsForBoundPoints_, procI)
                {
                    OPstream toSlave (Pstream::scheduled, jSlave);
                    toSlave << valueInNeibCellsForBoundPoints_[procI].size();

                    forAll(valueInNeibCellsForBoundPoints_[procI], pointI)
                    {
                        toSlave << valueInNeibCellsForBoundPoints_[procI][pointI].size();
                        forAll(valueInNeibCellsForBoundPoints_[procI][pointI], cellI)
                        {
                            toSlave << valueInNeibCellsForBoundPoints_[procI][pointI][cellI];
                        }
                    }
                }
            }
        }
        else
        {
            forAll (valueInNeibCellsForBoundPoints_, procI)
            {
                IPstream fromMaster (Pstream::scheduled, Pstream::masterNo());
                scalar lenBoundPoints = 0;
                fromMaster >> lenBoundPoints;

                List < List <scalar> > addValueInNeibCellsForBoundPoints(lenBoundPoints);
                forAll(addValueInNeibCellsForBoundPoints, pointI)
                {
                    scalar lenNeibCells = 0;
                    fromMaster >> lenNeibCells;
                    addValueInNeibCellsForBoundPoints[pointI].resize(lenNeibCells);
                    forAll(addValueInNeibCellsForBoundPoints[pointI], cellI)
                    {
                        fromMaster >> addValueInNeibCellsForBoundPoints[pointI][cellI];
                    }
                }
                valueInNeibCellsForBoundPoints_[procI] = addValueInNeibCellsForBoundPoints;
            }
        }

        // раздать valueInNeibCellsForBoundPoints_ - есть!
        // посчитать веса на "своих" патчах - есть!
        // раздать веса согласно хэш-таблице - может не нужно?

        List < List < labelList >  > neibCellsForFaceForEachPatch(mesh_.boundary().size());
        //List < List < List < point > > > valueInNeibCellsForFaceForEachPatch(mesh_.boundary().size());
        valueInNeibCellsForFaceForEachPatch_.resize(mesh_.boundary().size());

        forAll(mesh_.boundary(), patchIndex)
        {
            if(isType<processorFvPatch>(mesh_.boundary()[patchIndex]))
            {
                const fvPatch& p = mesh_.boundary()[patchIndex];
                const processorFvPatch& pp = refCast<const processorFvPatch>(p);
                List < labelList  > neibCellsForFaceFromPatch(pp.size());
                List < List < scalar > > valueInNeibCellsForFaceFromPatch(pp.size());
                valueInNeibCellsForFaceForEachPatch_[patchIndex].resize(pp.size());

                if(pp.myProcNo() < pp.neighbProcNo())
                {
                    labelList facesFromPatch(pp.size());
                    List < labelList > pointForFaceFromPatch(pp.size());

                    forAll(facesFromPatch, faceI)
                    {
                        facesFromPatch[faceI] = pp.start() + faceI;
                        pointForFaceFromPatch[faceI] = mesh_.faces()[facesFromPatch[faceI]];

                        forAll(pointForFaceFromPatch[faceI], pointI)
                        {
                            label currentPoint = pointForFaceFromPatch[faceI][pointI];

                            forAll(boundPointsLocal_[Pstream::myProcNo()], pointJ)
                            {
                                label boundPoint = boundPointsLocal_[Pstream::myProcNo()][pointJ];
                                if( currentPoint == boundPoint )
                                {
                                    forAll(neibCellsForBoundPoints_[Pstream::myProcNo()][pointJ], cellJ)
                                    {
                                        bool exist = false;
                                        forAll(neibCellsForFaceFromPatch[faceI], cellK)
                                        {
                                            if(neibCellsForBoundPoints_[Pstream::myProcNo()][pointJ][cellJ]==neibCellsForFaceFromPatch[faceI][cellK])
                                            {
                                                exist = true;
                                            }
                                        }
                                        if(!exist)
                                        {
                                            neibCellsForFaceFromPatch[faceI].append(neibCellsForBoundPoints_[Pstream::myProcNo()][pointJ][cellJ]);
                                            valueInNeibCellsForFaceFromPatch[faceI].append(valueInNeibCellsForBoundPoints_[Pstream::myProcNo()][pointJ][cellJ]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                neibCellsForFaceForEachPatch[patchIndex].append(neibCellsForFaceFromPatch);
                //valueInNeibCellsForFaceForEachPatch_[patchIndex].append(valueInNeibCellsForFaceFromPatch);
                valueInNeibCellsForFaceForEachPatch_[patchIndex] = valueInNeibCellsForFaceFromPatch;
            }
        }

    //}

//        Pout << "GdfForEachPatch_ " <<  GdfForEachPatch_ << endl;
//        Pout << "wf2ForEachPatch_ " <<  wf2ForEachPatch_ << endl;
//        Pout << "values: " << valueInNeibCellsForFaceForEachPatch_ << endl;
//        Pout << "values2: " << valueInNeibCellsForBoundPoints_ << endl;

        List < List < vector > > gradForFaceForEachPatch(mesh_.boundary().size());
        List < labelList > facesFromEachPatch(mesh_.boundary().size());

        forAll(mesh_.boundary(), patchI)
        {
            if(GdfForEachPatch_[patchI].size() > 0)
            {
                if(GdfForEachPatch_[patchI][0].size() > 0)
                {
                    if (isType<processorFvPatch>(mesh_.boundary()[patchI]))
                    {
                        const fvPatch& p = mesh_.boundary()[patchI];
                        const processorFvPatch& pp = refCast<const processorFvPatch>(p);

                        labelList facesFromPatch(pp.size());
                        forAll(facesFromPatch, faceI)
                        {
                            facesFromPatch[faceI] = pp.start() + faceI;
                        }
                        
                        facesFromEachPatch[patchI] = facesFromPatch;

                        gradForFaceForEachPatch[patchI].resize(pp.size());

                        forAll(GdfForEachPatch_[patchI], faceI)
                        {
                            label myFace = facesFromPatch[faceI];
                            vector gradF = vector::zero;

                            forAll(GdfForEachPatch_[patchI][faceI], cellI)
                            {
                                vector Gdf = GdfForEachPatch_[patchI][faceI][cellI];
                                scalar wf2 = wf2ForEachPatch_[patchI][faceI][cellI];
                                scalar sF1 = sF.boundaryField()[patchI][faceI];
                                gradF = gradF + wf2 * Gdf * (valueInNeibCellsForFaceForEachPatch_[patchI][faceI][cellI] - sF1);
                            }
                            
/*                            Pout << "for gradF" << endl;
                            Pout << "Gdf " << GdfForEachPatch_[patchI][faceI] << endl;
                            Pout << "wf2 " << wf2ForEachPatch_[patchI][faceI] << endl;
                            Pout << "value " << valueInNeibCellsForFaceForEachPatch_[patchI][faceI] << endl;
                            Pout << "sF[myFace] " << sF[myFace] << endl;
                            Pout << "myFace " << myFace << endl;*/

                            gradForFaceForEachPatch[patchI][faceI] = gradF;
//                            Pout << "gradF " << gradF << endl;
                        }
                    }
                }
            }
        }
        
        //Pout << "gradForFaceForEachPatch: " << gradForFaceForEachPatch << endl;

        // нужно только раздать
        for(int iProc=0; iProc<Pstream::nProcs(); iProc++)
        {
            for(int jProc=iProc+1; jProc<Pstream::nProcs(); jProc++)
            {
                Pair<label> patchPair;
                patchPair.first() = iProc;
                patchPair.second() = jProc;
                
                if(existingPatches_.found(patchPair))
                {
                
                    if(Pstream::myProcNo()==patchPair.first())
                    {
                        forAll(mesh_.boundary(), patchI)
                        {
                            if (isType<processorFvPatch>(mesh_.boundary()[patchI]))
                            {
                                const fvPatch& p = mesh_.boundary()[patchI];
                                const processorFvPatch& pp = refCast<const processorFvPatch>(p);
                                
                                if(pp.neighbProcNo()==patchPair.second())
                                {
                                    OPstream toSlave (Pstream::scheduled, patchPair.second());

                                    gradIF.boundaryField()[patchI] == gradForFaceForEachPatch[patchI];

                                    toSlave << gradForFaceForEachPatch[patchI];
                                }
                                
//                                Pout << "gradIF: (from first) " << gradForFaceForEachPatch[patchI] << endl;
                            }
                        
                        }
                        //Pout << "gradIF: (from first) " << gradIF << endl;
                        //бегаю по процессорным патчам, ищу, у кого сосед со вторым номером из пары --- тому и нужно раздать
                        //себе присваиваю, раздаю соседу
                    }
                    else
                    {
                        if(Pstream::myProcNo()==patchPair.second())
                        {
                            forAll(mesh_.boundary(), patchI)
                            {
                                if (isType<processorFvPatch>(mesh_.boundary()[patchI]))
                                {
                                    const fvPatch& p = mesh_.boundary()[patchI];
                                    const processorFvPatch& pp = refCast<const processorFvPatch>(p);

                                    if(pp.neighbProcNo()==patchPair.first())
                                    {
                                        IPstream fromSlave (Pstream::scheduled, patchPair.first());

                                        fromSlave >> gradIF.boundaryFieldRef()[patchI];
                                    }
                                    
                                    //Pout << "gradIF: (from second) " << gradIF.boundaryField()[patchI] << endl;
                                }
                            }
                            //бегаю по процессорным патчам, ищу, у кого сосед с первым номером из пары --- от того и нужно принимать
                            //принимаю и присваиваю себе
                            //Pout << "gradIF: (from second) " << gradIF << endl;
                        }
                    }
                }
            }
        }



    }

    return gradIF;
};

//- Calculate gradient of volume vector field on the faces.
//
// \param iVF      Internal vector field.
//                 Allowable values: constant reference to the volVectorField.
//
// \return         Gradient of iVF (tensor field) which was computed on the faces of mesh.
surfaceTensorField Foam::extendedFaceStencil::faceVectorGrad(const volVectorField& iVF)
{

//    Pout << "faceVectorGrad" << endl;

    surfaceVectorField gradComp0col = faceScalarGrad(iVF.component(0));
    surfaceVectorField gradComp1col = faceScalarGrad(iVF.component(1));
    surfaceVectorField gradComp2col = faceScalarGrad(iVF.component(2));

    surfaceTensorField gradIVF("divITF", 0*fvc::snGrad(iVF) * mesh_.Sf() / mesh_.magSf());

    forAll(mesh_.faces(), facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            gradIVF[facei].component(0) = gradComp0col[facei].component(0);
            gradIVF[facei].component(1) = gradComp1col[facei].component(0);
            gradIVF[facei].component(2) = gradComp2col[facei].component(0);

            gradIVF[facei].component(3) = gradComp0col[facei].component(1);
            gradIVF[facei].component(4) = gradComp1col[facei].component(1);
            gradIVF[facei].component(5) = gradComp2col[facei].component(1);

            gradIVF[facei].component(6) = gradComp0col[facei].component(2);
            gradIVF[facei].component(7) = gradComp1col[facei].component(2);
            gradIVF[facei].component(8) = gradComp2col[facei].component(2);
        }
    }

    return gradIVF;
};



//- Calculate divergence of volume vector field on the faces.
//
// \param iVF        Internal vector field.
//                   Allowable values: constant reference to the volVectorField.
//
// \return           Divergence of iVF (scalar field) which was computed on the faces of mesh.
surfaceScalarField Foam::extendedFaceStencil::faceVectorDiv(const volVectorField& iVF)
{
//    Pout << "faceVectorDiv" << endl;

    surfaceVectorField gradComp0 = faceScalarGrad(iVF.component(0));
    surfaceVectorField gradComp1 = faceScalarGrad(iVF.component(1));
    surfaceVectorField gradComp2 = faceScalarGrad(iVF.component(2));

    surfaceScalarField divIVF("divITF", 0*fvc::snGrad(iVF) & mesh_.Sf() / mesh_.magSf());

    forAll(mesh_.faces(), facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            divIVF[facei] = gradComp0[facei].component(0) + gradComp1[facei].component(1) + gradComp2[facei].component(2);
        }
    }

    return divIVF;
};

//- Calculate divergence of volume tensor field on the faces.
//
// \param iTF        Internal tensor field.
//                   Allowable values: constant reference to the volTensorField.
//
// \return           Divergence of iTF (vector field) which was computed on the faces of mesh.
surfaceVectorField Foam::extendedFaceStencil::faceTensorDiv(const volTensorField& iTF)
{
//    Pout << "faceTensorDiv" << endl;

    surfaceVectorField gradComp0 = faceScalarGrad(iTF.component(0));
    surfaceVectorField gradComp1 = faceScalarGrad(iTF.component(1));
    surfaceVectorField gradComp2 = faceScalarGrad(iTF.component(2));

    surfaceVectorField gradComp3 = faceScalarGrad(iTF.component(3));
    surfaceVectorField gradComp4 = faceScalarGrad(iTF.component(4));
    surfaceVectorField gradComp5 = faceScalarGrad(iTF.component(5));

    surfaceVectorField gradComp6 = faceScalarGrad(iTF.component(6));
    surfaceVectorField gradComp7 = faceScalarGrad(iTF.component(7));
    surfaceVectorField gradComp8 = faceScalarGrad(iTF.component(8));

    surfaceScalarField divComp0 = gradComp0.component(0) + gradComp3.component(1) + gradComp6.component(2);
    surfaceScalarField divComp1 = gradComp1.component(0) + gradComp4.component(1) + gradComp7.component(2);
    surfaceScalarField divComp2 = gradComp2.component(0) + gradComp5.component(1) + gradComp8.component(2);

    surfaceVectorField divITF("divITF", 0*fvc::snGrad(iTF.component(0)) * mesh_.Sf() / mesh_.magSf());

    forAll(mesh_.faces(), facei)
    {
        if (mesh_.isInternalFace(facei))
        {
            divITF[facei].component(0) = divComp0[facei];
            divITF[facei].component(1) = divComp1[facei];
            divITF[facei].component(2) = divComp2[facei];
        }
    }
    return divITF;
}

