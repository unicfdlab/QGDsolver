#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>

//- Find neighbour cells for each face (throught face points).
void Foam::extendedFaceStencil::findNeighbours()
{
    Pout << "Start findNeighbours()" << endl;
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

        //send local cell addressing to master process
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

        //send local point addressing to master process
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

        //redistribute cell addressing from master to slave processes
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

        //redistribute point addressing from master to slave processes
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
    
    Pout << "End findNeighbours()" << endl;
};

//
//END-OF-FILE
//

