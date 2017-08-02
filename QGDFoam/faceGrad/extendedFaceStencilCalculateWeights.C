#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>

//- Compute weights for least squares scheme for gradient calculation.
void Foam::extendedFaceStencil::calculateWeights()
{
    Pout << "Start calculateWeights()" << endl;

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
    
    Pout << "End calculateWeights()" << endl;
};

