#include "extendedFaceStencil.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <iostream>   // std::cout
#include <string>
#include <HashTable.H>


//- Calculate gradient of volume scalar function on the faces
//
// \param iF         Internal scalar field.
//                   Allowable values: constant reference to the volScalarField.
//
// \return           Gradient of iF (vector field) which was computed on the faces of mesh.
surfaceVectorField Foam::extendedFaceStencil::faceScalarGrad(const volScalarField& iF)
{

//    Pout << "faceScalarGrad" << endl;

    surfaceScalarField sF = linearInterpolate(iF);

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
    
    #warning "Add evaluation of gradient on external faces"

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
    
//    Pout << "Values from neighbouring patches transmitted to master" << endl;
    
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
    
//    Pout << "Information distributed across all other processors" << endl;

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
                        
                        gradForFaceForEachPatch[patchI][faceI] = gradF;
                    }
                }
            }
        }
    }
//    Pout << "Sending information to other processors" << endl;
//    Info << "existingPatches_ = " << existingPatches_ << endl;
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
//                Pout << "found patch pair (" << iProc << "," << jProc << ")" << endl;
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
//                                Pout << "Sending to " << patchPair.second() << endl;
                                gradIF.boundaryFieldRef()[patchI].operator = (gradForFaceForEachPatch[patchI]);
//                                Pout << "gradIF[] = " << gradIF.boundaryField()[patchI] << endl;
                                
                                OPstream toSlave (Pstream::scheduled, patchPair.second());
                                toSlave << gradForFaceForEachPatch[patchI];
//                                Pout << "Done" << endl;
                            }
                        }
                    }
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
//                                    Pout << "Recieving from " << patchPair.first() << endl;
                                    IPstream fromSlave (Pstream::scheduled, patchPair.first());
                                    
                                    fromSlave >> gradIF.boundaryFieldRef()[patchI];
//                                    Pout << "Done" << endl;
                                    
//                                    Pout << "gradIF[] = " << gradIF.boundaryField()[patchI] << endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
//    Pout << "... done" << endl;
    return gradIF;
};

//
//END-OF-FILE
//


