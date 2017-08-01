#ifndef extendedFaceStencil_H
#define extendedFaceStencil_H

#include "regIOobject.H"
#include "labelList.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "vector.H"
#include "List.H"

#include "processorFvPatch.H"

namespace Foam
{


class fvMesh;

/*---------------------------------------------------------------------------*\
                    Class extendedFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class extendedFaceStencil
:
    public regIOobject
{

public:
        // List of neighbour cells for each face
        labelListList neighbourCells_;
        
        List<List<label> > globalCellToProcessAddr_;
        List<List< List<label> > > globalPointToProcessAddr_;
        
        List < labelList > boundPointsGlobal_;
        List < labelList > boundPointsLocal_;
        List < List < labelList > > neibCellsForBoundPoints_;
        List < List < List< point > > > centersNeibCellsForBoundPoints_;
        List < List < List < point > > > centersNeibCellsForFaceForEachPatch_;
        List < List < List < scalar > > > valueInNeibCellsForBoundPoints_;
        List < List < List < scalar > > > valueInNeibCellsForFaceForEachPatch_;
        
        HashTable <bool, Pair<label>, Hash< Pair<label> > > existingPatches_;
        
//        List < List < vector > > gradForFaceForEachPatch;

protected:

private:
        // Mesh
        const fvMesh& mesh_;
        List< List<vector> > GdfAll_;
        List< scalarList > wf2All_;
        
        List < List< List<vector> > > GdfForEachPatch_;
        List < List< scalarList > > wf2ForEachPatch_;

public:

// Constructors

        //- Construct from IOobject. Optional flag for if IOobject is the
        //  top level regIOobject.
        extendedFaceStencil(const IOobject&, const bool isTime = false);

        //- Construct as copy
        extendedFaceStencil(const regIOobject&);
        
        //- Construct as copy, transferring registry registration to copy
        //  if registerCopy is true
        extendedFaceStencil(const regIOobject&, bool registerCopy);

        //- Destructor
        ~extendedFaceStencil();

        void rename(const word& newName);

        bool readData(Istream&);

        bool read();

        bool modified() const;

        bool readIfModified();

        bool writeData(Ostream&) const;

        bool writeObject
        (
            IOstream::streamFormat,
            IOstream::versionNumber,
            IOstream::compressionType
        ) const;

        bool write() const;

        void findNeighbours();

        void calculateWeights();

        surfaceVectorField faceScalarGrad(const volScalarField& iF);

        surfaceTensorField faceVectorGrad(const volVectorField& iVF);

        surfaceScalarField faceVectorDiv(const volVectorField& iVF);

        surfaceVectorField faceTensorDiv(const volTensorField& iTF);

};

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif