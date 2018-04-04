#ifndef leastSquaresBase_H
#define leastSquaresBase_H

#include "volFields.H"
#include "surfaceFields.H"
#include "surfaceMesh.H"

#include "fvCFD.H"

namespace Foam
{

namespace fvsc
{

class leastSquaresBase
{

private:

        const fvMesh& cMesh_;

protected:

        template <class T>
        using List2 = List<List<T> >;

        template <class T>
        using List3 = List2<List<T> >;
    
        template <class T>
        using labelHashTable = HashTable<T,label>;
    
        template <class T>
        using Triple = FixedList<T,3>;

protected:

        //Data for internal field
        
        // List of neighbour cells for each face
        labelListList neighbourCells_;
        
        // Weights
        List2<vector> GdfAll_;
        List2<scalar> wf2All_;
        
        // List of faces with degenerate stencil
        DynamicList<label> internalDegFaces_;
        
        //Data for patches
        /* For faces connecting 2 processors  */
        List<label> procPairs_; // if >= 0, connects processor through patch
        List<label> neigProcs_; // No of neighbouring processor
        labelHashTable<label> idProcPatchPairs_; //hash for processors, connected through the patch
        label nProcPatches_;
        
        /* Least square method weights at processor patches */
        List3<vector> procGdf_;
        List3<scalar> procWf2_;
        
        /* Addressing information */
        List3<label> myProcPatchCells_;
        List2<label> ownEnd_;
        List2<label> neiStart_;
        List2<label> neiEnd_;
        List2<label> corStart_;
        List2<label> corEnd_;
        
        /* Global information about domain decomposition */
        //labelHashTable <List<label> > pointProcs_; //list of processors containing point (global id)
        labelHashTable <List<label> > pointCells_; //list of cells containing cell (global id)
        labelHashTable <label> cellProc_;
        
        /* Addressing information for corner neighbouring processors */
        List2<Triple<label> > corAddr_;
        List2<label> corCellIds_;
        labelHashTable<label> corProcIds_;
        
        List<DynamicList<label> > procDegFaces_;

protected:

        void findNeighbours();

        void calculateWeights();

public:

    leastSquaresBase(const fvMesh& mesh);
    
    ~leastSquaresBase();

};

}

}

#endif

//
//END-OF-FILE
//


