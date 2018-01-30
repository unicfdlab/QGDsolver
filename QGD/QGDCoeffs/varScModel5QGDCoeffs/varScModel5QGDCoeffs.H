#include "QGDCoeffs.H"


#ifndef varScModel5QGDCoeffs_H
#define varScModel5QGDCoeffs_H

namespace Foam
{

class cellSet;

namespace qgd
{

class varScModel5QGDCoeffs : public QGDCoeffs
{

protected:

    scalar smoothCoeff_;
    scalar rC_;
    scalar minSc_;
    scalar maxSc_;
    scalar badQualitySc_;
    scalar qgdAspectRatioThreshold_;
    scalar constSc_;
    
    scalarField cqSc_;
    
    autoPtr<cellSet> constScCellSetPtr_;
    

public:

    //-
    TypeName ("varScModel5QGDCoeffs");
    
    //-
    varScModel5QGDCoeffs
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );
    
    //-
    ~varScModel5QGDCoeffs();
    
    //-
    void correct(const psiQGDThermo&);

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//


