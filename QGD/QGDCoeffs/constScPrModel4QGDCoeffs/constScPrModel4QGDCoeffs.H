#include "QGDCoeffs.H"

#ifndef constScPrModel4QGDCoeffs_H
#define constScPrModel4QGDCoeffs_H

namespace Foam
{

namespace qgd
{

class constScPrModel4QGDCoeffs : public QGDCoeffs
{

protected:

    scalar smoothCoeff_;

public:

    //-
    TypeName ("constScPrModel4QGDCoeffs");
    
    //-
    constScPrModel4QGDCoeffs
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );
    
    //-
    ~constScPrModel4QGDCoeffs();
    
    //-
    void correct(const psiQGDThermo&);

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//


