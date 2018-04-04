#include "QGDCoeffs.H"

#ifndef constScPrModel1QGDCoeffs_H
#define constScPrModel1QGDCoeffs_H

namespace Foam
{

namespace qgd
{

class constScPrModel1QGDCoeffs : public QGDCoeffs
{

public:

    //-
    TypeName ("constScPrModel1QGDCoeffs");

    //-
    constScPrModel1QGDCoeffs
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~constScPrModel1QGDCoeffs();

    //-
    void correct(const QGDThermo&);

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
