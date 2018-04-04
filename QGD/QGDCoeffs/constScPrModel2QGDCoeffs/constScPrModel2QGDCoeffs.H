#include "QGDCoeffs.H"

#ifndef constScPrModel2QGDCoeffs_H
#define constScPrModel2QGDCoeffs_H

namespace Foam
{

namespace qgd
{

class constScPrModel2QGDCoeffs : public QGDCoeffs
{

public:

    //-
    TypeName ("constScPrModel2QGDCoeffs");

    //-
    constScPrModel2QGDCoeffs
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );

    //-
    ~constScPrModel2QGDCoeffs();

    //-
    void correct(const QGDThermo&);

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//
