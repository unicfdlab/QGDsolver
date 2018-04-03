#include "QGDCoeffs.H"

#ifndef alphaH2bynuQHD_H
#define alphaH2bynuQHD_H

namespace Foam
{

namespace qgd
{

class alphaH2bynuQHD : public QGDCoeffs
{

public:

    //-
    TypeName ("alphaH2bynuQHD");
    
    //-
    alphaH2bynuQHD
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );
    
    //-
    ~alphaH2bynuQHD();
    
    //-
    void correct(const rhoQGDThermo&);
};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//


