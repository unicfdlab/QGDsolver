#include "QGDCoeffs.H"

#ifndef alphaHbyUQHD_H
#define alphaHbyUQHD_H

namespace Foam
{

namespace qgd
{

class alphaHbyUQHD : public QGDCoeffs
{

public:

    //-
    TypeName ("alphaHbyUQHD");
    
    //-
    alphaHbyUQHD
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );
    
    //-
    ~alphaHbyUQHD();
    
    //-
    void correct(const rhoQGDThermo&);

private:

    //-
    scalar uQGD_ ;

};


}//namespace qgd

}//namespace Foam

#endif

//
//END-OF-FILe
//


