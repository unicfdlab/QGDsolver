
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "regIOobject.H"
#include "volFields.H"

#ifndef QGDCoeffs_H
#define QGDCoeffs_H

namespace Foam
{

//forward declaration of rhoQGDThermo
class rhoQGDThermo;

namespace qgd
{

class QGDCoeffs : public regIOobject, public refCount
{

protected:
    
    //-
    const fvMesh& mesh_;
    
    //-
    const Time& runTime_;
    
    //-
    volScalarField muQGD_;
    
    //-
    volScalarField alphauQGD_;

    //-
    volScalarField beta_;

    //-
    volScalarField hQGD_;
    
    //-
    volScalarField aQGD_;
    
    //-
    volScalarField tauQGD_;
    
    //-
    volScalarField PrQGD_;
    
    //-
    volScalarField ScQGD_;

protected:
    
    //-
    virtual tmp<volScalarField> qgdLength(const fvMesh& mesh);

public:
 
    //-
    TypeName("QGDCoeffs");
    
    //-
    declareRunTimeSelectionTable
    (
        autoPtr,
        QGDCoeffs,
        dictionary,
        (
            const IOobject& io,
            const fvMesh& mesh,
            const dictionary& dict
        ),
        (io, mesh, dict)
    );
    
    //- Construct from components
    QGDCoeffs
    (
        const IOobject& io,
        const fvMesh& mesh,
        const dictionary& dict
    );
    
    //- Return a reference to the selected fvscStencil model
    static autoPtr<QGDCoeffs> New
    (
        const word& qgdCoeffsType,
        const fvMesh& mesh,
        const dictionary& dict
    );
    
    //-
    virtual ~QGDCoeffs();
    
    //-
    virtual bool writeData(Ostream&) const
    {
        return true;
    }
    
    //-
    virtual void correct(const rhoQGDThermo&);
    
    //-
    const volScalarField& hQGD() const;
    
    //-
    const volScalarField& tauQGD() const;
    
    //-
    const volScalarField& muQGD() const;
    
    //-
    const volScalarField& alphauQGD() const;

    //-
    const volScalarField& beta() const;

};

}

}

#endif

//
//END-OF-FILE
//


