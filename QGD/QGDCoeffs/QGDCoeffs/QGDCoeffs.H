
#include "fvMesh.H"
#include "runTimeSelectionTables.H"
#include "regIOobject.H"
#include "volFields.H"
#include "surfaceFields.H"

#ifndef QGDCoeffs_H
#define QGDCoeffs_H

namespace Foam
{

//forward declaration of psiQGDThermo
class psiQGDThermo;

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
    surfaceScalarField hQGDf_;

    //-
    volScalarField hQGD_;
    
    //-
    volScalarField aQGD_;
    
    //-
    volScalarField tauQGD_;
    
    //-
    surfaceScalarField tauQGDf_;
    
    //-
    volScalarField PrQGD_;
    
    //-
    volScalarField ScQGD_;

protected:
    
    //-
    virtual void updateQGDLength(const fvMesh&);

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
    virtual void correct(const psiQGDThermo&);
    
    //-
    const volScalarField& hQGD() const;
    
    //-
    const volScalarField& tauQGD() const;

    //-
    const surfaceScalarField& hQGDf() const;
    
    //-
    const surfaceScalarField& tauQGDf() const;

    //-
    const volScalarField& muQGD() const;
    
    //-
    const volScalarField& alphauQGD() const;

};

}

}

#endif

//
//END-OF-FILE
//


