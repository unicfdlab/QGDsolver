

#ifndef QGDTHERMO_H
#define QGDTHERMO_H

#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "Switch.H"

namespace Foam
{
class fvMesh;
class dictionary;
namespace qgd
{
  class QGDCoeffs;
}

class QGDThermo
{

  //-
  const fvMesh& mesh_;

  //-
  const dictionary& dict_;

  //-
  autoPtr<Foam::qgd::QGDCoeffs> qgdCoeffsPtr_;

  //-
  Switch implicitDiffusion_;

private:

  //-
  QGDThermo();

  //-
  QGDThermo(const QGDThermo& );

protected:

  //-
  qgd::QGDCoeffs& qgdCoeffs();

  //-
  void correctQGD(volScalarField& mu, volScalarField& alphau);

public:

  //- Construct from mesh and dictionary
  QGDThermo(const fvMesh& mesh, const dictionary& dict);

  virtual ~QGDThermo();

  //-
  virtual const volScalarField& c() const = 0;

  //-
  virtual const volScalarField& p() const = 0;

  //-
  virtual tmp<volScalarField> rho() const = 0;

  //-
  virtual tmp<volScalarField> mu() const = 0;

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

  //-
  Switch implicitDiffusion() const;

  //-
  virtual bool read();

};
}

#endif

//
//END-OF-FILE
//
