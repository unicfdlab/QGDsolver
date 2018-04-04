#include "leastSquaresBase.H"

Foam::fvsc::
leastSquaresBase::leastSquaresBase(const fvMesh& mesh)
:
cMesh_(mesh)
{
    this->findNeighbours();
    this->calculateWeights();
}

Foam::fvsc::
leastSquaresBase::~leastSquaresBase()
{
}

//
//END-OF-FILE
//

