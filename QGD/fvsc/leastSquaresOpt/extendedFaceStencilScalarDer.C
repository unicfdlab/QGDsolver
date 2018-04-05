#include "leastSquaresStencilOpt.H"
#include "polyMesh.H"
#include "fvMesh.H"
#include "word.H"
#include "IOstream.H"
#include "Ostream.H"
#include <HashTable.H>


void Foam::fvsc::leastSquaresOpt::faceScalarDer(const Field<scalar>& iF,const Field<scalar>& sF,int com, surfaceScalarField& rField)
{
    forAll(sF, facei)
    {
        rField[facei] = 0.0;
        forAll(GdfAll_[facei], i)
        {
            rField[facei] += wf2All_[facei][i]*GdfAll_[facei][i].component(com)*(iF[neighbourCells_[facei][i]] - sF[facei]);
        }
    }

};

void Foam::fvsc::leastSquaresOpt::faceScalarDer(const tmp<Field<scalar>>& tiF,const tmp<Field<scalar>>& tsF, int com, tmp<surfaceScalarField>& trField )
{
};


void Foam::fvsc::leastSquaresOpt::faceScalarDer(const List3<scalar>& procVfValues,const surfaceScalarField& sF,int derComp, surfaceScalarField& rField)
{   
    scalar gf = 0.0;
    forAll(procPairs_, patchI)
    {
        label procPatchId = procPairs_[patchI];
        if (procPatchId > -1)
        {
            fvsPatchScalarField& pgradf = rField.boundaryFieldRef()[procPatchId];
            const List2<scalar> & pvf = procVfValues[patchI];
            const List2<scalar> & pwf2= procWf2_[patchI];
            const List2<vector> & pgdf= procGdf_[patchI];

            forAll(pgradf, iFace)
            {
                gf = 0.0;
                forAll(procGdf_[patchI][iFace], i)
                {
                    gf += pwf2[iFace][i]*pgdf[iFace][i].component(derComp)*
                       (pvf[iFace][i] - sF.boundaryField()[procPatchId][iFace]);
                }
                rField[iFace] = gf;
            }
        }
    }
};

void Foam::fvsc::leastSquaresOpt::faceScalarDer(const tmp<List3<scalar>>& tprocVfValues,const tmp<surfaceScalarField>&, int derComp, tmp<surfaceScalarField>& trField )
{
};

