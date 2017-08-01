/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    QGDFoam

Description
    Quasi-Gasdynamic solver.
    Now only for orthogonal mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
//#include "Switch.H"
#include "vector.H"
#include "List.H"

labelListList neighbourCellsForEachFace(const fvMesh& mesh)
{
// Find all neighbour cells for each internal face
    // List of faces
    const faceList& faces=mesh.faces();

    // Or for all faces?
    labelListList neighbourCellsForFace(mesh.nInternalFaces());

    forAll(faces, facei)
    {
        if(mesh.isInternalFace(facei))
        {
          labelList neighbourCells;
          labelList pointsFacei = faces[facei];

          //forAll(pointsFacei, pointi)
          for (int k=0; k<pointsFacei.size(); k++)
          {
              label pointi = pointsFacei[k];
              // cells should be added to list if it wasn't contain in the list yet
              labelList neighbourCellsForPointI = mesh.pointCells()[pointi];
              //forAll(neighbourCellsForPointI, celli)
              for(int j=0; j<neighbourCellsForPointI.size(); j++)
              {
                  label celli = neighbourCellsForPointI[j];

                  bool contained = false;

                  //forAll(neighbourCells, ncell)
                  for(int f=0; f<neighbourCells.size(); f++)
                  {
                      label ncell=neighbourCells[f];
                      if(ncell==celli)
                        {
                            contained = true;
                        }
                  }

                  if(!contained)
                  {
                      neighbourCells.append(celli);
                  }
                }
            }
            neighbourCellsForFace[facei].append(neighbourCells);
        }
    }
// End Find all neighbour cells for each face
    return neighbourCellsForFace;
}

surfaceVectorField faceScalarGrad(const fvMesh& mesh, const volScalarField& iF)
{
    surfaceScalarField sF=fvc::interpolate(iF);

    surfaceVectorField gradIF("gradIF", fvc::snGrad(iF)  * mesh.Sf() / mesh.magSf());

    labelListList neighbourCellsForFace = neighbourCellsForEachFace(mesh);

// Find all neighbour cells for each internal face
   // List of faces
/*    const faceList& faces=mesh.faces();

    // Or for all faces?
    labelListList neighbourCellsForFace(mesh.nInternalFaces());

    forAll(faces, facei)
    {
        if(mesh.isInternalFace(facei))
        {
          labelList neighbourCells;
          labelList pointsFacei = faces[facei];

          //forAll(pointsFacei, pointi)
          for (int k=0; k<pointsFacei.size(); k++)
          {
              label pointi = pointsFacei[k];
              // cells should be added to list if it wasn't contain in the list yet
              labelList neighbourCellsForPointI = mesh.pointCells()[pointi];
              //forAll(neighbourCellsForPointI, celli)
              for(int j=0; j<neighbourCellsForPointI.size(); j++)
              {
                  label celli = neighbourCellsForPointI[j];

                  bool contained = false;

                  //forAll(neighbourCells, ncell)
                  for(int f=0; f<neighbourCells.size(); f++)
                  {
                      label ncell=neighbourCells[f];
                      if(ncell==celli)
                      {
                          contained = true;
                      }
                  }

                    if(!contained)
                    {
                        neighbourCells.append(celli);
                    }
                }
            }
            neighbourCellsForFace[facei].append(neighbourCells);
        }
    }*/
// End Find all neighbour cells for each face

// Find gradient of scalar field on the centers of internal faces
// 1. Compute vectors of distances between center of internal face and centers of its neighbour cells

    forAll(faces, facei)
    {
        vector gradF = vector::zero;

        if (mesh.isInternalFace(facei))
        {
            List<vector> df(neighbourCellsForFace[facei].size());
            scalarList wf2(neighbourCellsForFace[facei].size());
            symmTensor G(0);

            for(int i=0; i<neighbourCellsForFace[facei].size(); i++)
            {
                df[i] = mesh.cellCentres()[neighbourCellsForFace[facei][i]] - mesh.faceCentres()[facei];
                wf2[i] = 1/magSqr(df[i]);
                symmTensor addToG(0);
                addToG = sqr(df[i]);//*df[i];
                addToG = addToG * wf2[i];
                G += addToG;
            }

            symmTensor G0(0);

            if(mesh.nSolutionD()==1)
            {
                symmTensor G01(1, 0, 0, 1, 0, 1);
                symmTensor G02(sqr((Vector<label>::one + mesh.geometricD())/2));
                G0 = G01 - G02;
            }
            {
                if(mesh.nSolutionD()==2)
                {
                    G0 = sqr((Vector<label>::one - mesh.geometricD())/2);
                }
            };

            G = G + G0;
            G = inv(G);
            G = G - G0;

            for(int i=0; i<df.size(); i++)
            {
                df[i] = G&df[i];
                gradF = gradF + wf2[i]*df[i]*(iF[neighbourCellsForFace[facei][i]] - sF[facei]);//(Fi-Ff);
            }

            gradIF[facei] = gradF;
        }
    }

    return gradIF;
}

surfaceTensorField faceVectorGrad(const fvMesh& mesh, const volVectorField& iVF)
{
    surfaceVectorField gradComp0col = faceScalarGrad(mesh, iVF.component(0));
    surfaceVectorField gradComp1col = faceScalarGrad(mesh, iVF.component(1));
    surfaceVectorField gradComp2col = faceScalarGrad(mesh, iVF.component(2));

    surfaceTensorField gradIVF("divITF", 0*fvc::snGrad(iVF) * mesh.Sf() / mesh.magSf());

    forAll(mesh.faces(), facei)
    {
        if (mesh.isInternalFace(facei))
        {
            gradIVF[facei].component(0) = gradComp0col[facei].component(0);
            gradIVF[facei].component(1) = gradComp1col[facei].component(0);
            gradIVF[facei].component(2) = gradComp2col[facei].component(0);

            gradIVF[facei].component(3) = gradComp0col[facei].component(1);
            gradIVF[facei].component(4) = gradComp1col[facei].component(1);
            gradIVF[facei].component(5) = gradComp2col[facei].component(1);

            gradIVF[facei].component(6) = gradComp0col[facei].component(2);
            gradIVF[facei].component(7) = gradComp1col[facei].component(2);
            gradIVF[facei].component(8) = gradComp2col[facei].component(2);
        }
    }

    return gradIVF;
}

surfaceScalarField faceVectorDiv(const fvMesh& mesh, const volVectorField& iVF)
{
    surfaceVectorField gradComp0 = faceScalarGrad(mesh, iVF.component(0));
    surfaceVectorField gradComp1 = faceScalarGrad(mesh, iVF.component(1));
    surfaceVectorField gradComp2 = faceScalarGrad(mesh, iVF.component(2));

    surfaceScalarField divIVF("divITF", 0*fvc::snGrad(iVF) & mesh.Sf() / mesh.magSf());

    forAll(mesh.faces(), facei)
    {
        if (mesh.isInternalFace(facei))
        {
            divIVF[facei] = gradComp0[facei].component(0) + gradComp1[facei].component(1) + gradComp2[facei].component(2);
        }
    }

    return divIVF;
}

surfaceVectorField faceTensorDiv(const fvMesh& mesh, const volTensorField& iTF)
{
    surfaceVectorField gradComp0 = faceScalarGrad(mesh, iTF.component(0));
    surfaceVectorField gradComp1 = faceScalarGrad(mesh, iTF.component(1));
    surfaceVectorField gradComp2 = faceScalarGrad(mesh, iTF.component(2));

    surfaceVectorField gradComp3 = faceScalarGrad(mesh, iTF.component(3));
    surfaceVectorField gradComp4 = faceScalarGrad(mesh, iTF.component(4));
    surfaceVectorField gradComp5 = faceScalarGrad(mesh, iTF.component(5));

    surfaceVectorField gradComp6 = faceScalarGrad(mesh, iTF.component(6));
    surfaceVectorField gradComp7 = faceScalarGrad(mesh, iTF.component(7));
    surfaceVectorField gradComp8 = faceScalarGrad(mesh, iTF.component(8));

    surfaceScalarField divComp0 = gradComp0.component(0) + gradComp3.component(1) + gradComp6.component(2);
    surfaceScalarField divComp1 = gradComp1.component(0) + gradComp4.component(1) + gradComp7.component(2);
    surfaceScalarField divComp2 = gradComp2.component(0) + gradComp5.component(1) + gradComp8.component(2);

    surfaceVectorField divITF("divITF", 0*fvc::snGrad(iTF.component(0)) * mesh.Sf() / mesh.magSf());

    forAll(mesh.faces(), facei)
    {
        if (mesh.isInternalFace(facei))
        {
            divITF[facei].component(0) = divComp0[facei];
            divITF[facei].component(1) = divComp1[facei];
            divITF[facei].component(2) = divComp2[facei];
        }
    }
    return divITF;
}

