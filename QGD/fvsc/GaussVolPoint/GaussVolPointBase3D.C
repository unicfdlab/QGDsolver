
#include "GaussVolPointBase3D.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcSnGrad.H"
#include "volPointInterpolation.H"
#include "processorFvPatch.H"
#include "processorFvPatchFields.H"

Foam::fvsc::GaussVolPointBase3D::GaussVolPointBase3D(const fvMesh& mesh)
:
    volPoint_
    (
        volPointInterpolation::New(mesh)
    ),
    nfRef_(mesh.thisDb().lookupObject<surfaceVectorField>("nf")),
    bgfid_(mesh.boundary().size()),
    processorPatch_(mesh.boundary().size(), false),
    qf_(0),
    aqx_(0),
    aqy_(0),
    aqz_(0),
    vq_(0),
    bqf_(mesh.boundary().size()),
    baqx_(mesh.boundary().size()),
    baqy_(mesh.boundary().size()),
    baqz_(mesh.boundary().size()),
    bvq_(mesh.boundary().size()),
    bmv65_(mesh.boundary().size()),
    of_(0),
    bof_(mesh.boundary().size())
{
    forAll(bgfid_, iPatch)
    {
        if (isA<processorFvPatch>(mesh.boundary()[iPatch]))
        {
            processorPatch_[iPatch] = true;
        }
        const fvPatch& fvp = mesh.boundary()[iPatch];
        bgfid_[iPatch].resize(fvp.size());
        forAll(fvp, i)
        {
            bgfid_[iPatch][i] = mesh.boundary()[iPatch].start() + i;
        }
    }
    
    //sort quad faces and other
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();
    
    forAll(faces, i)
    {
        if (mesh.isInternalFace(i))
        {
            if (faces[i].size() == 4)
            {
                qf_.append(i);
            }
            else
            {
                of_.append(i);
            }
        }
    }
    forAll(bgfid_, iPatch)
    {
        const fvPatch& fvp = mesh.boundary()[iPatch];
        forAll(fvp, i)
        {
            if (faces[bgfid_[iPatch][i]].size() == 4)
            {
                bqf_[iPatch].append(i);
            }
            else
            {
                bof_[iPatch].append(i);
            }
        }
    }
    
    //calculate weights
    label facei = -1;
    aqx_.resize(qf_.size());
    aqy_.resize(qf_.size());
    aqz_.resize(qf_.size());
    vq_.resize(qf_.size());
    label own = -1;
    label nei = -1;
    label p1 = -1, p2 = -1, p3 = -1, p4 = -1;
    const scalar OneBySix = (1.0 / 6.0);
    forAll(qf_, i)
    {
        facei = qf_[i];
        own = mesh.owner()[facei];
        nei = mesh.neighbour()[facei];
        p1  = faces[facei][0];
        p2  = faces[facei][1];
        p3  = faces[facei][2];
        p4  = faces[facei][3];
        //p5 - is nei
        //p6 - is own
        
        vq_[i] = (points[p3] - points[p1]) & 
        (
            (points[p4] - points[p2]) ^ (mesh.C()[own] - mesh.C()[nei])
        );
        vq_[i] *= OneBySix;
        
        /* coefficient for X */
        aqx_[i].resize(6);
        aqx_[i][0] = OneBySix*((mesh.C()[nei].y() - mesh.C()[own].y())*(points[p2].z() - points[p4].z())
            - (mesh.C()[nei].z() - mesh.C()[own].z())*(points[p2].y() - points[p4].y()));
        aqx_[i][1] = OneBySix*((mesh.C()[nei].y() - mesh.C()[own].y())*(points[p3].z() - points[p1].z())
            - (mesh.C()[nei].z() - mesh.C()[own].z())*(points[p3].y() - points[p1].y()));
        aqx_[i][5] = OneBySix*((points[p1].y() - points[p3].y())*(points[p2].z() - points[p4].z())
            - (points[p1].z() - points[p3].z())*(points[p2].y() - points[p4].y()));
        
        aqx_[i][2] = -aqx_[i][0];
        aqx_[i][3] = -aqx_[i][1];
        aqx_[i][4] = -aqx_[i][5];
        
        /* coefficient for Y */
        aqy_[i].resize(6);
        aqy_[i][0] = OneBySix*((mesh.C()[nei].z() - mesh.C()[own].z())*(points[p2].x() - points[p4].x())
            - (mesh.C()[nei].x() - mesh.C()[own].x())*(points[p2].z() - points[p4].z()));
        aqy_[i][1] = OneBySix*((mesh.C()[nei].z() - mesh.C()[own].z())*(points[p3].x() - points[p1].x())
            - (mesh.C()[nei].x() - mesh.C()[own].x())*(points[p3].z() - points[p1].z()));
        aqy_[i][5] = OneBySix*((points[p1].z() - points[p3].z())*(points[p2].x() - points[p4].x())
            - (points[p1].x() - points[p3].x())*(points[p2].z() - points[p4].z()));
        
        aqy_[i][2] = -aqy_[i][0];
        aqy_[i][3] = -aqy_[i][1];
        aqy_[i][4] = -aqy_[i][5];
        
        /* coefficient for Z */
        aqz_[i].resize(6);
        aqz_[i][0] = OneBySix*((mesh.C()[nei].x() - mesh.C()[own].x())*(points[p2].y() - points[p4].y())
            - (mesh.C()[nei].y() - mesh.C()[own].y())*(points[p2].x() - points[p4].x()));
        aqz_[i][1] = OneBySix*((mesh.C()[nei].x() - mesh.C()[own].x())*(points[p3].y() - points[p1].y())
            - (mesh.C()[nei].y() - mesh.C()[own].y())*(points[p3].x() - points[p1].x()));
        aqz_[i][5] = OneBySix*((points[p1].x() - points[p3].x())*(points[p2].y() - points[p4].y())
            - (points[p1].y() - points[p3].y())*(points[p2].x() - points[p4].x()));
        
        aqz_[i][2] = -aqz_[i][0];
        aqz_[i][3] = -aqz_[i][1];
        aqz_[i][4] = -aqz_[i][5];
    }
    forAll(bqf_, iPatch)
    {
        baqx_[iPatch].resize(bqf_[iPatch].size());
        baqy_[iPatch].resize(bqf_[iPatch].size());
        baqz_[iPatch].resize(bqf_[iPatch].size());
        bvq_[iPatch].resize(bqf_[iPatch].size());
        
        vectorField v6(mesh.boundary()[iPatch].size(), vector::zero);
        vectorField v5(mesh.boundary()[iPatch].size(), vector::zero);
        
        v6 = mesh.boundary()[iPatch].Cn();
        if (processorPatch_[iPatch])
        {
            v5 = refCast<const processorFvPatch>(mesh.boundary()[iPatch]).
                    procPolyPatch().neighbFaceCellCentres();
        }
        else
        {
            v5 = v6 + 2.0*
            (
                mesh.boundary()[iPatch].Cf()
                -
                v6
            );
        }
        
        bmv65_[iPatch] = mag
        (
            v5 - v6
        );
        label gFaceId = -1;
        forAll(bqf_[iPatch], k)
        {
            facei = bqf_[iPatch][k];
            gFaceId = bgfid_[iPatch][facei];
            
            p1  = faces[gFaceId][0];
            p2  = faces[gFaceId][1];
            p3  = faces[gFaceId][2];
            p4  = faces[gFaceId][3];
            //p5 - is nei and stored in v5
            //p6 - is own and stored in v6
            
            bvq_[iPatch][k] = (points[p3] - points[p1]) & 
            (
                (points[p4] - points[p2]) ^ (v6[facei] - v5[facei])
            );
            bvq_[iPatch][k] *= OneBySix;
            
            /* coefficient for X */
            baqx_[iPatch][k].resize(6);
            baqx_[iPatch][k][0] = OneBySix*((v5[facei].y() - v6[facei].y())*(points[p2].z() - points[p4].z())
                - (v5[facei].z() - v6[facei].z())*(points[p2].y() - points[p4].y()));
            baqx_[iPatch][k][1] = OneBySix*((v5[facei].y() - v6[facei].y())*(points[p3].z() - points[p1].z())
                - (v5[facei].z() - v6[facei].z())*(points[p3].y() - points[p1].y()));
            baqx_[iPatch][k][5] = OneBySix*((points[p1].y() - points[p3].y())*(points[p2].z() - points[p4].z())
                - (points[p1].z() - points[p3].z())*(points[p2].y() - points[p4].y()));
            
            baqx_[iPatch][k][2] = -baqx_[iPatch][k][0];
            baqx_[iPatch][k][3] = -baqx_[iPatch][k][1];
            baqx_[iPatch][k][4] = -baqx_[iPatch][k][5];

            /* coefficient for Y */
            baqy_[iPatch][k].resize(6);
            baqy_[iPatch][k][0] = OneBySix*((v5[facei].z() - v6[facei].z())*(points[p2].x() - points[p4].x())
                - (v5[facei].x() - v6[facei].x())*(points[p2].z() - points[p4].z()));
            baqy_[iPatch][k][1] = OneBySix*((v5[facei].z() - v6[facei].z())*(points[p3].x() - points[p1].x())
                - (v5[facei].x() - v6[facei].x())*(points[p3].z() - points[p1].z()));
            baqy_[iPatch][k][5] = OneBySix*((points[p1].z() - points[p3].z())*(points[p2].x() - points[p4].x())
                - (points[p1].x() - points[p3].x())*(points[p2].z() - points[p4].z()));
            
            baqy_[iPatch][k][2] = -baqy_[iPatch][k][0];
            baqy_[iPatch][k][3] = -baqy_[iPatch][k][1];
            baqy_[iPatch][k][4] = -baqy_[iPatch][k][5];
            
            /* coefficient for Z */
            baqz_[iPatch][k].resize(6);
            baqz_[iPatch][k][0] = OneBySix*((v5[facei].x() - v6[facei].x())*(points[p2].y() - points[p4].y())
                - (v5[facei].y() - v6[facei].y())*(points[p2].x() - points[p4].x()));
            baqz_[iPatch][k][1] = OneBySix*((v5[facei].x() - v6[facei].x())*(points[p3].y() - points[p1].y())
                - (v5[facei].y() - v6[facei].y())*(points[p3].x() - points[p1].x()));
            baqz_[iPatch][k][5] = OneBySix*((points[p1].x() - points[p3].x())*(points[p2].y() - points[p4].y())
                - (points[p1].y() - points[p3].y())*(points[p2].x() - points[p4].x()));
            
            baqz_[iPatch][k][2] = -baqz_[iPatch][k][0];
            baqz_[iPatch][k][3] = -baqz_[iPatch][k][1];
            baqz_[iPatch][k][4] = -baqz_[iPatch][k][5];
        }
    }
};


Foam::fvsc::GaussVolPointBase3D::~GaussVolPointBase3D()
{
}

#define VEC_CMPT(V,CMPT)\
    V[CMPT]

#define SCA_CMPT(V,CMPT)\
    V

#define dfdxif(vf,pf,dfdxfield,aqi,icmpt,ocmpt,iop,oop)                 \
{                                                                       \
    label celll = -1;                                                   \
    label facei = -1;                                                   \
    scalar dfdxface = 0.0;                                              \
    forAll(qf_, i)                                                      \
    {                                                                   \
        facei = qf_[i];                                                 \
        celll = vf.mesh().neighbour()[facei];                           \
        dfdxface =                                                      \
            iop(vf.primitiveField()[celll],icmpt) * aqi[i][4];          \
        celll = vf.mesh().owner()[facei];                               \
        dfdxface +=                                                     \
            iop(vf.primitiveField()[celll],icmpt) * aqi[i][5];          \
                                                                        \
        forAll(faces[facei], k)                                         \
        {                                                               \
            dfdxface +=                                                 \
                iop(pf[faces[facei][k]],icmpt) * aqi[i][k];             \
        }                                                               \
        oop(dfdxfield.primitiveFieldRef()[facei],ocmpt)                 \
            += (dfdxface / vq_[i]);                                     \
    }                                                                   \
}

#define dfdxbf(vf,pf,patchi,dfdxfield,baqi,icmpt,ocmpt,iop,oop)         \
{                                                                       \
    label qfacei = -1;                                                  \
    label gFaceId = -1;                                                 \
    scalar dfdxface = 0.0;                                              \
    forAll(bqf_[patchi], i)                                             \
    {                                                                   \
        qfacei = bqf_[patchi][i];                                       \
        gFaceId = bgfid_[patchi][qfacei];                               \
        dfdxface =                                                      \
            iop(psi5[patchi][qfacei],icmpt) * baqi[patchi][i][4];       \
        dfdxface +=                                                     \
            iop(psi6[qfacei],icmpt) * baqi[patchi][i][5];               \
        forAll(faces[gFaceId], k)                                       \
        {                                                               \
            dfdxface+=                                                  \
                    iop(pf[faces[gFaceId][k]],icmpt) *                  \
                    baqi[patchi][i][k];                                 \
        }                                                               \
        oop(dfdxfield.boundaryFieldRef()[patchi][qfacei],ocmpt) +=      \
        dfdxface / bvq_[patchi][i];                                     \
    }                                                                   \
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfIF
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceScalarField& divf,
    const surfaceScalarField& dfdn
)
{
    dfdxif(vf,pf,divf,aqx_,0,0,VEC_CMPT,SCA_CMPT) //X
    dfdxif(vf,pf,divf,aqy_,1,0,VEC_CMPT,SCA_CMPT) //Y
    dfdxif(vf,pf,divf,aqz_,2,0,VEC_CMPT,SCA_CMPT) //Z

    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            divf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfBF
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceScalarField& divf,
    const surfaceScalarField& dfdn
)
{
    List<List<vector> > psi5 (vf.boundaryField().size());
    forAll(psi5, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psi5[iPatch] = refCast<const processorFvPatchField<vector> >
                (vf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psi5[iPatch] = vf.boundaryField()[iPatch] + 
                vf.boundaryField()[iPatch].snGrad()
                *
                bmv65_[iPatch]*0.5;
        }
        
        vectorField psi6 (vf.boundaryField()[iPatch].patchInternalField());
        
        dfdxbf(vf,pf,iPatch,divf,baqx_,0,0,VEC_CMPT,SCA_CMPT) //X
        dfdxbf(vf,pf,iPatch,divf,baqy_,1,0,VEC_CMPT,SCA_CMPT) //Y
        dfdxbf(vf,pf,iPatch,divf,baqz_,2,0,VEC_CMPT,SCA_CMPT) //Z
        
        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                divf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfIF
(
    const volTensorField& tf,
    const pointTensorField& pf,
    const faceList& faces,
    surfaceVectorField& divf,
    const surfaceVectorField& dfdn
)
{
    //X
    dfdxif(tf,pf,divf,aqx_,0,0,VEC_CMPT,VEC_CMPT) // dT_xx / dx
    dfdxif(tf,pf,divf,aqy_,3,0,VEC_CMPT,VEC_CMPT) // dT_yx / dy
    dfdxif(tf,pf,divf,aqz_,6,0,VEC_CMPT,VEC_CMPT) // dT_zx / dz
    //Y
    dfdxif(tf,pf,divf,aqx_,1,1,VEC_CMPT,VEC_CMPT) // dT_xy / dx
    dfdxif(tf,pf,divf,aqy_,4,1,VEC_CMPT,VEC_CMPT) // dT_yy / dy
    dfdxif(tf,pf,divf,aqz_,7,1,VEC_CMPT,VEC_CMPT) // dT_zy / dz
    //Z
    dfdxif(tf,pf,divf,aqx_,2,2,VEC_CMPT,VEC_CMPT) // dT_xz / dx
    dfdxif(tf,pf,divf,aqy_,5,2,VEC_CMPT,VEC_CMPT) // dT_yz / dy
    dfdxif(tf,pf,divf,aqz_,8,2,VEC_CMPT,VEC_CMPT) // dT_zz / dz
    //
    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            divf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcDivfBF
(
    const volTensorField& tf,
    const pointTensorField& pf,
    const faceList& faces,
    surfaceVectorField& divf,
    const surfaceVectorField& dfdn
)
{
    List<List<tensor> > psi5 (tf.boundaryField().size());
    forAll(psi5, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psi5[iPatch] = refCast<const processorFvPatchField<tensor> >
                (tf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psi5[iPatch] = tf.boundaryField()[iPatch] + 
                tf.boundaryField()[iPatch].snGrad()
                *
                bmv65_[iPatch]*0.5;
        }
        
        tensorField psi6 (tf.boundaryField()[iPatch].patchInternalField());
        
        dfdxbf(tf,pf,iPatch,divf,baqx_,0,0,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,baqy_,3,0,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,baqz_,6,0,VEC_CMPT,VEC_CMPT) //d/dz
        
        dfdxbf(tf,pf,iPatch,divf,baqx_,1,1,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,baqy_,4,1,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,baqz_,7,1,VEC_CMPT,VEC_CMPT) //d/dz
        
        dfdxbf(tf,pf,iPatch,divf,baqx_,2,2,VEC_CMPT,VEC_CMPT) //d/dx
        dfdxbf(tf,pf,iPatch,divf,baqy_,5,2,VEC_CMPT,VEC_CMPT) //d/dy
        dfdxbf(tf,pf,iPatch,divf,baqz_,8,2,VEC_CMPT,VEC_CMPT) //d/dz
        
        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                divf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfIF
(
    const volScalarField& sf,
    const pointScalarField& pf,
    const faceList& faces,
    surfaceVectorField& gradf,
    const surfaceVectorField& dfdn
)
{
    //quad faces
    dfdxif(sf,pf,gradf,aqx_,0,0,SCA_CMPT,VEC_CMPT) //X
    dfdxif(sf,pf,gradf,aqy_,0,1,SCA_CMPT,VEC_CMPT) //Y
    dfdxif(sf,pf,gradf,aqz_,0,2,SCA_CMPT,VEC_CMPT) //Z
    
    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            gradf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfBF
(
    const volScalarField& sf,
    const pointScalarField& pf,
    const faceList& faces,
    surfaceVectorField& gradf,
    const surfaceVectorField& dfdn
)
{
    List<List<scalar> > psi5 (sf.boundaryField().size());
    forAll(psi5, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psi5[iPatch] = refCast<const processorFvPatchField<scalar> >
                (sf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psi5[iPatch] = sf.boundaryField()[iPatch] + 
                sf.boundaryField()[iPatch].snGrad()
                *
                bmv65_[iPatch]*0.5;
        }
        
        scalarField psi6 (sf.boundaryField()[iPatch].patchInternalField());
        
        dfdxbf(sf,pf,iPatch,gradf,baqx_,0,0,SCA_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,baqy_,0,1,SCA_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,baqz_,0,2,SCA_CMPT,VEC_CMPT) //Z
        
        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                gradf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfIF
(
    const volVectorField& vf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceTensorField& gradf,
    const surfaceTensorField& dfdn
)
{
    //quad faces
    dfdxif(vf,pf,gradf,aqx_,0,0,VEC_CMPT,VEC_CMPT) //X
    dfdxif(vf,pf,gradf,aqx_,1,1,VEC_CMPT,VEC_CMPT) //Y
    dfdxif(vf,pf,gradf,aqx_,2,2,VEC_CMPT,VEC_CMPT) //Z
    
    dfdxif(vf,pf,gradf,aqy_,0,3,VEC_CMPT,VEC_CMPT) //X
    dfdxif(vf,pf,gradf,aqy_,1,4,VEC_CMPT,VEC_CMPT) //Y
    dfdxif(vf,pf,gradf,aqy_,2,5,VEC_CMPT,VEC_CMPT) //Z
    
    dfdxif(vf,pf,gradf,aqz_,0,6,VEC_CMPT,VEC_CMPT) //X
    dfdxif(vf,pf,gradf,aqz_,1,7,VEC_CMPT,VEC_CMPT) //Y
    dfdxif(vf,pf,gradf,aqz_,2,8,VEC_CMPT,VEC_CMPT) //Z

    //other faces
    {
        label facei = -1;
        forAll(of_, i)
        {
            facei = of_[i];
            gradf.primitiveFieldRef()[facei] = 
                dfdn.primitiveField()[facei];
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::calcGradfBF
(
    const volVectorField& sf,
    const pointVectorField& pf,
    const faceList& faces,
    surfaceTensorField& gradf,
    const surfaceTensorField& dfdn
)
{
    List<List<vector> > psi5 (sf.boundaryField().size());
    forAll(psi5, iPatch)
    {
        if (processorPatch_[iPatch])
        {
            psi5[iPatch] = refCast<const processorFvPatchField<vector> >
                (sf.boundaryField()[iPatch]).patchNeighbourField();
        }
        else
        {
            psi5[iPatch] = sf.boundaryField()[iPatch] + 
                sf.boundaryField()[iPatch].snGrad()
                *
                bmv65_[iPatch]*0.5;
        }
        
        vectorField psi6 (sf.boundaryField()[iPatch].patchInternalField());
        
        dfdxbf(sf,pf,iPatch,gradf,baqx_,0,0,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,baqx_,1,1,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,baqx_,2,2,VEC_CMPT,VEC_CMPT) //Z

        dfdxbf(sf,pf,iPatch,gradf,baqy_,0,3,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,baqy_,1,4,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,baqy_,2,5,VEC_CMPT,VEC_CMPT) //Z

        dfdxbf(sf,pf,iPatch,gradf,baqz_,0,6,VEC_CMPT,VEC_CMPT) //X
        dfdxbf(sf,pf,iPatch,gradf,baqz_,1,7,VEC_CMPT,VEC_CMPT) //Y
        dfdxbf(sf,pf,iPatch,gradf,baqz_,2,8,VEC_CMPT,VEC_CMPT) //Z

        //for other faces - apply surface normal derivative
        {
            label facei = -1;
            forAll(bof_[iPatch], i)
            {
                facei = bof_[iPatch][i];
                gradf.boundaryFieldRef()[iPatch][facei] = 
                    dfdn.boundaryField()[iPatch][facei];
            }
        }
    }
}

void Foam::fvsc::GaussVolPointBase3D::faceGrad(const volScalarField& sf, surfaceVectorField& gradf)
{
    pointScalarField pF
    (
        volPoint_.interpolate
        (
            sf
        )
    );
    //
    const surfaceVectorField& nf = nfRef_();
    surfaceVectorField dfdn
    (
        nf * fvc::snGrad(sf)
    );
    //
    const faceList& faces = sf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcGradfIF(sf, pF, faces, gradf, dfdn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcGradfBF(sf, pF, faces, gradf, dfdn);
};

void Foam::fvsc::GaussVolPointBase3D::faceGrad(const volVectorField& vf, surfaceTensorField& gradf)
{
    pointVectorField pF
    (
        volPoint_.interpolate
        (
            vf
        )
    );
    //
    const surfaceVectorField& nf = nfRef_();
    surfaceTensorField dfdn
    (
        nf * fvc::snGrad(vf)
    );
    //
    const faceList& faces = vf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcGradfIF(vf, pF, faces, gradf, dfdn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcGradfBF(vf, pF, faces, gradf, dfdn);
};

void Foam::fvsc::GaussVolPointBase3D::faceDiv(const volVectorField& vf, surfaceScalarField& divf)
{
    pointVectorField pF
    (
        volPoint_.interpolate
        (
            vf
        )
    );
    
    surfaceScalarField divfn (nfRef_() & fvc::snGrad(vf));
    const faceList& faces = vf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcDivfIF(vf, pF, faces, divf, divfn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcDivfBF(vf, pF, faces, divf, divfn);
};

void Foam::fvsc::GaussVolPointBase3D::faceDiv(const volTensorField& tf, surfaceVectorField& divf)
{
    pointTensorField pF
    (
        volPoint_.interpolate
        (
            tf
        )
    );
    
    surfaceVectorField divfn (nfRef_() & fvc::snGrad(tf));
    const faceList& faces = tf.mesh().faces();
    /*
     *
     * Calculate grad at internal faces
     *
     */
    calcDivfIF(tf, pF, faces, divf, divfn);
    /*
     *
     * Calculate grad at boundary faces
     *
     */
    calcDivfBF(tf, pF, faces, divf, divfn);
}

//
//END-OF-FILE
//


