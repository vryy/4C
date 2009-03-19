/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_evaluate.cpp
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef D_SOLID3
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_sh8p8.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "Epetra_Time.h"
#include "Teuchos_TimeMonitor.hpp"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices6VoigtTo2Tensor(
  const int*& voigt6row,
  const int*& voigt6col,
  const bool transpose
  )
{
  const int Voigt6Row[NUMSTR_] = {0,1,2, 0,1,2};
  const int Voigt6Col[NUMSTR_] = {0,1,2, 1,2,0};

  if (transpose)
  {
    voigt6row = &(Voigt6Col[0]);
    voigt6col = &(Voigt6Row[0]);
  }
  else
  {
    voigt6row = &(Voigt6Row[0]);
    voigt6col = &(Voigt6Col[0]);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices9VoigtTo2Tensor(
  const int*& voigt9row,
  const int*& voigt9col,
  const bool transpose
  )
{
  // 9-Voigt C-index                      0 1 2  3 4 5  6 7 8
  const int Voigt9Row[NUMDFGR_] = {0,1,2, 0,1,2, 0,2,1};
  const int Voigt9Col[NUMDFGR_] = {0,1,2, 1,2,0, 2,1,0};

  if (transpose)
  {
    voigt9row = &(Voigt9Col[0]);
    voigt9col = &(Voigt9Row[0]);
  }
  else
  {
    voigt9row = &(Voigt9Row[0]);
    voigt9col = &(Voigt9Col[0]);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices2TensorTo9Voigt(
  const int*& voigt3x3
  )
{
  // tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
  // C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
  // Access : 3*i+j
  // 9-Voigt C-indices    0   3   6   8   1   4   5   7   2
  const int Voigt3x3[NUMDFGR_] = {0,3,6, 8,1,4, 5,7,2};

  voigt3x3 = &(Voigt3x3[0]);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Indices2TensorTo6Voigt(
  const int*& voigt3x3
  )
{
  // tensor indices ij = 11, 12, 13, 21, 22, 23, 31, 32, 33
  // C indices           00, 01, 02, 10, 11, 12, 20, 21, 22
  // Access : 3*i+j
  // 6-Voigt C-indices    0   3   5   3   1   4   5   4   2
  const int Voigt3x3[NUMDFGR_] = {0,3,5, 3,1,4, 5,4,2};

  voigt3x3 = &(Voigt3x3[0]);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector9Voigt(
  LINALG::Matrix<NUMDFGR_,1>& fvct,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,
  const bool transpose
  )
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col,transpose);
    
  for (int ij=0; ij<NUMDFGR_; ++ij)
  {
    const int i = voigt9row[ij];
    const int j = voigt9col[ij];
    fvct(ij,0) = fmat(i,j);  // F_ij
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector6Voigt(
  LINALG::Matrix<NUMSTR_,1>& bvct,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& bmat,
  const VoigtType outvoigt6
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);
    
  for (int ij=0; ij<NUMSTR_; ++ij)
  {
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    if (ij < NUMDIM_)
      bvct(ij) = bmat(i,j);  // B_ij
    else
      if (outvoigt6 == voigt6_strain)
        bvct(ij) = bmat(i,j) + bmat(j,i);  // B_ij+B_ji
      else
        bvct(ij) = bmat(i,j);  // B_ij
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Vector6VoigtToMatrix2Tensor(
  LINALG::Matrix<NUMDIM_,NUMDIM_>& bmat,
  const LINALG::Matrix<NUMSTR_,1>& bvct,
  const VoigtType invoigt6
  )
{
  const int* voigt3x3sym = NULL;
  Indices2TensorTo6Voigt(voigt3x3sym);  // access is via (i,j) -> 3*i+j  

  for (int i=0; i<NUMDIM_; ++ i)
  {
    for (int j=0; j<NUMDIM_; ++j)
    {
      const int ij = voigt3x3sym[NUMDIM_*i+j];
      if (i == j)
        bmat(i,j) = bvct(ij);
      else
        if (invoigt6 == voigt6_strain)
          bmat(i,j) = 0.5*bvct(ij);
        else
          bmat(i,j) = bvct(ij);
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector9VoigtDiffByItself(
  LINALG::Matrix<NUMDFGR_,NUMDFGR_>& invfderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& invfmat,
  const bool transpose
  )
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col);

  // VERIFIED

  for (int ij=0; ij<NUMDFGR_; ++ij)
  {
    const int i = voigt9row[ij];
    const int j = voigt9col[ij];
    for (int kl=0; kl<NUMDFGR_; ++kl)
    {
      const int k = voigt9row[kl];
      const int l = voigt9col[kl];
      if (transpose)
        invfderf(ij,kl) = -invfmat(j,k)*invfmat(l,i);
      else
        invfderf(ij,kl) = -invfmat(i,k)*invfmat(l,j);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_,NUMSTR_>& invfderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& invfmat
  )
{
  const int voigt6row[NUMSTR_] = {0,1,2, 0,1,2};
  const int voigt6col[NUMSTR_] = {0,1,2, 1,2,0};

  // VERIFIED

//  cout << endl;
  for (int ij=0; ij<NUMSTR_; ++ij)
  {
//    cout << "[";
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMSTR_; ++kl)
    {
      const int k = voigt6row[kl];
      const int l = voigt6col[kl];
      invfderf(ij,kl) = -0.5*(invfmat(i,k)*invfmat(l,j) + invfmat(i,l)*invfmat(k,j));
//     cout << "ct["<<i+1<<","<<k+1<<"]*ct["<<l+1<<","<<j+1<<"]+ct["<<i+1<<","<<l+1<<"]*ct["<<k+1<<","<<j+1<<"]";
      if (ij >= NUMDIM_)
      {
        invfderf(ij,kl) += -0.5*(invfmat(j,k)*invfmat(l,i) + invfmat(j,l)*invfmat(k,i));
//        cout << "+ct["<<j+1<<","<<k+1<<"]*ct["<<l+1<<","<<i+1<<"]+ct["<<j+1<<","<<l+1<<"]*ct["<<k+1<<","<<i+1<<"]";
      }
//      cout << ", ";
    }
//    cout << "]," << endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector6VoigtTwiceDiffByItself(
  LINALG::Matrix<NUMSTR_,NUMSTR_*NUMSTR_>& invbvdderb,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& ibt
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  for (int ij=0; ij<NUMSTR_; ++ij)
  {
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMSTR_; ++kl)
    {
      const int k = voigt6row[kl];
      const int l = voigt6col[kl];
      for (int mn=0; mn<NUMSTR_; ++mn)
      {
        const int m = voigt6row[mn];
        const int n = voigt6col[mn];
        const int klmn = NUMSTR_*kl + mn;
        invbvdderb(ij,klmn) = 0.25*(
          ( ibt(i,m)*ibt(n,k) + ibt(i,n)*ibt(m,k) )*ibt(l,j)
          + ibt(i,k)*( ibt(l,m)*ibt(n,j) + ibt(l,n)*ibt(m,j) )
          + ( ibt(i,m)*ibt(n,l) + ibt(i,n)*ibt(m,l) )*ibt(k,j)
          + ibt(i,l)*( ibt(k,m)*ibt(n,j) + ibt(k,n)*ibt(m,j) )
          );
        if (ij >= NUMDIM_)  // swap i 'n' j 
          invbvdderb(ij,klmn) += 0.25*(
            ( ibt(j,m)*ibt(n,k) + ibt(j,n)*ibt(m,k) )*ibt(l,i)
            + ibt(j,k)*( ibt(l,m)*ibt(n,i) + ibt(l,n)*ibt(m,i) )
            + ( ibt(j,m)*ibt(n,l) + ibt(j,n)*ibt(m,l) )*ibt(k,i)
            + ibt(j,l)*( ibt(k,m)*ibt(n,i) + ibt(k,n)*ibt(m,i) )
            );
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_,NUMSTR_>& sqfderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  // VERIFIED

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // (F.F)_{,F} with F^T=F
//  cout << endl;
  for (int ij=0; ij<NUMSTR_; ++ij)
  {
//    cout << "[";
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMSTR_; ++kl)
    {
      const int k = voigt6row[kl];
      const int l = voigt6col[kl];
      sqfderf(ij,kl) = id(i,k)*fmat(l,j) + id(j,l)*fmat(i,k);
//      cout << "id["<<i+1<<","<<k+1<<"]*St["<<l+1<<","<<j+1<<"]+id["<<j+1<<","<<l+1<<"]*St["<<i+1<<","<<k+1<<"]";
      if (ij >= NUMDIM_)
      {
        sqfderf(ij,kl) += id(j,k)*fmat(l,i) + id(i,l)*fmat(j,k);
//        cout << "+id["<<j+1<<","<<k+1<<"]*St["<<l+1<<","<<i+1<<"]+id["<<i+1<<","<<l+1<<"]*St["<<j+1<<","<<k+1<<"]";
      }
//      cout << ", ";
    }
//    cout << "]," << endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector9VoigtDiffByItself(
  LINALG::Matrix<NUMDFGR_,NUMDFGR_>& sqfderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,
  const bool transpose
  )
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col);

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // (F^T.F)_{,F}
//  cout << endl;
  for (int ij=0; ij<NUMDFGR_; ++ij)
  {
//    cout << "[";
    const int i = voigt9row[ij];
    const int j = voigt9col[ij];
    for (int kl=0; kl<NUMDFGR_; ++kl)
    {
      const int k = voigt9row[kl];
      const int l = voigt9col[kl];
//      cout << "i=" << i << ", j=" << j << ", k=" << k << ", l=" << l << endl;
      if (transpose)  // swap indices of fmat
        sqfderf(ij,kl) = id(i,k)*fmat(j,l) + id(j,l)*fmat(k,i);
      else
        sqfderf(ij,kl) = id(i,k)*fmat(l,j) + id(j,l)*fmat(i,k);
//      cout << "id["<<i+1<<","<<k+1<<"]*St["<<l+1<<","<<j+1<<"]+id["<<j+1<<","<<l+1<<"]*St["<<i+1<<","<<k+1<<"]";
//      cout << ", ";
    }
//    cout << "]," << endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtTwiceDiffByItself(
  LINALG::Matrix<NUMSTR_,NUMSTR_*NUMSTR_>& sqfdderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // (F^T.F)_{,FF} with F^T=F
  for (int ij=0; ij<NUMSTR_; ++ij)
  {
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMSTR_; ++kl)
    {
      const int k = voigt6row[kl];
      const int l = voigt6col[kl];
      for (int mn=0; mn<NUMSTR_; ++mn)
      {
        const int m = voigt6row[mn];
        const int n = voigt6col[mn];
        const int klmn = NUMSTR_*kl + mn;
        sqfdderf(ij,klmn) = id(k,m)*id(i,n)*id(j,l) + id(i,l)*id(k,m)*id(j,n);
        if (ij >= NUMDIM_)
          sqfdderf(ij,klmn) += id(k,m)*id(j,n)*id(i,l) + id(j,l)*id(k,m)*id(i,n);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToMatrix6x9Voigt(
  LINALG::Matrix<NUMSTR_,NUMDFGR_>& bm,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& bt,
  const bool transpose
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col);

  for (int ij=0; ij<NUMSTR_; ++ij)
  {
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMDFGR_; ++kl)
    {
      const int k = voigt9row[kl];
      const int l = voigt9col[kl];
      if (j == l)
        bm(ij,kl) = (transpose) ? bt(k,i) : bt(i,k);
      else if ( (ij >= NUMDIM_) and (i == l) )
        bm(ij,kl) = (transpose) ? bt(k,j) : bt(j,k);
      else
        bm(ij,kl) = 0.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToLeftRightProductMatrix6x6Voigt(
  LINALG::Matrix<NUMSTR_,NUMSTR_>& bm,  ///< (out) 6x6 Voigt matrix
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& bt,  ///< (in) 3x3 matrix of 2-tensor
  const bool transpose, ///< 3x3 input matrix is transposed
  const VoigtType outvoigt6,  ///< 6-Voigt vector layout on rows of 6x6 matrix
  const VoigtType invoigt6  ///< 6-Voigt vector layout on columns of 6x6 matrix
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  for (int ab=0; ab<NUMSTR_; ++ab)
  {
    const int a = voigt6row[ab];
    const int b = voigt6col[ab];
    for (int AB=0; AB<NUMSTR_; ++AB)
    {
      const int A = voigt6row[AB];
      const int B = voigt6col[AB];
      if (transpose)
      {
        bm(AB,ab) = bt(A,a)*bt(B,b);
        if (ab >= NUMSTR_) bm(AB,ab) += bt(A,b)*bt(B,a);
      }
      else
      {
        bm(AB,ab) = bt(a,A)*bt(b,B);
        if (ab >= NUMSTR_) bm(AB,ab) += bt(b,A)*bt(a,B);
      }
      if ( (invoigt6 == voigt6_stress) and (ab >= NUMSTR_) )
        bm(AB,ab) *= 2.0;
      if ( (outvoigt6 == voigt6_stress) and (AB >= NUMSTR_) )
        bm(AB,ab) *= 0.5;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::ExtractDispAndPres(
  std::vector<double>& mystat,
  LINALG::Matrix<NUMDISP_,1>& mydisp,
  LINALG::Matrix<NUMPRES_,1>& mypres
  )
{
  for (int inod=0; inod<NUMNOD_; ++inod)
  {
    for (int idis=0; idis<NODDISP_; ++idis)
      mydisp(idis+(inod*NODDISP_),0) = mystat[idis+(inod*NODDOF_)];
    mypres(inod,0) = mystat[NODDISP_+0+(inod*NODDOF_)];
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementMatrix(
  LINALG::Matrix<NUMDOF_,NUMDOF_>* mat,
  const LINALG::Matrix<NUMDISP_,NUMDISP_>* matdd,
  const LINALG::Matrix<NUMDISP_,NUMPRES_>* matdp,
  const LINALG::Matrix<NUMDISP_,NUMPRES_>* matpd,
  const LINALG::Matrix<NUMPRES_,NUMPRES_>* matpp
)
{
  const int d2dp[NUMDISP_] = {0,1,2,  4,5,6,  8,9,10,   12,13,14,   16,17,18,   20,21,22,   24,25,26,   28,29,30  };
  const int p2dp[NUMPRES_] = {      3,      7,       11,         15,         19,         23,         27,        31};
  for (int i=0; i<NUMDISP_; ++i)
  {
    const int I = d2dp[i];
    for (int j=0; j<NUMDISP_; ++j)
    {
      const int J = d2dp[j];
      if (matdd != NULL)
        (*mat)(I,J) = (*matdd)(i,j);
      else
        (*mat)(I,J) = 0.0;
    }

    for (int l=0; l<NUMPRES_; ++l)
    {
      const int L = p2dp[l];
      if (matdp != NULL)
        (*mat)(I,L) = (*matdp)(i,l);
      else
        (*mat)(I,L) = 0.0;
    }
  }
  for (int k=0; k<NUMPRES_; ++k)
  {
    const int K = p2dp[k];
    for (int j=0; j<NUMDISP_; ++j)
    {
      const int J = d2dp[j];
      if (matpd != NULL)
        (*mat)(K,J) = (*matpd)(k,j);
      else if (matdp != NULL)
        (*mat)(K,J) = (*matdp)(j,k);
      else
        (*mat)(K,J) = 0.0;
    }
    for (int l=0; l<NUMPRES_; ++l)
    {
      const int L = p2dp[l];
      if (matpp != NULL)
        (*mat)(K,L) = (*matpp)(k,l);
      else
        (*mat)(K,L) = 0.0;
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementVector(
  LINALG::Matrix<NUMDOF_,1>* vct,
  const LINALG::Matrix<NUMDISP_,1>* vctd,
  const LINALG::Matrix<NUMPRES_,1>* vctp
)
{
  const int d2dp[NUMDISP_] = {0,1,2,  4,5,6,  8,9,10,   12,13,14,   16,17,18,   20,21,22,   24,25,26,   28,29,30  };
  const int p2dp[NUMPRES_] = {      3,      7,       11,         15,         19,         23,         27,        31};
  vct->Clear();
  if (vctd != NULL)
  {
    for (int i=0; i<NUMDISP_; ++i)
    {
      const int I = d2dp[i];
      (*vct)(I,0) = (*vctd)(i,0);
    }
  }
  if (vctp != NULL)
  {
    for (int k=0; k<NUMPRES_; ++k)
    {
      const int K = p2dp[k];
      (*vct)(K,0) = (*vctp)(k,0);
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AssembleVolume(
  Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
  const double& volume  ///< current element volume
  )
{
  double totvol = params.get<double>("volume");
  params.set("volume",totvol+volume);
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
