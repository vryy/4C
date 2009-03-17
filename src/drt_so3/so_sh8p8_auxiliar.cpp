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
  const int Voigt6Row[NUMSTR_SOSH8P8] = {0,1,2, 0,1,2};
  const int Voigt6Col[NUMSTR_SOSH8P8] = {0,1,2, 1,2,0};

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
  const int Voigt9Row[NUMDFGR_SOSH8P8] = {0,1,2, 0,1,2, 0,2,1};
  const int Voigt9Col[NUMDFGR_SOSH8P8] = {0,1,2, 1,2,0, 2,1,0};

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
  const int Voigt3x3[NUMDFGR_SOSH8P8] = {0,3,6, 8,1,4, 5,7,2};

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
  const int Voigt3x3[NUMDFGR_SOSH8P8] = {0,3,5, 3,1,4, 5,4,2};

  voigt3x3 = &(Voigt3x3[0]);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector9Voigt(
  LINALG::Matrix<NUMDFGR_SOSH8P8,1>& fvct,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat,
  const bool transpose
  )
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col,transpose);
    
  for (int I=0; I<NUMDFGR_SOSH8P8; ++I)
  {
    const int i = voigt9row[I];
    const int j = voigt9col[I];
    fvct(I,0) = fmat(i,j);  // F_ij
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToVector6Voigt(
  LINALG::Matrix<NUMSTR_SOSH8P8,1>& bvct,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bmat,
  const VoigtType outvoigt6
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);
    
  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigt6row[I];
    const int j = voigt6col[I];
    if (I < NUMDIM_SOSH8P8)
      bvct(I) = bmat(i,j);  // B_ij
    else
      if (outvoigt6 == voigt6_strain)
        bvct(I) = bmat(i,j) + bmat(j,i);  // B_ij+B_ji
      else
        bvct(I) = bmat(i,j);  // B_ij
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Vector6VoigtToMatrix2Tensor(
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bmat,
  const LINALG::Matrix<NUMSTR_SOSH8P8,1>& bvct,
  const VoigtType invoigt6
  )
{
  const int* voigt3x3sym = NULL;
  Indices2TensorTo6Voigt(voigt3x3sym);  // access is via (i,j) -> 3*i+j  

  for (int i=0; i<NUMDIM_SOSH8P8; ++ i)
  {
    for (int j=0; j<NUMDIM_SOSH8P8; ++j)
    {
      const int I = voigt3x3sym[NUMDIM_SOSH8P8*i+j];
      if (i==j)
        bmat(i,j) = bvct(I);
      else
        if (invoigt6 == voigt6_strain)
          bmat(i,j) = 0.5*bvct(I);
        else
          bmat(i,j) = bvct(I);
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector9VoigtDiffByItself(
  LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDFGR_SOSH8P8>& invfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invfmat,
  const bool transpose
  )
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col);

  // VERIFIED

  for (int I=0; I<NUMDFGR_SOSH8P8; ++I)
  {
    const int i = voigt9row[I];
    const int j = voigt9col[I];
    for (int K=0; K<NUMDFGR_SOSH8P8; ++K)
    {
      const int k = voigt9row[K];
      const int l = voigt9col[K];
      if (transpose)
        invfderf(I,K) = -invfmat(j,k)*invfmat(l,i);
      else
        invfderf(I,K) = -invfmat(i,k)*invfmat(l,j);
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::InvVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& invfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& invfmat
  )
{
  const int voigt6row[NUMSTR_SOSH8P8] = {0,1,2, 0,1,2};
  const int voigt6col[NUMSTR_SOSH8P8] = {0,1,2, 1,2,0};

  // VERIFIED

//  cout << endl;
  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
//    cout << "[";
    const int i = voigt6row[I];
    const int j = voigt6col[I];
    for (int K=0; K<NUMSTR_SOSH8P8; ++K)
    {
      const int k = voigt6row[K];
      const int l = voigt6col[K];
      invfderf(I,K) = -0.5*(invfmat(i,k)*invfmat(l,j) + invfmat(i,l)*invfmat(k,j));
//     cout << "ct["<<i+1<<","<<k+1<<"]*ct["<<l+1<<","<<j+1<<"]+ct["<<i+1<<","<<l+1<<"]*ct["<<k+1<<","<<j+1<<"]";
      if (I >= NUMDIM_SOSH8P8)
      {
        invfderf(I,K) += -0.5*(invfmat(j,k)*invfmat(l,i) + invfmat(j,l)*invfmat(k,i));
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
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8*NUMSTR_SOSH8P8>& invbvdderb,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigt6row[I];
    const int j = voigt6col[I];
    for (int K=0; K<NUMSTR_SOSH8P8; ++K)
    {
      const int k = voigt6row[K];
      const int l = voigt6col[K];
      for (int M=0; M<NUMSTR_SOSH8P8; ++M)
      {
        const int m = voigt6row[M];
        const int n = voigt6col[M];
        const int N = NUMSTR_SOSH8P8*K + M;
        invbvdderb(I,N) = 0.25*(
          ( bt(i,m)*bt(n,k) + bt(i,n)*bt(m,k) )*bt(l,j)
          + bt(i,k)*( bt(l,m)*bt(n,j) + bt(l,n)*bt(m,j) )
          + ( bt(i,m)*bt(n,l) + bt(i,n)*bt(m,l) )*bt(k,j)
          + bt(i,l)*( bt(k,m)*bt(n,j) + bt(k,n)*bt(m,j) )
          );
        if (I >= NUMDIM_SOSH8P8)  // swap i 'n' j 
          invbvdderb(I,N) += 0.25*(
            ( bt(j,m)*bt(n,k) + bt(j,n)*bt(m,k) )*bt(l,i)
            + bt(j,k)*( bt(l,m)*bt(n,i) + bt(l,n)*bt(m,i) )
            + ( bt(j,m)*bt(n,l) + bt(j,n)*bt(m,l) )*bt(k,i)
            + bt(j,l)*( bt(k,m)*bt(n,i) + bt(k,n)*bt(m,i) )
            );
      }
    }
  }

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& sqfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  // VERIFIED

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> id(true);
  for (int i=0; i<NUMDIM_SOSH8P8; ++i) id(i,i) = 1.0;

  // (F.F)_{,F} with F^T=F
//  cout << endl;
  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
//    cout << "[";
    const int i = voigt6row[I];
    const int j = voigt6col[I];
    for (int K=0; K<NUMSTR_SOSH8P8; ++K)
    {
      const int k = voigt6row[K];
      const int l = voigt6col[K];
      sqfderf(I,K) = id(i,k)*fmat(l,j) + id(j,l)*fmat(i,k);
//      cout << "id["<<i+1<<","<<k+1<<"]*St["<<l+1<<","<<j+1<<"]+id["<<j+1<<","<<l+1<<"]*St["<<i+1<<","<<k+1<<"]";
      if (I >= NUMDIM_SOSH8P8)
      {
        sqfderf(I,K) += id(j,k)*fmat(l,i) + id(i,l)*fmat(j,k);
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
  LINALG::Matrix<NUMDFGR_SOSH8P8,NUMDFGR_SOSH8P8>& sqfderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat,
  const bool transpose
  )
{
  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col);

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> id(true);
  for (int i=0; i<NUMDIM_SOSH8P8; ++i) id(i,i) = 1.0;

  // (F^T.F)_{,F}
//  cout << endl;
  for (int I=0; I<NUMDFGR_SOSH8P8; ++I)
  {
//    cout << "[";
    const int i = voigt9row[I];  // A
    const int j = voigt9col[I];  // B
    for (int K=0; K<NUMDFGR_SOSH8P8; ++K)
    {
      const int k = voigt9row[K];  // C
      const int l = voigt9col[K];  // D
//      cout << "i=" << i << ", j=" << j << ", k=" << k << ", l=" << l << endl;
      if (transpose)
        sqfderf(I,K) = fmat(i,k)*id(j,l) + id(i,l)*fmat(j,k);
//        sqfderf(I,K) = id(i,k)*fmat(j,l) + id(j,l)*fmat(k,i);
      else
        sqfderf(I,K) = fmat(k,i)*id(j,l) + id(i,l)*fmat(k,j);
//        sqfderf(I,K) = id(i,k)*fmat(l,j) + id(j,l)*fmat(i,k);
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
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8*NUMSTR_SOSH8P8>& sqfdderf,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& fmat
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8> id(true);
  for (int i=0; i<NUMDIM_SOSH8P8; ++i) id(i,i) = 1.0;

  // (F^T.F)_{,FF} with F^T=F
  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigt6row[I];
    const int j = voigt6col[I];
    for (int K=0; K<NUMSTR_SOSH8P8; ++K)
    {
      const int k = voigt6row[K];
      const int l = voigt6col[K];
      for (int M=0; M<NUMSTR_SOSH8P8; ++M)
      {
        const int m = voigt6row[M];
        const int n = voigt6col[M];
        const int N = NUMSTR_SOSH8P8*K + M;
        sqfdderf(I,N) = id(k,m)*id(i,n)*id(j,l) + id(i,l)*id(k,m)*id(j,n);
        if (I >= NUMDIM_SOSH8P8)
          sqfdderf(I,K) += id(k,m)*id(j,n)*id(i,l) + id(j,l)*id(k,m)*id(i,n);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToMatrix6x9Voigt(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMDFGR_SOSH8P8>& bm,
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt,
  const bool transpose
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  const int* voigt9row = NULL;
  const int* voigt9col = NULL;
  Indices9VoigtTo2Tensor(voigt9row,voigt9col);

  for (int I=0; I<NUMSTR_SOSH8P8; ++I)
  {
    const int i = voigt6row[I];
    const int j = voigt6col[I];
    for (int K=0; K<NUMDFGR_SOSH8P8; ++K)
    {
      const int k = voigt9row[K];
      const int l = voigt9col[K];
      if (j == l)
        bm(I,K) = (transpose) ? bt(k,i) : bt(i,k);
      else if ( (I >= NUMDIM_SOSH8P8) and (i == l) )
        bm(I,K) = (transpose) ? bt(k,j) : bt(j,k);
      else
        bm(I,K) = 0.0;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::Matrix2TensorToLeftRightProductMatrix6x6Voigt(
  LINALG::Matrix<NUMSTR_SOSH8P8,NUMSTR_SOSH8P8>& bm,  ///< (out) 6x6 Voigt matrix
  const LINALG::Matrix<NUMDIM_SOSH8P8,NUMDIM_SOSH8P8>& bt,  ///< (in) 3x3 matrix of 2-tensor
  const bool transpose, ///< 3x3 input matrix is transposed
  const VoigtType outvoigt6,  ///< 6-Voigt vector layout on rows of 6x6 matrix
  const VoigtType invoigt6  ///< 6-Voigt vector layout on columns of 6x6 matrix
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  for (int ab=0; ab<NUMSTR_SOSH8P8; ++ab)
  {
    const int a = voigt6row[ab];
    const int b = voigt6col[ab];
    for (int AB=0; AB<NUMSTR_SOSH8P8; ++AB)
    {
      const int A = voigt6row[AB];
      const int B = voigt6col[AB];
      if (transpose)
      {
        bm(AB,ab) = bt(A,a)*bt(B,b);
        if (ab >= NUMSTR_SOSH8P8) bm(AB,ab) += bt(A,b)*bt(B,a);
      }
      else
      {
        bm(AB,ab) = bt(a,A)*bt(b,B);
        if (ab >= NUMSTR_SOSH8P8) bm(AB,ab) += bt(b,A)*bt(a,B);
      }
      if ( (invoigt6 == voigt6_stress) and (ab >= NUMSTR_SOSH8P8) )
        bm(AB,ab) *= 2.0;
      if ( (outvoigt6 == voigt6_stress) and (AB >= NUMSTR_SOSH8P8) )
        bm(AB,ab) *= 0.5;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::ExtractDispAndPres(
  std::vector<double>& mystat,
  LINALG::Matrix<NUMDISP_SOSH8P8,1>& mydisp,
  LINALG::Matrix<NUMPRES_SOSH8P8,1>& mypres
  )
{
  for (int inod=0; inod<NUMNOD_SOSH8P8; ++inod)
  {
    for (int idis=0; idis<NODDISP_SOSH8P8; ++idis)
      mydisp(idis+(inod*NODDISP_SOSH8P8),0) = mystat[idis+(inod*NODDOF_SOSH8P8)];
    mypres(inod,0) = mystat[3+(inod*NODDOF_SOSH8P8)];
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::BuildElementMatrix(
  LINALG::Matrix<NUMDOF_SOSH8P8,NUMDOF_SOSH8P8>* mat,
  const LINALG::Matrix<NUMDISP_SOSH8P8,NUMDISP_SOSH8P8>* matdd,
  const LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* matdp,
  const LINALG::Matrix<NUMDISP_SOSH8P8,NUMPRES_SOSH8P8>* matpd,
  const LINALG::Matrix<NUMPRES_SOSH8P8,NUMPRES_SOSH8P8>* matpp
)
{
  const int d2dp[NUMDISP_SOSH8P8] = {0,1,2,  4,5,6,  8,9,10,   12,13,14,   16,17,18,   20,21,22,   24,25,26,   28,29,30  };
  const int p2dp[NUMPRES_SOSH8P8] = {      3,      7,       11,         15,         19,         23,         27,        31};
  for (int i=0; i<NUMDISP_SOSH8P8; ++i)
  {
    const int I = d2dp[i];
    for (int j=0; j<NUMDISP_SOSH8P8; ++j)
    {
      const int J = d2dp[j];
      if (matdd != NULL)
        (*mat)(I,J) = (*matdd)(i,j);
      else
        (*mat)(I,J) = 0.0;
    }

    for (int l=0; l<NUMPRES_SOSH8P8; ++l)
    {
      const int L = p2dp[l];
      if (matdp != NULL)
        (*mat)(I,L) = (*matdp)(i,l);
      else
        (*mat)(I,L) = 0.0;
    }
  }
  for (int k=0; k<NUMPRES_SOSH8P8; ++k)
  {
    const int K = p2dp[k];
    for (int j=0; j<NUMDISP_SOSH8P8; ++j)
    {
      const int J = d2dp[j];
      if (matpd != NULL)
        (*mat)(K,J) = (*matpd)(k,j);
      else if (matdp != NULL)
        (*mat)(K,J) = (*matdp)(j,k);
      else
        (*mat)(K,J) = 0.0;
    }
    for (int l=0; l<NUMPRES_SOSH8P8; ++l)
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
  LINALG::Matrix<NUMDOF_SOSH8P8,1>* vct,
  const LINALG::Matrix<NUMDISP_SOSH8P8,1>* vctd,
  const LINALG::Matrix<NUMPRES_SOSH8P8,1>* vctp
)
{
  const int d2dp[NUMDISP_SOSH8P8] = {0,1,2,  4,5,6,  8,9,10,   12,13,14,   16,17,18,   20,21,22,   24,25,26,   28,29,30  };
  const int p2dp[NUMPRES_SOSH8P8] = {      3,      7,       11,         15,         19,         23,         27,        31};
  vct->Clear();
  if (vctd != NULL)
  {
    for (int i=0; i<NUMDISP_SOSH8P8; ++i)
    {
      const int I = d2dp[i];
      (*vct)(I,0) = (*vctd)(i,0);
    }
  }
  if (vctp != NULL)
  {
    for (int k=0; k<NUMPRES_SOSH8P8; ++k)
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
