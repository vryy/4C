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
  const static int Voigt6Row[NUMSTR_] = {0,1,2, 0,1,2};
  const static int Voigt6Col[NUMSTR_] = {0,1,2, 1,2,0};

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
  // 9-Voigt C-index               0 1 2  3 4 5  6 7 8
  const static int Voigt9Row[NUMDFGR_] = {0,1,2, 0,1,2, 0,2,1};
  const static int Voigt9Col[NUMDFGR_] = {0,1,2, 1,2,0, 2,1,0};

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
  const static int Voigt3x3[NUMDFGR_] = {0,3,6, 8,1,4, 5,7,2};

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
  const static int Voigt3x3[NUMDFGR_] = {0,3,5, 3,1,4, 5,4,2};

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
#if 0
        invfderf(ij,kl) += -0.5*(invfmat(j,k)*invfmat(l,i) + invfmat(j,l)*invfmat(k,i));
#else
        invfderf(ij,kl) *= 2.0;
#endif
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

  // VERIFIED

//  cout << endl;
//  cout << endl;
  for (int ij=0; ij<NUMSTR_; ++ij)
  {
//    cout << "[\n";
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMSTR_; ++kl)
    {
      const int k = voigt6row[kl];
      const int l = voigt6col[kl];
      for (int mn=kl; mn<NUMSTR_; ++mn)
      {
        const int m = voigt6row[mn];
        const int n = voigt6col[mn];
        const int klmn = NUMSTR_*kl + mn;
        const int mnkl = NUMSTR_*mn + kl;
        invbvdderb(ij,klmn) = 0.25*(
            ( ibt(i,m)*ibt(n,k) + ibt(i,n)*ibt(m,k) )*ibt(l,j)
            + ibt(i,k)*( ibt(l,m)*ibt(n,j) + ibt(l,n)*ibt(m,j) )
            + ( ibt(i,m)*ibt(n,l) + ibt(i,n)*ibt(m,l) )*ibt(k,j)
            + ibt(i,l)*( ibt(k,m)*ibt(n,j) + ibt(k,n)*ibt(m,j) )
          );
//        cout << ""
//             << "(ct["<<i+1<<","<<m+1<<"]*ct["<<n+1<<","<<k+1<<"]+ct["<<i+1<<","<<n+1<<"]*ct["<<m+1<<","<<k+1<<"])*ct["<<l+1<<","<<j+1<<"]"
//             << "+ct["<<i+1<<","<<k+1<<"]*(ct["<<l+1<<","<<m+1<<"]*ct["<<n+1<<","<<j+1<<"]+ct["<<l+1<<","<<n+1<<"]*ct["<<m+1<<","<<j+1<<"])"
//             << "+(ct["<<i+1<<","<<m+1<<"]*ct["<<n+1<<","<<l+1<<"]+ct["<<i+1<<","<<n+1<<"]*ct["<<m+1<<","<<l+1<<"])*ct["<<k+1<<","<<j+1<<"]"
//             << "+ct["<<i+1<<","<<l+1<<"]*(ct["<<k+1<<","<<m+1<<"]*ct["<<n+1<<","<<j+1<<"]+ct["<<k+1<<","<<n+1<<"]*ct["<<m+1<<","<<j+1<<"])"
//             << "";
        if (ij >= NUMDIM_)  // swap 'i' and 'j' 
        {
#if 0
          invbvdderb(ij,klmn) += 0.25*(
              ( ibt(j,m)*ibt(n,k) + ibt(j,n)*ibt(m,k) )*ibt(l,i)
              + ibt(j,k)*( ibt(l,m)*ibt(n,i) + ibt(l,n)*ibt(m,i) )
              + ( ibt(j,m)*ibt(n,l) + ibt(j,n)*ibt(m,l) )*ibt(k,i)
              + ibt(j,l)*( ibt(k,m)*ibt(n,i) + ibt(k,n)*ibt(m,i) )
            );
#else
          invbvdderb(ij,klmn) *= 2.0;
#endif
//          cout << ""
//               << "+(ct["<<j+1<<","<<m+1<<"]*ct["<<n+1<<","<<k+1<<"]+ct["<<j+1<<","<<n+1<<"]*ct["<<m+1<<","<<k+1<<"])*ct["<<l+1<<","<<i+1<<"]"
//               << "+ct["<<j+1<<","<<k+1<<"]*(ct["<<l+1<<","<<m+1<<"]*ct["<<n+1<<","<<i+1<<"]+ct["<<l+1<<","<<n+1<<"]*ct["<<m+1<<","<<i+1<<"])"
//               << "+(ct["<<j+1<<","<<m+1<<"]*ct["<<n+1<<","<<l+1<<"]+ct["<<j+1<<","<<n+1<<"]*ct["<<m+1<<","<<l+1<<"])*ct["<<k+1<<","<<i+1<<"]"
//               << "+ct["<<j+1<<","<<l+1<<"]*(ct["<<k+1<<","<<m+1<<"]*ct["<<n+1<<","<<i+1<<"]+ct["<<k+1<<","<<n+1<<"]*ct["<<m+1<<","<<i+1<<"])"
//               << "";
        }
        if (mn != kl)
          invbvdderb(ij,mnkl) = invbvdderb(ij,klmn);
//        cout << ",\n";
      }
//      cout << "";
    }
//    cout << "],\n";
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtDiffByItself(
  LINALG::Matrix<NUMSTR_,NUMSTR_>& sqfderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat,
  const VoigtType outvoigt6
  )
{
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  // VERIFIED

#if 0
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
      if ( (outvoigt6 == voigt6_strain) and (ij >= NUMDIM_) )
      {
        sqfderf(ij,kl) += id(j,k)*fmat(l,i) + id(i,l)*fmat(j,k);
//        cout << "+id["<<j+1<<","<<k+1<<"]*St["<<l+1<<","<<i+1<<"]+id["<<i+1<<","<<l+1<<"]*St["<<j+1<<","<<k+1<<"]";
      }
//      cout << ", ";
    }
//    cout << "]," << endl;
  }
#else
  if (outvoigt6 != voigt6_strain)
    dserror("Can only produce row of strain-like type");
  sqfderf(0,0) = 2.0*fmat(0,0);
  sqfderf(1,0) = 0.0;
  sqfderf(2,0) = 0.0;
  sqfderf(3,0) = fmat(1,0)+fmat(0,1);
  sqfderf(4,0) = 0.0;
  sqfderf(5,0) = fmat(2,0)+fmat(0,2);
  sqfderf(0,1) = 0.0;
  sqfderf(1,1) = 2.0*fmat(1,1);
  sqfderf(2,1) = 0.0;
  sqfderf(3,1) = fmat(1,0)+fmat(0,1);
  sqfderf(4,1) = fmat(2,1)+fmat(1,2);
  sqfderf(5,1) = 0.0;
  sqfderf(0,2) = 0.0;
  sqfderf(1,2) = 0.0;
  sqfderf(2,2) = 2.0*fmat(2,2);
  sqfderf(3,2) = 0.0;
  sqfderf(4,2) = fmat(2,1)+fmat(1,2);
  sqfderf(5,2) = fmat(2,0)+fmat(0,2);
  sqfderf(0,3) = fmat(0,1);
  sqfderf(1,3) = fmat(0,1);
  sqfderf(2,3) = 0.0;
  sqfderf(3,3) = fmat(1,1)+fmat(0,0);
  sqfderf(4,3) = fmat(0,2);
  sqfderf(5,3) = fmat(2,1);
  sqfderf(0,4) = 0.0;
  sqfderf(1,4) = fmat(1,2);
  sqfderf(2,4) = fmat(1,2);
  sqfderf(3,4) = fmat(0,2);
  sqfderf(4,4) = fmat(2,2)+fmat(1,1);
  sqfderf(5,4) = fmat(1,0);
  sqfderf(0,5) = fmat(2,0);
  sqfderf(1,5) = 0.0;
  sqfderf(2,5) = fmat(2,0);
  sqfderf(3,5) = fmat(2,1);
  sqfderf(4,5) = fmat(1,0);
  sqfderf(5,5) = fmat(2,2)+fmat(0,0);
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*
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
*/

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtTwiceDiffByItself(
  LINALG::Matrix<NUMSTR_,NUMSTR_*NUMSTR_>& sqfdderf,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& fmat
  )
{
#if 0
  const int* voigt6row = NULL;
  const int* voigt6col = NULL;
  Indices6VoigtTo2Tensor(voigt6row,voigt6col);

  // identity 2-tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> id(true);
  for (int i=0; i<NUMDIM_; ++i) id(i,i) = 1.0;

  // VERIFIED

  // (F^T.F)_{,FF} with F^T=F
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
      for (int mn=0; mn<NUMSTR_; ++mn)
      {
        const int m = voigt6row[mn];
        const int n = voigt6col[mn];
        const int klmn = NUMSTR_*kl + mn;
        double sqfdderf_ijklmn = 0.25*(id(i,k)*id(l,m)*id(j,n)+id(j,l)*id(i,m)*id(k,n)
                                       +id(i,k)*id(l,n)*id(j,m)+id(j,l)*id(i,n)*id(k,m)  // swap 'm' and 'n'
                                       +id(i,l)*id(k,m)*id(j,n)+id(j,k)*id(i,m)*id(l,n)  // swap 'k' and 'l'
                                       +id(i,l)*id(k,n)*id(j,m)+id(j,k)*id(i,n)*id(l,m));  // swap 'm' and 'n' as well as 'k' and 'l'
        if (ij >= NUMDIM_)  // swap 'i' and 'j'
        {
          sqfdderf_ijklmn += 0.25*(id(j,k)*id(l,m)*id(i,n)+id(i,l)*id(j,m)*id(k,n)  // swap 'i' and 'j'
                                   +id(j,k)*id(l,n)*id(i,m)+id(i,l)*id(j,n)*id(k,m)  // swap 'i' and 'j' as well as 'm' and 'n'
                                   +id(j,l)*id(k,m)*id(i,n)+id(i,k)*id(j,m)*id(l,n)  // swap 'i' and 'j' as well as 'k' and 'l'
                                   +id(j,l)*id(k,n)*id(i,m)+id(i,k)*id(j,n)*id(l,m) );  // swap 'i' and 'j' as well as 'm' and 'n' as well as 'k' and 'l'
        }
        sqfdderf(ij,klmn) = sqfdderf_ijklmn;
//        cout << sqfdderf_ijklmn;
//        cout << ", ";
      }
//      cout << "\n";
    }
//    cout << "],\n";       
  }
#else
  sqfdderf.Clear();
  sqfdderf(0,0) = 2.0;
  sqfdderf(0,21) = 0.5;
  sqfdderf(0,35) = 0.5;

  sqfdderf(1,7) = 2.0;
  sqfdderf(1,21) = 0.5;
  sqfdderf(1,28) = 0.5;

  sqfdderf(2,14) = 2.0;
  sqfdderf(2,28) = 0.5;
  sqfdderf(2,35) = 0.5;

  sqfdderf(3,3) = 1.0;
  sqfdderf(3,9) = 1.0;
  sqfdderf(3,18) = 1.0;
  sqfdderf(3,19) = 1.0;
  sqfdderf(3,29) = 0.5;
  sqfdderf(3,34) = 0.5;

  sqfdderf(4,9) = 1.0;
  sqfdderf(4,16) = 1.0;
  sqfdderf(4,23) = 0.5;
  sqfdderf(4,25) = 1.0;
  sqfdderf(4,26) = 1.0;
  sqfdderf(4,33) = 0.5;

  sqfdderf(5,5) = 1.0;
  sqfdderf(5,17) = 1.0;
  sqfdderf(5,22) = 0.5;
  sqfdderf(5,27) = 0.5;
  sqfdderf(5,30) = 1.0;
  sqfdderf(5,32) = 1.0;
#endif

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::SqVector6VoigtTwiceDiffByItself(
  int* isqfdderf,  //[NUMSTR_*6];
  LINALG::Matrix<NUMSTR_,6>& sqfdderf
  )
{
  isqfdderf[NUMSTR_*0+0] = 0;  sqfdderf(0,0) = 2.0;
  isqfdderf[NUMSTR_*0+1] = 21;  sqfdderf(0,1) = 0.5;
  isqfdderf[NUMSTR_*0+2] = 35;  sqfdderf(0,2) = 0.5;
  isqfdderf[NUMSTR_*0+3] = -1;  sqfdderf(0,3) = 0.0;  // dummy
  isqfdderf[NUMSTR_*0+4] = -1;  sqfdderf(0,4) = 0.0;  // dummy
  isqfdderf[NUMSTR_*0+5] = -1;  sqfdderf(0,5) = 0.0;  // dummy

  isqfdderf[NUMSTR_*1+0] = 7;   sqfdderf(1,0) = 2.0;
  isqfdderf[NUMSTR_*1+1] = 21;  sqfdderf(1,1) = 0.5;
  isqfdderf[NUMSTR_*1+2] = 28;  sqfdderf(1,2) = 0.5;
  isqfdderf[NUMSTR_*1+3] = -1;  sqfdderf(1,3) = 0.0;  // dummy
  isqfdderf[NUMSTR_*1+4] = -1;  sqfdderf(1,4) = 0.0;  // dummy
  isqfdderf[NUMSTR_*1+5] = -1;  sqfdderf(1,5) = 0.0;  // dummy

  isqfdderf[NUMSTR_*2+0] = 14;  sqfdderf(2,0) = 2.0;
  isqfdderf[NUMSTR_*2+1] = 28;  sqfdderf(2,1) = 0.5;
  isqfdderf[NUMSTR_*2+2] = 35;  sqfdderf(2,2) = 0.5;
  isqfdderf[NUMSTR_*2+3] = -1;  sqfdderf(2,3) = 0.0;  // dummy
  isqfdderf[NUMSTR_*2+4] = -1;  sqfdderf(2,4) = 0.0;  // dummy
  isqfdderf[NUMSTR_*2+5] = -1;  sqfdderf(2,5) = 0.0;  // dummy

  isqfdderf[NUMSTR_*3+0] = 3;  sqfdderf(3,0) = 1.0;
  isqfdderf[NUMSTR_*3+1] = 9;  sqfdderf(3,1) = 1.0;
  isqfdderf[NUMSTR_*3+2] = 18;  sqfdderf(3,2) = 1.0;
  isqfdderf[NUMSTR_*3+3] = 19;  sqfdderf(3,3) = 1.0;
  isqfdderf[NUMSTR_*3+4] = 29;  sqfdderf(3,4) = 0.5;
  isqfdderf[NUMSTR_*3+5] = 34;  sqfdderf(3,5) = 0.5;

  isqfdderf[NUMSTR_*4+0] = 9;  sqfdderf(4,0) = 1.0;
  isqfdderf[NUMSTR_*4+1] = 16;  sqfdderf(4,1) = 1.0;
  isqfdderf[NUMSTR_*4+2] = 23;  sqfdderf(4,2) = 0.5;
  isqfdderf[NUMSTR_*4+3] = 25;  sqfdderf(4,3) = 1.0;
  isqfdderf[NUMSTR_*4+4] = 26;  sqfdderf(4,4) = 1.0;
  isqfdderf[NUMSTR_*4+5] = 33;  sqfdderf(4,5) = 0.5;

  isqfdderf[NUMSTR_*5+0] = 5;  sqfdderf(5,0) = 1.0;
  isqfdderf[NUMSTR_*5+1] = 17;  sqfdderf(5,1) = 1.0;
  isqfdderf[NUMSTR_*5+2] = 22;  sqfdderf(5,2) = 0.5;
  isqfdderf[NUMSTR_*5+3] = 27;  sqfdderf(5,3) = 0.5;
  isqfdderf[NUMSTR_*5+4] = 30;  sqfdderf(5,4) = 1.0;
  isqfdderf[NUMSTR_*5+5] = 32;  sqfdderf(5,5) = 1.0;

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

  // VERIFIED

//  cout << endl;
  for (int ij=0; ij<NUMSTR_; ++ij)
  {
//    cout << "[";
    const int i = voigt6row[ij];
    const int j = voigt6col[ij];
    for (int kl=0; kl<NUMDFGR_; ++kl)
    {
      const int k = voigt9row[kl];
      const int l = voigt9col[kl];
      if (j == l)
        if (transpose)
        {
          bm(ij,kl) = bt(k,i);
//      cout << "bt["<<k+1<<","<<i+1<<"]";
        }
        else
          bm(ij,kl) = bt(i,k);
      else if ( (ij >= NUMDIM_) and (i == l) )
        if (transpose)
        {
          bm(ij,kl) = bt(k,j);
//      cout << "bt["<<k+1<<","<<j+1<<"]";
        }
        else
          bm(ij,kl) = bt(j,k);
      else
      {
        bm(ij,kl) = 0.0;
//      cout << "0";
      }
//      cout << ", ";
    }
//    cout << "]," << endl;    
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
void DRT::ELEMENTS::So_sh8p8::StretchTensor(
  double* detut,  // determinant of material stretch tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_>* ut,  // material stretch tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_>* invut,  // inverse material stretch tensor
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& ct  // right Cauchy-Green tensor
  )
{
  if ( (ut == NULL) and (invut == NULL) )
    dserror("Senseless call: You do not want to compute anything");

  // pi
//  const double PI = 3.14159265358979323846;

  // set identity tensor
  LINALG::Matrix<NUMDIM_,NUMDIM_> it(true);
  for (int i=0; i<NUMDIM_; ++i) it(i,i) = 1.0;

  // squared right Cauchy-Green deformation tensor
  // C^2 = C . C
  LINALG::Matrix<NUMDIM_,NUMDIM_> c2t;
  c2t.MultiplyTN(ct,ct);

  // invariants of right Cauchy-Green tensor
  // 1st principal invariant: I_C = tr(C)
  const double ci = ct(0,0) + ct(1,1) + ct(2,2);
  // 2nd principal invariant: II_C = 1/2 ( tr(C)^2 - tr(C^2) )
  const double c2i = c2t(0,0) + c2t(1,1) + c2t(2,2);
  const double cii = 0.5*(ci*ci - c2i);
  // 3rd principal invariant: III_C = det(C)
  const double ciii = ct.Determinant();

  // determination of I_U acc. to [1]
  double ui = 0.0;
  {
    // auxiliar variables to get trace of material stretch tensor U
    const double xi = (2.0*ci*ci*ci - 9.0*ci*cii + 27.0*ciii)/27.0;

    double eta = (4.0*cii*cii*cii - ci*ci*cii*cii + 4.0*ci*ci*ci*ciii
                  - 18.0*ci*cii*ciii + 27.0*ciii*ciii)/27.0;
    if (eta < 0.0)
      if (fabs(eta) < EPS6)
        eta = 0.0;
      else
        dserror("Trouble with negative eta=%g", eta);

    // const double zeta = -2.0*ci/3.0
    //                   + pow(xi+sqrt(eta), 1.0/3.0)
    //                   + pow(xi-sqrt(eta), 1.0/3.0);
    double zeta = -2.0*ci/3.0;
    const double xiplussqrteta = 32.0*(xi+sqrt(eta));
    if (xiplussqrteta < 0.0)
      zeta -= pow(fabs(xiplussqrteta), (1.0/3.0));
    else
      zeta += pow(xiplussqrteta, (1.0/3.0));
    const double ximinussqrteta = 32.0*(xi-sqrt(eta));
    if (ximinussqrteta < 0.0)
      zeta -= pow(fabs(ximinussqrteta), (1.0/3.0));
    else
      zeta += pow(ximinussqrteta, (1.0/3.0));

    // invariants of material stretch tensor U
    // 1st invariant: I_U = tr(U)
    if (fabs(zeta+2.0*ci) < EPS12)
    {
      ui = sqrt(ci + 2.0*sqrt(cii));
    }
    else
    {
      const double aux = sqrt(2.0*ci+zeta);
      ui = 0.5*(aux + sqrt(2.0*ci-zeta+16.0*sqrt(ciii)/aux));
    }
  }

  // 2nd and 3rd invariant of material stretch tensor U
  // 2nd invariant: II_U = 1/2 * (I_U^2 - I_C)
  const double uii = 0.5*(ui*ui - ci);  // OR // sqrt(cii + 2.0*sqrt(ciii)*ui);
  // 3rd invariant: III_U = det(U) = sqrt(III_C)
  const double uiii = sqrt(ciii);

  // inverse of material stretch tensor U^{-1}
  // Hoger & Carlson [1] identified
  //     U^{-1} = [ III_U^2*(III_U+I_U*I_C)
  //                + I_U^2*(I_U*III_C + III_U*II_C) ]^{-1}
  //            * [ I_U*(I_U*II_U - III_U)*C^2
  //                - (I_U*II_U - III_U)*(III_U +I_U*I_C)*C
  //                + { II_U*III_U*(III_U+I_U*I_C)
  //                    + I_U^2*(II_U*II_C+III_C) }*I ]
  //            = 1/denom * [ pc2*C^2 + pc*C + pi*I ]
  if (invut != NULL)
  {
    const double denom = uiii*uiii*(uiii + ui*ci)
                       + ui*ui*(ui*ciii + uiii*cii);
    const double pc2 = ui*(ui*uii - uiii);
    const double pc = -(ui*uii - uiii)*(uiii + ui*ci);
    const double pi = uii*uiii*(uiii + ui*ci)
                    + ui*ui*(uii*cii + ciii);
    for (int j=0; j<NUMDIM_; j++)
    {
      for (int i=0; i<NUMDIM_; i++)
      {
        const double invut_ij = (pc2*c2t(i,j) + pc*ct(i,j) + pi*it(i,j))/denom;
        (*invut)(i,j) = invut_ij;
      }
    }
  }

  // material stretch tensor U
  // Hoger & Carlson [1] wrote
  //     U = [ II_U * { II_U * (II_U + I_C) + II_C } + III_C ]^{-1}
  //       * [ -(I_U*II_U - III_U)*C^2
  //           + (I_U*II_U - III_U)*(II_U + I_C)*C
  //           + { I_U*III_U + III_U * ( II_U*(II_U+I_C)+II_C ) }*I ]
  //       = 1/denom * [ pc2*C^2 + pc*C + pi*I ]
  //
  // alternative:
  // U could be calculated based on R later: U = R^T . F
  if (ut != NULL)
  {
    const double denom = uii*(uii*(uii+ci) + cii) + ciii;
    const double pc2 = -(ui*uii - uiii);
    const double pc = (ui*uii - uiii)*(uii + ci);
    const double pi = ui*ciii + uiii*(uii*(uii+ci) + cii);
    for (int j=0; j<NUMDIM_; j++)
    {
      for (int i=0; i<NUMDIM_; i++)
      {
        const double ut_ij = (pc2*c2t(i,j) + pc*ct(i,j) + pi*it(i,j))/denom;
        (*ut)(i,j) = ut_ij;
      }
    }
  }

  // determinat of right stretch tensor
  if (detut != NULL)
  {
    *detut = uiii;
  }

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_sh8p8::SymSpectralDecompJacIter(
  LINALG::Matrix<NUMDIM_,NUMDIM_>& ew,
  LINALG::Matrix<NUMDIM_,NUMDIM_>& ev,
  const LINALG::Matrix<NUMDIM_,NUMDIM_>& at,
  const double itertol,
  const int itermax
  )
{
  // sum of all entries (moduli) in #at
  double asum = 0.0;

  // initialise eigenvalue tensor and eigenvector tensor
#if 0
  asum = at.Norm1();
  ew.Update(at);
  ev.Clear();
  for (int idim=0; idim<NUMDIM_; idim++) ev(idim,idim) = 1.0;
#else
  {
    asum = 0.0;
    for (int jdim=0; jdim<NUMDIM_; jdim++)    
    {
      for (int idim=0; idim<NUMDIM_; idim++)
      {
        asum += fabs(at(idim,jdim));
        ew(idim,jdim) = at(idim,jdim);
        ev(idim,jdim) = 0.0;
      }
      ev(jdim,jdim) = 1.0;
    }
  }
#endif

  // check for trivial problem
  if (asum < EPS12)
  {
    ew.Clear();
    return 0;
  }

  // scale sum of at compenents to achieve relative convergence check
  asum /= (double) (NUMDIM_ * NUMDIM_);

  // reduce ew to diagonal (the eigenvalues)
  double itercnt = 0;  // initialise iteration index <i>
  while (itercnt < itermax)
  {
    double vsum = 0.0;  // sum of all subtriangluar entries
    // loop lower triangle
    for (int jdim=1; jdim<NUMDIM_; jdim++)
    {
      for (int idim=0; idim<jdim; idim++)
      {
        // sum of all triag entries
        vsum += ew(idim,jdim);

        // rotation angle th
        // 2*th = atan(2*evt(idim,jdim)/(ew[idim,idim]-ew[jdim,jdim])
        const double th = 0.5*atan2(2.0*ew(idim,jdim),ew(idim,idim)-ew(jdim,jdim));
        const double sith = sin(th);  // sine of rotation angle
        const double coth = cos(th);  // cosine of rotation angle

        // this defines the rotation matrix,
        // e.g.
        //
        //             [ T_{idim,idim}  0  T_{idim,jdim} ]
        //   T^<i+1> = [             0  1              0 ]
        //             [ T_{jdim,idim}  0  T_{jdim,jdim} ]
        //
        //       [ cos(th)  0  -sin(th) ]
        //     = [       0  1         0 ]
        //       [ sin(th)  0   cos(th) ]
        
        // update eigenvector matrix by right-multiplying with T
        // T is mostly 0 thus it is more efficient to do explicitly
        //    ev^<i+1> = ev^<i> . T^<i+1>
        for (int kdim=0; kdim<NUMDIM_; kdim++)
        {
          const double evki = ev(kdim,idim);
          ev(kdim,idim) = coth*evki + sith*ev(kdim,jdim);
          ev(kdim,jdim) = -sith*evki + coth*ev(kdim,jdim);
        }

        // update eigenvalue tensor by right-multiplying with T and
        // left-multiplying with transposed T
        //    ew^<i+1> = transposed(T^<i+1>) . ew^<i> . T^<i+1>
        // modify "idim" and "jdim" columns
        for (int kdim=0; kdim<NUMDIM_; kdim++)
        {
          const double ewki = ew(kdim,idim);
          ew(kdim,idim) = coth*ewki + sith*ew(kdim,jdim);
          ew(kdim,jdim) = -sith*ewki + coth*ew(kdim,jdim);
        }
        // modify diagonal terms
        ew(idim,idim) = coth*ew(idim,idim) + sith*ew(jdim,idim);
        ew(jdim,jdim) = -sith*ew(idim,jdim) + coth*ew(jdim,jdim);
        ew(idim,jdim) = 0.0;
        // make symmetric
        for (int kdim=0; kdim<NUMDIM_; kdim++)
        {
          ew(idim,kdim) = ew(kdim,idim);
          ew(jdim,kdim) = ew(kdim,jdim);
        }
      }
    }

    // check convergence
    if (fabs(vsum)/asum < itertol)
    {
      break;
    }
    // increment iteration index
    itercnt += 1;
  }

  // check if iteration loop diverged
  int err = 1;
  if (itercnt == itermax)
  {
    err = 1;  // failed
    dserror("Divergent spectral decomposition (Jacobi's iterative method)!");
  }
  else
  {
    err = 0;  // passed
  }

  return err;
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

  // k_dd
  for (int j=0; j<NUMDISP_; ++j) {
    const int J = d2dp[j];
    for (int i=0; i<NUMDISP_; ++i) {
      const int I = d2dp[i];
      if (matdd != NULL)
        (*mat)(I,J) = (*matdd)(i,j);
      else
        (*mat)(I,J) = 0.0;
    }
  }

  // k_dp
  if (matdp != NULL) {
    for (int l=0; l<NUMPRES_; ++l) {
      const int L = p2dp[l];
      for (int i=0; i<NUMDISP_; ++i) {
        const int I = d2dp[i];
        (*mat)(I,L) = (*matdp)(i,l);
      }
    }
  }
  else {
    for (int l=0; l<NUMPRES_; ++l) {
      const int L = p2dp[l];
      for (int i=0; i<NUMDISP_; ++i) {
        const int I = d2dp[i];
        (*mat)(I,L) = 0.0;
      }
    }
  }

  // k_pd
  if (matpd != NULL) {
    for (int j=0; j<NUMDISP_; ++j) {
      const int J = d2dp[j];
      for (int k=0; k<NUMPRES_; ++k) {
        const int K = p2dp[k];
        (*mat)(K,J) = (*matpd)(k,j);
      }
    }
  }
  else if (matdp != NULL) {
    for (int j=0; j<NUMDISP_; ++j) {
      const int J = d2dp[j];
      for (int k=0; k<NUMPRES_; ++k) {
        const int K = p2dp[k];
        (*mat)(K,J) = (*matdp)(j,k);
      }
    }
  }
  else {
    for (int j=0; j<NUMDISP_; ++j) {
      const int J = d2dp[j];
      for (int k=0; k<NUMPRES_; ++k) {
        const int K = p2dp[k];
        (*mat)(K,J) = 0.0;
      }
    }
  }

  // k_pp
  if (matpp != NULL) {
    for (int l=0; l<NUMPRES_; ++l) {
      const int L = p2dp[l];
      for (int k=0; k<NUMPRES_; ++k) {
        const int K = p2dp[k];
        (*mat)(K,L) = (*matpp)(k,l);
      }
    }
  }
  else {
    for (int l=0; l<NUMPRES_; ++l) {
      const int L = p2dp[l];
      for (int k=0; k<NUMPRES_; ++k) {
        const int K = p2dp[k];
        (*mat)(K,L) = 0.0;
      }
    }
  }

  // done
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
  
  // r_d
  if (vctd != NULL) {
    for (int i=0; i<NUMDISP_; ++i) {
      const int I = d2dp[i];
      (*vct)(I,0) = (*vctd)(i,0);
    }
  }

  // r_p
  if (vctp != NULL) {
    for (int k=0; k<NUMPRES_; ++k) {
      const int K = p2dp[k];
      (*vct)(K,0) = (*vctp)(k,0);
    }
  }

  // What shall we do with a drunken sailor?
  return; 
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::AssembleVolume(
  Teuchos::ParameterList& params,  ///< parameter list for in 'n' out
  const double& volume  ///< current element volume
  )
{
  const double totvol = params.get<double>("volume");
  params.set("volume",totvol+volume);
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
