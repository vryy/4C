/*!----------------------------------------------------------------------
\file inversedesign_evaluate.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)

#include "inversedesign.H"
#include "so_hex8.H"
#include "so_weg6.H"
#include "../drt_mat/material.H"
#include "so_tet4.H"
#include "../drt_lib/drt_dserror.H"
#include "Epetra_SerialDenseSolver.h"


/*----------------------------------------------------------------------*
 |  Lambda tensor                                             (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildLambda(LINALG::SerialDenseMatrix& L, 
                                           const LINALG::SerialDenseMatrix& F) const

{
#if 0 // this is for theory only
  double Lambda4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int p=0; p<3; ++p)
      for (int q=0; q<3; ++q)
        for (int m=0; m<3; ++m)
          Lambda4[k][m][p][q] = F(k,p)*F(q,m);
#endif


#if 1
  L(0,0) = F(0,0)*F(0,0);
  L(0,1) = F(0,1)*F(0,0);
  L(0,2) = F(0,2)*F(0,0);
  L(0,3) = F(0,0)*F(1,0);
  L(0,4) = F(0,1)*F(1,0);
  L(0,5) = F(0,2)*F(1,0);
  L(0,6) = F(0,0)*F(2,0);
  L(0,7) = F(0,1)*F(2,0);
  L(0,8) = F(0,2)*F(2,0);
  
  L(1,0) = F(1,0)*F(0,0);
  L(1,1) = F(1,1)*F(0,0);
  L(1,2) = F(1,2)*F(0,0);
  L(1,3) = F(1,0)*F(1,0);
  L(1,4) = F(1,1)*F(1,0);
  L(1,5) = F(1,2)*F(1,0);
  L(1,6) = F(1,0)*F(2,0);
  L(1,7) = F(1,1)*F(2,0);
  L(1,8) = F(1,2)*F(2,0);
  
  L(2,0) = F(2,0)*F(0,0);
  L(2,1) = F(2,1)*F(0,0);
  L(2,2) = F(2,2)*F(0,0);
  L(2,3) = F(2,0)*F(1,0);
  L(2,4) = F(2,1)*F(1,0);
  L(2,5) = F(2,2)*F(1,0);
  L(2,6) = F(2,0)*F(2,0);
  L(2,7) = F(2,1)*F(2,0);
  L(2,8) = F(2,2)*F(2,0);
  
  L(3,0) = F(0,0)*F(0,1);
  L(3,1) = F(0,1)*F(0,1);
  L(3,2) = F(0,2)*F(0,1);
  L(3,3) = F(0,0)*F(1,1);
  L(3,4) = F(0,1)*F(1,1);
  L(3,5) = F(0,2)*F(1,1);
  L(3,6) = F(0,0)*F(2,1);
  L(3,7) = F(0,1)*F(2,1);
  L(3,8) = F(0,2)*F(2,1);
  
  L(4,0) = F(1,0)*F(0,1);
  L(4,1) = F(1,1)*F(0,1);
  L(4,2) = F(1,2)*F(0,1);
  L(4,3) = F(1,0)*F(1,1);
  L(4,4) = F(1,1)*F(1,1);
  L(4,5) = F(1,2)*F(1,1);
  L(4,6) = F(1,0)*F(2,1);
  L(4,7) = F(1,1)*F(2,1);
  L(4,8) = F(1,2)*F(2,1);
  
  L(5,0) = F(2,0)*F(0,1);
  L(5,1) = F(2,1)*F(0,1);
  L(5,2) = F(2,2)*F(0,1);
  L(5,3) = F(2,0)*F(1,1);
  L(5,4) = F(2,1)*F(1,1);
  L(5,5) = F(2,2)*F(1,1);
  L(5,6) = F(2,0)*F(2,1);
  L(5,7) = F(2,1)*F(2,1);
  L(5,8) = F(2,2)*F(2,1);
  
  L(6,0) = F(0,0)*F(0,2);
  L(6,1) = F(0,1)*F(0,2);
  L(6,2) = F(0,2)*F(0,2);
  L(6,3) = F(0,0)*F(1,2);
  L(6,4) = F(0,1)*F(1,2);
  L(6,5) = F(0,2)*F(1,2);
  L(6,6) = F(0,0)*F(2,2);
  L(6,7) = F(0,1)*F(2,2);
  L(6,8) = F(0,2)*F(2,2);
  
  L(7,0) = F(1,0)*F(0,2);
  L(7,1) = F(1,1)*F(0,2);
  L(7,2) = F(1,2)*F(0,2);
  L(7,3) = F(1,0)*F(1,2);
  L(7,4) = F(1,1)*F(1,2);
  L(7,5) = F(1,2)*F(1,2);
  L(7,6) = F(1,0)*F(2,2);
  L(7,7) = F(1,1)*F(2,2);
  L(7,8) = F(1,2)*F(2,2);
  
  L(8,0) = F(2,0)*F(0,2);
  L(8,1) = F(2,1)*F(0,2);
  L(8,2) = F(2,2)*F(0,2);
  L(8,3) = F(2,0)*F(1,2);
  L(8,4) = F(2,1)*F(1,2);
  L(8,5) = F(2,2)*F(1,2);
  L(8,6) = F(2,0)*F(2,2);
  L(8,7) = F(2,1)*F(2,2);
  L(8,8) = F(2,2)*F(2,2);
#else
  //               k  m  p  q
  L(0,0) = Lambda4[0][0][0][0];
  L(0,1) = Lambda4[0][0][1][0];
  L(0,2) = Lambda4[0][0][2][0];
  L(0,3) = Lambda4[0][0][0][1];
  L(0,4) = Lambda4[0][0][1][1];
  L(0,5) = Lambda4[0][0][2][1];
  L(0,6) = Lambda4[0][0][0][2];
  L(0,7) = Lambda4[0][0][1][2];
  L(0,8) = Lambda4[0][0][2][2];
  
  L(1,0) = Lambda4[1][0][0][0];
  L(1,1) = Lambda4[1][0][1][0];
  L(1,2) = Lambda4[1][0][2][0];
  L(1,3) = Lambda4[1][0][0][1];
  L(1,4) = Lambda4[1][0][1][1];
  L(1,5) = Lambda4[1][0][2][1];
  L(1,6) = Lambda4[1][0][0][2];
  L(1,7) = Lambda4[1][0][1][2];
  L(1,8) = Lambda4[1][0][2][2];
  
  L(2,0) = Lambda4[2][0][0][0];
  L(2,1) = Lambda4[2][0][1][0];
  L(2,2) = Lambda4[2][0][2][0];
  L(2,3) = Lambda4[2][0][0][1];
  L(2,4) = Lambda4[2][0][1][1];
  L(2,5) = Lambda4[2][0][2][1];
  L(2,6) = Lambda4[2][0][0][2];
  L(2,7) = Lambda4[2][0][1][2];
  L(2,8) = Lambda4[2][0][2][2];
  
  L(3,0) = Lambda4[0][1][0][0];
  L(3,1) = Lambda4[0][1][1][0];
  L(3,2) = Lambda4[0][1][2][0];
  L(3,3) = Lambda4[0][1][0][1];
  L(3,4) = Lambda4[0][1][1][1];
  L(3,5) = Lambda4[0][1][2][1];
  L(3,6) = Lambda4[0][1][0][2];
  L(3,7) = Lambda4[0][1][1][2];
  L(3,8) = Lambda4[0][1][2][2];
  
  L(4,0) = Lambda4[1][1][0][0];
  L(4,1) = Lambda4[1][1][1][0];
  L(4,2) = Lambda4[1][1][2][0];
  L(4,3) = Lambda4[1][1][0][1];
  L(4,4) = Lambda4[1][1][1][1];
  L(4,5) = Lambda4[1][1][2][1];
  L(4,6) = Lambda4[1][1][0][2];
  L(4,7) = Lambda4[1][1][1][2];
  L(4,8) = Lambda4[1][1][2][2];
  
  L(5,0) = Lambda4[2][1][0][0];
  L(5,1) = Lambda4[2][1][1][0];
  L(5,2) = Lambda4[2][1][2][0];
  L(5,3) = Lambda4[2][1][0][1];
  L(5,4) = Lambda4[2][1][1][1];
  L(5,5) = Lambda4[2][1][2][1];
  L(5,6) = Lambda4[2][1][0][2];
  L(5,7) = Lambda4[2][1][1][2];
  L(5,8) = Lambda4[2][1][2][2];
  
  L(6,0) = Lambda4[0][2][0][0];
  L(6,1) = Lambda4[0][2][1][0];
  L(6,2) = Lambda4[0][2][2][0];
  L(6,3) = Lambda4[0][2][0][1];
  L(6,4) = Lambda4[0][2][1][1];
  L(6,5) = Lambda4[0][2][2][1];
  L(6,6) = Lambda4[0][2][0][2];
  L(6,7) = Lambda4[0][2][1][2];
  L(6,8) = Lambda4[0][2][2][2];
  
  L(7,0) = Lambda4[1][2][0][0];
  L(7,1) = Lambda4[1][2][1][0];
  L(7,2) = Lambda4[1][2][2][0];
  L(7,3) = Lambda4[1][2][0][1];
  L(7,4) = Lambda4[1][2][1][1];
  L(7,5) = Lambda4[1][2][2][1];
  L(7,6) = Lambda4[1][2][0][2];
  L(7,7) = Lambda4[1][2][1][2];
  L(7,8) = Lambda4[1][2][2][2];
  
  L(8,0) = Lambda4[2][2][0][0];
  L(8,1) = Lambda4[2][2][1][0];
  L(8,2) = Lambda4[2][2][2][0];
  L(8,3) = Lambda4[2][2][0][1];
  L(8,4) = Lambda4[2][2][1][1];
  L(8,5) = Lambda4[2][2][2][1];
  L(8,6) = Lambda4[2][2][0][2];
  L(8,7) = Lambda4[2][2][1][2];
  L(8,8) = Lambda4[2][2][2][2];
#endif
  
  return; 
} // DRT::ELEMENTS::InvDesign::BuildLambda

/*----------------------------------------------------------------------*
 |  IF tensor                                                 (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildIF(LINALG::SerialDenseMatrix& IF, 
                                       const LINALG::SerialDenseMatrix& F) const

{
#if 0 // this is for theory reference only
  double IF4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int l=0; l<3; ++l)
      for (int m=0; m<3; ++m)
        for (int n=0; n<3; ++n)
          IF4[k][l][m][n] = 0.5 * ( F(k,m)*F(l,n) + F(k,n)*F(l,m) );
#endif
  
#if 1
  IF(0,0) = ( F(0,0)*F(0,0) + F(0,0)*F(0,0) );
  IF(0,1) = ( F(0,1)*F(0,1) + F(0,1)*F(0,1) );
  IF(0,2) = ( F(0,2)*F(0,2) + F(0,2)*F(0,2) );
  IF(0,3) = ( F(0,0)*F(0,1) + F(0,1)*F(0,0) ) * 2.0;
  IF(0,4) = ( F(0,1)*F(0,2) + F(0,2)*F(0,1) ) * 2.0;
  IF(0,5) = ( F(0,2)*F(0,0) + F(0,0)*F(0,2) ) * 2.0;

  IF(1,0) = ( F(1,0)*F(1,0) + F(1,0)*F(1,0) );
  IF(1,1) = ( F(1,1)*F(1,1) + F(1,1)*F(1,1) );
  IF(1,2) = ( F(1,2)*F(1,2) + F(1,2)*F(1,2) );
  IF(1,3) = ( F(1,0)*F(1,1) + F(1,1)*F(1,0) ) * 2.0;
  IF(1,4) = ( F(1,1)*F(1,2) + F(1,2)*F(1,1) ) * 2.0;
  IF(1,5) = ( F(1,2)*F(1,0) + F(1,0)*F(1,2) ) * 2.0;

  IF(2,0) = ( F(2,0)*F(2,0) + F(2,0)*F(2,0) );
  IF(2,1) = ( F(2,1)*F(2,1) + F(2,1)*F(2,1) );
  IF(2,2) = ( F(2,2)*F(2,2) + F(2,2)*F(2,2) );
  IF(2,3) = ( F(2,0)*F(2,1) + F(2,1)*F(2,0) ) * 2.0;
  IF(2,4) = ( F(2,1)*F(2,2) + F(2,2)*F(2,1) ) * 2.0;
  IF(2,5) = ( F(2,2)*F(2,0) + F(2,0)*F(2,2) ) * 2.0;

  IF(3,0) = ( F(0,0)*F(1,0) + F(0,0)*F(1,0) );
  IF(3,1) = ( F(0,1)*F(1,1) + F(0,1)*F(1,1) );
  IF(3,2) = ( F(0,2)*F(1,2) + F(0,2)*F(1,2) );
  IF(3,3) = ( F(0,0)*F(1,1) + F(0,1)*F(1,0) ) * 2.0;
  IF(3,4) = ( F(0,1)*F(1,2) + F(0,2)*F(1,1) ) * 2.0;
  IF(3,5) = ( F(0,2)*F(1,0) + F(0,0)*F(1,2) ) * 2.0;

  IF(4,0) = ( F(1,0)*F(2,0) + F(1,0)*F(2,0) );
  IF(4,1) = ( F(1,1)*F(2,1) + F(1,1)*F(2,1) );
  IF(4,2) = ( F(1,2)*F(2,2) + F(1,2)*F(2,2) );
  IF(4,3) = ( F(1,0)*F(2,1) + F(1,1)*F(2,0) ) * 2.0;
  IF(4,4) = ( F(1,1)*F(2,2) + F(1,2)*F(2,1) ) * 2.0;
  IF(4,5) = ( F(1,2)*F(2,0) + F(1,0)*F(2,2) ) * 2.0;

  IF(5,0) = ( F(2,0)*F(0,0) + F(2,0)*F(0,0) );
  IF(5,1) = ( F(2,1)*F(0,1) + F(2,1)*F(0,1) );
  IF(5,2) = ( F(2,2)*F(0,2) + F(2,2)*F(0,2) );
  IF(5,3) = ( F(2,0)*F(0,1) + F(2,1)*F(0,0) ) * 2.0;
  IF(5,4) = ( F(2,1)*F(0,2) + F(2,2)*F(0,1) ) * 2.0;
  IF(5,5) = ( F(2,2)*F(0,0) + F(2,0)*F(0,2) ) * 2.0;

  IF.Scale(0.5);
  
#else
  //            k  l  m  n
  IF(0,0) = IF4[0][0][0][0];
  IF(0,1) = IF4[0][0][1][1];
  IF(0,2) = IF4[0][0][2][2];
  IF(0,3) = IF4[0][0][0][1] * 2.0;
  IF(0,4) = IF4[0][0][1][2] * 2.0;
  IF(0,5) = IF4[0][0][2][0] * 2.0;

  IF(1,0) = IF4[1][1][0][0];
  IF(1,1) = IF4[1][1][1][1];
  IF(1,2) = IF4[1][1][2][2];
  IF(1,3) = IF4[1][1][0][1] * 2.0;
  IF(1,4) = IF4[1][1][1][2] * 2.0;
  IF(1,5) = IF4[1][1][2][0] * 2.0;

  IF(2,0) = IF4[2][2][0][0];
  IF(2,1) = IF4[2][2][1][1];
  IF(2,2) = IF4[2][2][2][2];
  IF(2,3) = IF4[2][2][0][1] * 2.0;
  IF(2,4) = IF4[2][2][1][2] * 2.0;
  IF(2,5) = IF4[2][2][2][0] * 2.0;

  IF(3,0) = IF4[0][1][0][0];
  IF(3,1) = IF4[0][1][1][1];
  IF(3,2) = IF4[0][1][2][2];
  IF(3,3) = IF4[0][1][0][1] * 2.0;
  IF(3,4) = IF4[0][1][1][2] * 2.0;
  IF(3,5) = IF4[0][1][2][0] * 2.0;

  IF(4,0) = IF4[1][2][0][0];
  IF(4,1) = IF4[1][2][1][1];
  IF(4,2) = IF4[1][2][2][2];
  IF(4,3) = IF4[1][2][0][1] * 2.0;
  IF(4,4) = IF4[1][2][1][2] * 2.0;
  IF(4,5) = IF4[1][2][2][0] * 2.0;

  IF(5,0) = IF4[2][0][0][0];
  IF(5,1) = IF4[2][0][1][1];
  IF(5,2) = IF4[2][0][2][2];
  IF(5,3) = IF4[2][0][0][1] * 2.0;
  IF(5,4) = IF4[2][0][1][2] * 2.0;
  IF(5,5) = IF4[2][0][2][0] * 2.0;
#endif

  return; 
} // DRT::ELEMENTS::InvDesign::BuildIF

/*----------------------------------------------------------------------*
 |  Theta tensor                                              (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildTheta(LINALG::SerialDenseMatrix& Theta, 
                                       const LINALG::SerialDenseMatrix& F) const

{
#if 0 // for theory reference only
  Epetra_SerialDenseMatrix K(3,3);
  K(0,0) = 1.0;
  K(1,1) = 1.0;
  K(2,2) = 1.0;
  
  double Theta4[3][3][3][3];
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
        for (int l=0; l<3; ++l)
          Theta4[i][j][k][l] = 0.5 * ( K(i,l)*F(k,j) + K(j,l)*F(k,i) );
#endif
  
#if 1
  Theta(0,0) = ( F(0,0) );
  Theta(0,1) = ( F(1,0) );
  Theta(0,2) = ( F(2,0) );
  Theta(0,3) = 0.0;
  Theta(0,4) = 0.0;
  Theta(0,5) = 0.0;
  Theta(0,6) = 0.0;
  Theta(0,7) = 0.0;
  Theta(0,8) = 0.0;
  
  Theta(1,0) = 0.0;
  Theta(1,1) = 0.0;
  Theta(1,2) = 0.0;
  Theta(1,3) = ( F(0,1) );
  Theta(1,4) = ( F(1,1) );
  Theta(1,5) = ( F(2,1) );
  Theta(1,6) = 0.0;
  Theta(1,7) = 0.0;
  Theta(1,8) = 0.0;

  Theta(2,0) = 0.0;
  Theta(2,1) = 0.0;
  Theta(2,2) = 0.0;
  Theta(2,3) = 0.0;
  Theta(2,4) = 0.0;
  Theta(2,5) = 0.0;
  Theta(2,6) = ( F(0,2) );
  Theta(2,7) = ( F(1,2) );
  Theta(2,8) = ( F(2,2) );

  Theta(3,0) = ( F(0,1) ); 
  Theta(3,1) = ( F(1,1) );
  Theta(3,2) = ( F(2,1) );
  Theta(3,3) = ( F(0,0) );
  Theta(3,4) = ( F(1,0) );
  Theta(3,5) = ( F(2,0) );
  Theta(3,6) = 0.0;
  Theta(3,7) = 0.0;
  Theta(3,8) = 0.0;

  Theta(4,0) = 0.0; 
  Theta(4,1) = 0.0;
  Theta(4,2) = 0.0;
  Theta(4,3) = ( F(0,2) );
  Theta(4,4) = ( F(1,2) );
  Theta(4,5) = ( F(2,2) );
  Theta(4,6) = ( F(0,1) );
  Theta(4,7) = ( F(1,1) );
  Theta(4,8) = ( F(2,1) );

  Theta(5,0) = ( F(0,2) );
  Theta(5,1) = ( F(1,2) );
  Theta(5,2) = ( F(2,2) );
  Theta(5,3) = 0.0;
  Theta(5,4) = 0.0;
  Theta(5,5) = 0.0;
  Theta(5,6) = ( F(0,0) );
  Theta(5,7) = ( F(1,0) );
  Theta(5,8) = ( F(2,0) );
#else
  //                  i  j  k  l
  Theta(0,0) = Theta4[0][0][0][0];
  Theta(0,1) = Theta4[0][0][1][0];
  Theta(0,2) = Theta4[0][0][2][0];
  Theta(0,3) = Theta4[0][0][0][1];
  Theta(0,4) = Theta4[0][0][1][1];
  Theta(0,5) = Theta4[0][0][2][1];
  Theta(0,6) = Theta4[0][0][0][2];
  Theta(0,7) = Theta4[0][0][1][2];
  Theta(0,8) = Theta4[0][0][2][2];
  
  Theta(1,0) = Theta4[1][1][0][0];
  Theta(1,1) = Theta4[1][1][1][0];
  Theta(1,2) = Theta4[1][1][2][0];
  Theta(1,3) = Theta4[1][1][0][1];
  Theta(1,4) = Theta4[1][1][1][1];
  Theta(1,5) = Theta4[1][1][2][1];
  Theta(1,6) = Theta4[1][1][0][2];
  Theta(1,7) = Theta4[1][1][1][2];
  Theta(1,8) = Theta4[1][1][2][2];

  Theta(2,0) = Theta4[2][2][0][0];
  Theta(2,1) = Theta4[2][2][1][0];
  Theta(2,2) = Theta4[2][2][2][0];
  Theta(2,3) = Theta4[2][2][0][1];
  Theta(2,4) = Theta4[2][2][1][1];
  Theta(2,5) = Theta4[2][2][2][1];
  Theta(2,6) = Theta4[2][2][0][2];
  Theta(2,7) = Theta4[2][2][1][2];
  Theta(2,8) = Theta4[2][2][2][2];

  Theta(3,0) = Theta4[0][1][0][0] * 2.0; 
  Theta(3,1) = Theta4[0][1][1][0] * 2.0;
  Theta(3,2) = Theta4[0][1][2][0] * 2.0;
  Theta(3,3) = Theta4[0][1][0][1] * 2.0;
  Theta(3,4) = Theta4[0][1][1][1] * 2.0;
  Theta(3,5) = Theta4[0][1][2][1] * 2.0;
  Theta(3,6) = Theta4[0][1][0][2] * 2.0;
  Theta(3,7) = Theta4[0][1][1][2] * 2.0;
  Theta(3,8) = Theta4[0][1][2][2] * 2.0;

  Theta(4,0) = Theta4[1][2][0][0] * 2.0; 
  Theta(4,1) = Theta4[1][2][1][0] * 2.0;
  Theta(4,2) = Theta4[1][2][2][0] * 2.0;
  Theta(4,3) = Theta4[1][2][0][1] * 2.0;
  Theta(4,4) = Theta4[1][2][1][1] * 2.0;
  Theta(4,5) = Theta4[1][2][2][1] * 2.0;
  Theta(4,6) = Theta4[1][2][0][2] * 2.0;
  Theta(4,7) = Theta4[1][2][1][2] * 2.0;
  Theta(4,8) = Theta4[1][2][2][2] * 2.0;

  Theta(5,0) = Theta4[2][0][0][0] * 2.0;
  Theta(5,1) = Theta4[2][0][1][0] * 2.0;
  Theta(5,2) = Theta4[2][0][2][0] * 2.0;
  Theta(5,3) = Theta4[2][0][0][1] * 2.0;
  Theta(5,4) = Theta4[2][0][1][1] * 2.0;
  Theta(5,5) = Theta4[2][0][2][1] * 2.0;
  Theta(5,6) = Theta4[2][0][0][2] * 2.0;
  Theta(5,7) = Theta4[2][0][1][2] * 2.0;
  Theta(5,8) = Theta4[2][0][2][2] * 2.0;
#endif

  return; 
} // DRT::ELEMENTS::InvDesign::BuildTheta


/*----------------------------------------------------------------------*
 |  Ypsilon tensor                                            (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildYpsilon(LINALG::SerialDenseMatrix& Y, 
                                      const LINALG::SerialDenseMatrix& F,
                                      const LINALG::SerialDenseMatrix& S) const

{

  
#if 0  // for analogy to theory only
  Epetra_SerialDenseMatrix K(3,3);
  K(0,0) = 1.0;
  K(1,1) = 1.0;
  K(2,2) = 1.0;
  double Y4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int l=0; l<3; ++l)
      for (int p=0; p<3; ++p)
        for (int q=0; q<3; ++q)
        {
          Y4[k][l][p][q] = 0.0;
          for (int m=0; m<3; ++m)
            for (int n=0; n<3; ++n)
                Y4[k][l][p][q] += 
                  0.5 * (K(k,p) * ( K(m,q)*F(l,n) + K(n,q)*F(l,m) ) + K(l,p) * ( K(n,q)*F(k,m) + K(m,q)*F(k,n) )) * S(m,n);
         }
#endif                
                
                
#if 1
  // init Y to zero
  const int tdim = Y.M()*Y.N();
  for (int i=0; i<tdim; ++i) 
    Y.A()[i] = 0.0;
    
  Y(0,0)= (4.0*F(0,0))*S(0,0)   + ( 2.0*F(0,1) ) * S(0,1) + ( 2.0*F(0,2) ) * S(0,2) + ( 2.0*F(0,1) ) * S(1,0)                                                     + ( 2.0*F(0,2) ) * S(2,0)                                                     ;
  Y(0,3)=                       + ( 2.0*F(0,0) ) * S(0,1)                           + ( 2.0*F(0,0) ) * S(1,0) + ( 4.0*F(0,1) ) * S(1,1) + ( 2.0*F(0,2) ) * S(1,2)                           + ( 2.0*F(0,2) ) * S(2,1)                           ;
  Y(0,6)=                                                 + ( 2.0*F(0,0) ) * S(0,2)                                                     + ( 2.0*F(0,1) ) * S(1,2) + ( 2.0*F(0,0) ) * S(2,0) + ( 2.0*F(0,1) ) * S(2,1) + ( 4.0*F(0,2) ) * S(2,2) ;
  Y(1,1)= (4.0*F(1,0))*S(0,0)   + ( 2.0*F(1,1) ) * S(0,1) + ( 2.0*F(1,2) ) * S(0,2) + ( 2.0*F(1,1) ) * S(1,0)                                                     + ( 2.0*F(1,2) ) * S(2,0)                                                     ;
  Y(1,4)=                       + ( 2.0*F(1,0) ) * S(0,1)                           + ( 2.0*F(1,0) ) * S(1,0) + ( 4.0*F(1,1) ) * S(1,1) + ( 2.0*F(1,2) ) * S(1,2)                           + ( 2.0*F(1,2) ) * S(2,1)                           ;
  Y(1,7)=                                                 + ( 2.0*F(1,0) ) * S(0,2)                                                     + ( 2.0*F(1,1) ) * S(1,2) + ( 2.0*F(1,0) ) * S(2,0) + ( 2.0*F(1,1) ) * S(2,1) + ( 4.0*F(1,2) ) * S(2,2) ;
  Y(2,2)= (4.0*F(2,0))*S(0,0)   + ( 2.0*F(2,1) ) * S(0,1) + ( 2.0*F(2,2) ) * S(0,2) + ( 2.0*F(2,1) ) * S(1,0)                                                     + ( 2.0*F(2,2) ) * S(2,0)                                                     ;
  Y(2,5)=                       + ( 2.0*F(2,0) ) * S(0,1)                           + ( 2.0*F(2,0) ) * S(1,0) + ( 4.0*F(2,1) ) * S(1,1) + ( 2.0*F(2,2) ) * S(1,2)                           + ( 2.0*F(2,2) ) * S(2,1)                           ;
  Y(2,8)=                                                 + ( 2.0*F(2,0) ) * S(0,2)                                                     + ( 2.0*F(2,1) ) * S(1,2) + ( 2.0*F(2,0) ) * S(2,0) + ( 2.0*F(2,1) ) * S(2,1) + ( 4.0*F(2,2) ) * S(2,2) ;
  Y(3,0)= (2.0*F(1,0)) * S(0,0) + F(1,1) * S(0,1)         + F(1,2) * S(0,2)         + F(1,1) * S(1,0)                                                             + F(1,2) * S(2,0)                                                             ;
  Y(3,1)= (2.0*F(0,0)) * S(0,0) + F(0,1) * S(0,1)         + F(0,2) * S(0,2)         + F(0,1) * S(1,0)                                                             + F(0,2) * S(2,0)                                                             ;
  Y(3,3)=                       + F(1,0) * S(0,1)                                   + F(1,0) * S(1,0)         + ( 2.0*F(1,1) ) * S(1,1) + F(1,2) * S(1,2)                                   + F(1,2) * S(2,1)                                   ;
  Y(3,4)=                       + F(0,0) * S(0,1)                                   + F(0,0) * S(1,0)         + ( 2.0*F(0,1) ) * S(1,1) + F(0,2) * S(1,2)                                   + F(0,2) * S(2,1)                                   ;
  Y(3,6)=                                                 + F(1,0) * S(0,2)                                                             + F(1,1) * S(1,2)         + F(1,0) * S(2,0)         + F(1,1) * S(2,1)         + ( 2.0*F(1,2) ) * S(2,2) ;
  Y(3,7)=                                                 + F(0,0) * S(0,2)                                                             + F(0,1) * S(1,2)         + F(0,0) * S(2,0)         + F(0,1) * S(2,1)         + ( 2.0*F(0,2) ) * S(2,2) ;
  Y(4,1)= (2.0*F(2,0)) * S(0,0) + F(2,1) * S(0,1)         + F(2,2) * S(0,2)         + F(2,1) * S(1,0)                                                             + F(2,2) * S(2,0)                                                             ;
  Y(4,2)= (2.0*F(1,0)) * S(0,0) + F(1,1) * S(0,1)         + F(1,2) * S(0,2)         + F(1,1) * S(1,0)                                                             + F(1,2) * S(2,0)                                                             ;
  Y(4,4)=                       + F(2,0) * S(0,1)                                   + F(2,0) * S(1,0)         + ( 2.0*F(2,1) ) * S(1,1) + F(2,2) * S(1,2)                                   + F(2,2) * S(2,1)                                   ;
  Y(4,5)=                       + F(1,0) * S(0,1)                                   + F(1,0) * S(1,0)         + ( 2.0*F(1,1) ) * S(1,1) + F(1,2) * S(1,2)                                   + F(1,2) * S(2,1)                                   ;
  Y(4,7)=                                                 + F(2,0) * S(0,2)                                                             + F(2,1) * S(1,2)         + F(2,0) * S(2,0)         + F(2,1) * S(2,1)         + ( 2.0*F(2,2) ) * S(2,2) ;
  Y(4,8)=                                                 + F(1,0) * S(0,2)                                                             + F(1,1) * S(1,2)         + F(1,0) * S(2,0)         + F(1,1) * S(2,1)         + ( 2.0*F(1,2) ) * S(2,2) ;
  Y(5,0)= (2.0*F(2,0)) * S(0,0) + F(2,1) * S(0,1)         + F(2,2) * S(0,2)         + F(2,1) * S(1,0)                                                             + F(2,2) * S(2,0)                                                             ;
  Y(5,2)= (2.0*F(0,0)) * S(0,0) + F(0,1) * S(0,1)         + F(0,2) * S(0,2)         + F(0,1) * S(1,0)                                                             + F(0,2) * S(2,0)                                                             ;
  Y(5,3)=                       + F(2,0) * S(0,1)                                   + F(2,0) * S(1,0)         + ( 2.0*F(2,1) ) * S(1,1) + F(2,2) * S(1,2)                                   + F(2,2) * S(2,1)                                   ;
  Y(5,5)=                       + F(0,0) * S(0,1)                                   + F(0,0) * S(1,0)         + ( 2.0*F(0,1) ) * S(1,1) + F(0,2) * S(1,2)                                   + F(0,2) * S(2,1)                                   ;
  Y(5,6)=                                                 + F(2,0) * S(0,2)                                                             + F(2,1) * S(1,2)         + F(2,0) * S(2,0)         + F(2,1) * S(2,1)         + ( 2.0*F(2,2) ) * S(2,2) ;
  Y(5,8)=                                                 + F(0,0) * S(0,2)                                                             + F(0,1) * S(1,2)         + F(0,0) * S(2,0)         + F(0,1) * S(2,1)         + ( 2.0*F(0,2) ) * S(2,2) ;
  Y.Scale(0.5);

#else
  //          k  l  p  q
  Y(0,0) = Y4[0][0][0][0];
  Y(0,1) = Y4[0][0][1][0];
  Y(0,2) = Y4[0][0][2][0];
  Y(0,3) = Y4[0][0][0][1];
  Y(0,4) = Y4[0][0][1][1];
  Y(0,5) = Y4[0][0][2][1];
  Y(0,6) = Y4[0][0][0][2];
  Y(0,7) = Y4[0][0][1][2];
  Y(0,8) = Y4[0][0][2][2];
  
  Y(1,0) = Y4[1][1][0][0];
  Y(1,1) = Y4[1][1][1][0];
  Y(1,2) = Y4[1][1][2][0];
  Y(1,3) = Y4[1][1][0][1];
  Y(1,4) = Y4[1][1][1][1];
  Y(1,5) = Y4[1][1][2][1];
  Y(1,6) = Y4[1][1][0][2];
  Y(1,7) = Y4[1][1][1][2];
  Y(1,8) = Y4[1][1][2][2];

  Y(2,0) = Y4[2][2][0][0];
  Y(2,1) = Y4[2][2][1][0];
  Y(2,2) = Y4[2][2][2][0];
  Y(2,3) = Y4[2][2][0][1];
  Y(2,4) = Y4[2][2][1][1];
  Y(2,5) = Y4[2][2][2][1];
  Y(2,6) = Y4[2][2][0][2];
  Y(2,7) = Y4[2][2][1][2];
  Y(2,8) = Y4[2][2][2][2];

  Y(3,0) = Y4[0][1][0][0];
  Y(3,1) = Y4[0][1][1][0];
  Y(3,2) = Y4[0][1][2][0];
  Y(3,3) = Y4[0][1][0][1];
  Y(3,4) = Y4[0][1][1][1];
  Y(3,5) = Y4[0][1][2][1];
  Y(3,6) = Y4[0][1][0][2];
  Y(3,7) = Y4[0][1][1][2];
  Y(3,8) = Y4[0][1][2][2];

  Y(4,0) = Y4[1][2][0][0];
  Y(4,1) = Y4[1][2][1][0];
  Y(4,2) = Y4[1][2][2][0];
  Y(4,3) = Y4[1][2][0][1];
  Y(4,4) = Y4[1][2][1][1];
  Y(4,5) = Y4[1][2][2][1];
  Y(4,6) = Y4[1][2][0][2];
  Y(4,7) = Y4[1][2][1][2];
  Y(4,8) = Y4[1][2][2][2];

  Y(5,0) = Y4[2][0][0][0];
  Y(5,1) = Y4[2][0][1][0];
  Y(5,2) = Y4[2][0][2][0];
  Y(5,3) = Y4[2][0][0][1];
  Y(5,4) = Y4[2][0][1][1];
  Y(5,5) = Y4[2][0][2][1];
  Y(5,6) = Y4[2][0][0][2];
  Y(5,7) = Y4[2][0][1][2];
  Y(5,8) = Y4[2][0][2][2];
#endif                

  return; 
} // DRT::ELEMENTS::InvDesign::BuildYpsilon


/*----------------------------------------------------------------------*
 |  hex8 integration method                                    (public) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::soh8_nlnstiffmass(
      DRT::ELEMENTS::So_hex8*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    ///< element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     ///< element mass matrix
      Epetra_SerialDenseVector* force,          ///< element internal force vector
      Epetra_SerialDenseMatrix* elestress,      ///< stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy,         ///< stress output option
      const bool                euler_almansi)  ///< strain output option
{
   const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
  
  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8,NUMDIM_SOH8);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xcurr(i,0) = ele->Nodes()[i]->X()[0];
    xcurr(i,1) = ele->Nodes()[i]->X()[1];
    xcurr(i,2) = ele->Nodes()[i]->X()[2];

    xrefe(i,0) = xcurr(i,0) + disp[i*NODDOF_SOH8+0];
    xrefe(i,1) = xcurr(i,1) + disp[i*NODDOF_SOH8+1];
    xrefe(i,2) = xcurr(i,2) + disp[i*NODDOF_SOH8+2];
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) 
  {
    //----------------------------------- get inverse of Jacobian mapping
    const double              detj = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];
    
    //------------------------------ compute derivs wrt to spatial coords
    LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8,NUMNOD_SOH8);
    n_xyz.Multiply('N','N',1.0,invj,int_hex8.deriv_gp[gp],0.0);
  
    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::SerialDenseMatrix f(NUMDIM_SOH8,NUMDIM_SOH8);
    f.Multiply('T','T',1.0,xrefe,n_xyz,0.0);
    
    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::SerialDenseMatrix F(f);
    const double j = LINALG::NonsymInverse3x3(F);
    
    //------------------------------------ build F^T as vector 9x1
    LINALG::SerialDenseVector FTvec(9);
    FTvec(0) = F(0,0);
    FTvec(1) = F(0,1);
    FTvec(2) = F(0,2);
    FTvec(3) = F(1,0);
    FTvec(4) = F(1,1);
    FTvec(5) = F(1,2);
    FTvec(6) = F(2,0);
    FTvec(7) = F(2,1);
    FTvec(8) = F(2,2);
    
    //--------------------------------------------- build operator Lambda
    LINALG::SerialDenseMatrix Lambda(9,9);
    BuildLambda(Lambda,F);
    
    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::SerialDenseMatrix IF(6,6);
    BuildIF(IF,F);
    
    //---------------------------------------------- build operator Theta
    LINALG::SerialDenseMatrix Theta(6,9);
    BuildTheta(Theta,F);
    
    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8,NUMDIM_SOH8);
    cauchygreen.Multiply('T','N',1.0,F,F,0.0);
    
    //-- Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain(NUMSTR_SOH8);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::SerialDenseMatrix B(NUMSTR_SOH8,NUMDOF_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i) 
    {
      B(0,NODDOF_SOH8*i+0) = n_xyz(0,i);
      B(0,NODDOF_SOH8*i+1) = 0.0;
      B(0,NODDOF_SOH8*i+2) = 0.0;

      B(1,NODDOF_SOH8*i+0) = 0.0;
      B(1,NODDOF_SOH8*i+1) = n_xyz(1,i);
      B(1,NODDOF_SOH8*i+2) = 0.0;

      B(2,NODDOF_SOH8*i+0) = 0.0;
      B(2,NODDOF_SOH8*i+1) = 0.0;
      B(2,NODDOF_SOH8*i+2) = n_xyz(2,i);

      B(3,NODDOF_SOH8*i+0) = n_xyz(1,i);
      B(3,NODDOF_SOH8*i+1) = n_xyz(0,i);
      B(3,NODDOF_SOH8*i+2) = 0.0;

      B(4,NODDOF_SOH8*i+0) = 0.0;
      B(4,NODDOF_SOH8*i+1) = n_xyz(2,i);
      B(4,NODDOF_SOH8*i+2) = n_xyz(1,i);

      B(5,NODDOF_SOH8*i+0) = n_xyz(2,i);
      B(5,NODDOF_SOH8*i+1) = 0.0;
      B(5,NODDOF_SOH8*i+2) = n_xyz(0,i);
    }
    
    //--------------------------- build N_x operator (wrt spatial config)
    Epetra_SerialDenseMatrix N_x(9,NUMDOF_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i) 
    {
      N_x(0,3*i+0) = n_xyz(0,i);
      N_x(1,3*i+1) = n_xyz(0,i);
      N_x(2,3*i+2) = n_xyz(0,i);
      
      N_x(3,3*i+0) = n_xyz(1,i);
      N_x(4,3*i+1) = n_xyz(1,i);
      N_x(5,3*i+2) = n_xyz(1,i);
      
      N_x(6,3*i+0) = n_xyz(2,i);
      N_x(7,3*i+1) = n_xyz(2,i);
      N_x(8,3*i+2) = n_xyz(2,i);
    }
    
    //------------------------------------------------- call material law
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress(NUMSTR_SOH8);
    double density;
    ele->soh8_mat_sel(&stress,&cmat,&density,&glstrain,&F,gp,params);

    //------------------------------------------- compute cauchy stresses
    Epetra_SerialDenseVector cstress(NUMSTR_SOH8);
    cstress.Multiply('N','N',j,IF,stress,0.0);
    
    //--------------------------------------- output strains and stresses
    if (elestrain)
    {
      if (!euler_almansi)
      {
        for (int i = 0; i < 3; ++i) 
          (*elestrain)(gp,i) = glstrain(i);
        for (int i = 3; i < 6; ++i) 
          (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
      else
      {
        // rewriting Green-Lagrange strains in matrix format
        LINALG::SerialDenseMatrix gl(NUMDIM_SOH8,NUMDIM_SOH8);
        gl(0,0) = glstrain(0);
        gl(0,1) = 0.5*glstrain(3);
        gl(0,2) = 0.5*glstrain(5);
        gl(1,0) = gl(0,1);
        gl(1,1) = glstrain(1);
        gl(1,2) = 0.5*glstrain(4);
        gl(2,0) = gl(0,2);
        gl(2,1) = gl(1,2);
        gl(2,2) = glstrain(2);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOH8,NUMDIM_SOH8);
        LINALG::SerialDenseMatrix euler_almansi(NUMDIM_SOH8,NUMDIM_SOH8);
        temp.Multiply('N','N',1.0,gl,f,0.0);
        euler_almansi.Multiply('T','N',1.0,f,temp,0.0);

        (*elestrain)(gp,0) = euler_almansi(0,0);
        (*elestrain)(gp,1) = euler_almansi(1,1);
        (*elestrain)(gp,2) = euler_almansi(2,2);
        (*elestrain)(gp,3) = euler_almansi(0,1);
        (*elestrain)(gp,4) = euler_almansi(1,2);
        (*elestrain)(gp,5) = euler_almansi(0,2);
      }
    }
    
    if (elestress)
    {
      if (!cauchy)
        for (int i = 0; i < 6; ++i) 
          (*elestress)(gp,i) = stress(i);
      else
        for (int i = 0; i < 6; ++i) 
          (*elestress)(gp,i) = cstress(i);
    }
        
    //-------------------------------------------- build operator Ypsilon
    LINALG::SerialDenseMatrix S(3,3);
    S(0,0) = stress(0);
    S(0,1) = stress(3);
    S(0,2) = stress(5);
    S(1,0) = stress(3);
    S(1,1) = stress(1);
    S(1,2) = stress(4);
    S(2,0) = stress(5);
    S(2,1) = stress(4);
    S(2,2) = stress(2);
    LINALG::SerialDenseMatrix Ypsilon(6,9);
    BuildYpsilon(Ypsilon,F,S);
    
    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * int_hex8.weights(gp);

    //------------------------------------------- assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T','N',intfac,B,cstress,1.0);
    }
    
    //*******************************************************************
    if (stiffmatrix)
    {
      LINALG::SerialDenseMatrix sum(6,9);
      LINALG::SerialDenseMatrix tmp2(1,9);

      //================================== 1.) B^T * cstress * F^T * N_x
      sum.Multiply('N','T',1.0,cstress,FTvec,0.0);
      
      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      sum.Multiply('N','N',-j,Ypsilon,Lambda,1.0);
      
      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
      LINALG::SerialDenseMatrix sixnine(6,9);
      LINALG::SerialDenseMatrix sixsix(6,6);
      sixnine.Multiply('N','N',1.0,Theta,Lambda,0.0);
      sixsix.Multiply('N','N',1.0,IF,cmat,0.0);
      sum.Multiply('N','N',-j,sixsix,sixnine,1.0);

      //================ put everything together: K = intfac * B*T * sum * N_x
      LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
      tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
      (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);

    } // if (stiffmatrix)
    //*******************************************************************

    
    //*******************************************************************
    // Strictly, inverse design analysis is stationary and should not have
    // a mass term. Loosely, if no Dirichlet-BCs are present a small
    // mass term might be used. Note that we use density by unit deformed
    // volume here!
    if (massmatrix)
    {
      const double fac = density * detj * int_hex8.weights(gp);
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) 
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) 
        {
          const double massfactor = (int_hex8.shapefct_gp[gp])(inod) 
                                  * (int_hex8.shapefct_gp[gp])(jnod)
                                  * fac;     
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
    } // if (massmatrix)
    //*******************************************************************

  } // for (int gp=0; gp<NUMGPT_SOH8; ++gp) 
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  
  
} // DRT::ELEMENTS::InvDesign::soh8_nlnstiffmass



/*----------------------------------------------------------------------*
 |  compute and store material configuration                  (private) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::soh8_StoreMaterialConfiguration(
                                              DRT::ELEMENTS::So_hex8* ele, 
                                              const vector<double>& disp)
{
  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8,3);  // material coord. of element
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8,3);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xcurr(i,0) = ele->Nodes()[i]->X()[0];
    xcurr(i,1) = ele->Nodes()[i]->X()[1];
    xcurr(i,2) = ele->Nodes()[i]->X()[2];

    xrefe(i,0) = xcurr(i,0) + disp[i*NODDOF_SOH8+0];
    xrefe(i,1) = xcurr(i,1) + disp[i*NODDOF_SOH8+1];
    xrefe(i,2) = xcurr(i,2) + disp[i*NODDOF_SOH8+2];
  }
  
  LINALG::SerialDenseMatrix invJ(3,3);
  LINALG::SerialDenseMatrix N_XYZ(3,NUMNOD_SOH8);
  LINALG::SerialDenseMatrix F(3,3);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) 
  {
    // compute invJ and detJ
    invJ.Multiply('N','N',1.0,int_hex8.deriv_gp[gp],xrefe,0.0);
    detJ_[gp] = LINALG::NonsymInverse3x3(invJ);
    
    // compute F
    N_XYZ.Multiply('N','N',1.0,invJ,int_hex8.deriv_gp[gp],0.0);
    F.Multiply('T','T',1.0,xcurr,N_XYZ,0.0);
    
    // put stuff into storage
    MatrixtoStorage(gp,invJ,JHistory());
    MatrixtoStorage(gp,F,FHistory());
  }  

  return;
}


/*----------------------------------------------------------------------*
 |  weg6 integration method                                    (public) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::sow6_nlnstiffmass(
      DRT::ELEMENTS::So_weg6*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    ///< element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     ///< element mass matrix
      Epetra_SerialDenseVector* force,          ///< element internal force vector
      Epetra_SerialDenseMatrix* elestress,      ///< element stresses
      Epetra_SerialDenseMatrix* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy,         ///< stress output option
      const bool                euler_almansi)  ///< strain output option
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  Epetra_SerialDenseMatrix* shapefct; //[NUMNOD_WEG6][NUMGPT_WEG6]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_WEG6*NUMDIM][NUMNOD_WEG6]
/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_WEG6]
  ele->sow6_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::SerialDenseMatrix xrefe(NUMNOD_WEG6,NUMDIM_WEG6);  // material coord. of element
  LINALG::SerialDenseMatrix xcurr(NUMNOD_WEG6,NUMDIM_WEG6);  // current  coord. of element
  for (int i=0; i<NUMNOD_WEG6; ++i)
  {
    xcurr(i,0) = ele->Nodes()[i]->X()[0];
    xcurr(i,1) = ele->Nodes()[i]->X()[1];
    xcurr(i,2) = ele->Nodes()[i]->X()[2];

    xrefe(i,0) = xcurr(i,0) + disp[i*NODDOF_WEG6+0];
    xrefe(i,1) = xcurr(i,1) + disp[i*NODDOF_WEG6+1];
    xrefe(i,2) = xcurr(i,2) + disp[i*NODDOF_WEG6+2];
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  for (int gp=0; gp<NUMGPT_WEG6; ++gp)
  {
    // get submatrix of deriv at actual gp
    LINALG::SerialDenseMatrix deriv_gp(NUMDIM_WEG6,NUMGPT_WEG6);
    for (int m=0; m<NUMDIM_WEG6; ++m)
      for (int n=0; n<NUMGPT_WEG6; ++n)
        deriv_gp(m,n)=(*deriv)(NUMDIM_WEG6*gp+m,n);

    const double detj              = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];
    
    //------------------------------ compute derivs wrt to spatial coords
    LINALG::SerialDenseMatrix n_xyz(NUMDIM_WEG6,NUMNOD_WEG6);
    n_xyz.Multiply('N','N',1.0,invj,deriv_gp,0.0);
    
    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::SerialDenseMatrix f(NUMDIM_WEG6,NUMDIM_WEG6);
    f.Multiply('T','T',1.0,xrefe,n_xyz,0.0);
    
    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::SerialDenseMatrix F(f);
    const double j = LINALG::NonsymInverse3x3(F);

    //------------------------------------ build F^T as vector 9x1
    LINALG::SerialDenseVector FTvec(9);
    FTvec(0) = F(0,0);
    FTvec(1) = F(0,1);
    FTvec(2) = F(0,2);
    FTvec(3) = F(1,0);
    FTvec(4) = F(1,1);
    FTvec(5) = F(1,2);
    FTvec(6) = F(2,0);
    FTvec(7) = F(2,1);
    FTvec(8) = F(2,2);
    
    //--------------------------------------------- build operator Lambda
    LINALG::SerialDenseMatrix Lambda(9,9);
    BuildLambda(Lambda,F);
    
    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::SerialDenseMatrix IF(6,6);
    BuildIF(IF,F);
    
    //---------------------------------------------- build operator Theta
    LINALG::SerialDenseMatrix Theta(6,9);
    BuildTheta(Theta,F);

    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_WEG6,NUMDIM_WEG6);
    cauchygreen.Multiply('T','N',1.0,F,F,0.0);
    
    //-- Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain(NUMSTR_WEG6);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);
            
    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::SerialDenseMatrix B(NUMSTR_WEG6,NUMDOF_WEG6);
    for (int i=0; i<NUMNOD_WEG6; ++i) 
    {
      B(0,NODDOF_WEG6*i+0) = n_xyz(0,i);
      B(0,NODDOF_WEG6*i+1) = 0.0;
      B(0,NODDOF_WEG6*i+2) = 0.0;

      B(1,NODDOF_WEG6*i+0) = 0.0;
      B(1,NODDOF_WEG6*i+1) = n_xyz(1,i);
      B(1,NODDOF_WEG6*i+2) = 0.0;

      B(2,NODDOF_WEG6*i+0) = 0.0;
      B(2,NODDOF_WEG6*i+1) = 0.0;
      B(2,NODDOF_WEG6*i+2) = n_xyz(2,i);

      B(3,NODDOF_WEG6*i+0) = n_xyz(1,i);
      B(3,NODDOF_WEG6*i+1) = n_xyz(0,i);
      B(3,NODDOF_WEG6*i+2) = 0.0;

      B(4,NODDOF_WEG6*i+0) = 0.0;
      B(4,NODDOF_WEG6*i+1) = n_xyz(2,i);
      B(4,NODDOF_WEG6*i+2) = n_xyz(1,i);

      B(5,NODDOF_WEG6*i+0) = n_xyz(2,i);
      B(5,NODDOF_WEG6*i+1) = 0.0;
      B(5,NODDOF_WEG6*i+2) = n_xyz(0,i);
    }
    
    //--------------------------- build N_x operator (wrt spatial config)
    Epetra_SerialDenseMatrix N_x(9,NUMDOF_WEG6);
    for (int i=0; i<NUMNOD_WEG6; ++i) 
    {
      N_x(0,3*i+0) = n_xyz(0,i);
      N_x(1,3*i+1) = n_xyz(0,i);
      N_x(2,3*i+2) = n_xyz(0,i);
      
      N_x(3,3*i+0) = n_xyz(1,i);
      N_x(4,3*i+1) = n_xyz(1,i);
      N_x(5,3*i+2) = n_xyz(1,i);
      
      N_x(6,3*i+0) = n_xyz(2,i);
      N_x(7,3*i+1) = n_xyz(2,i);
      N_x(8,3*i+2) = n_xyz(2,i);
    }
    
    //------------------------------------------------- call material law
    Epetra_SerialDenseMatrix cmat(NUMSTR_WEG6,NUMSTR_WEG6);
    Epetra_SerialDenseVector stress(NUMSTR_WEG6);
    double density;
    ele->sow6_mat_sel(&stress,&cmat,&density,&glstrain,&F,gp,params);
    
    //------------------------------------------- compute cauchy stresses
    Epetra_SerialDenseVector cstress(NUMSTR_WEG6);
    cstress.Multiply('N','N',j,IF,stress,0.0);

    //--------------------------------------- output strains and stresses
    if (elestrain)
    {
      if (!euler_almansi)
      {
        for (int i = 0; i < 3; ++i) 
          (*elestrain)(gp,i) = glstrain(i);
        for (int i = 3; i < 6; ++i) 
          (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
      else
      {
        // rewriting Green-Lagrange strains in matrix format
        LINALG::SerialDenseMatrix gl(NUMDIM_WEG6,NUMDIM_WEG6);
        gl(0,0) = glstrain(0);
        gl(0,1) = 0.5*glstrain(3);
        gl(0,2) = 0.5*glstrain(5);
        gl(1,0) = gl(0,1);
        gl(1,1) = glstrain(1);
        gl(1,2) = 0.5*glstrain(4);
        gl(2,0) = gl(0,2);
        gl(2,1) = gl(1,2);
        gl(2,2) = glstrain(2);

        LINALG::SerialDenseMatrix temp(NUMDIM_WEG6,NUMDIM_WEG6);
        LINALG::SerialDenseMatrix euler_almansi(NUMDIM_WEG6,NUMDIM_WEG6);
        temp.Multiply('N','N',1.0,gl,f,0.0);
        euler_almansi.Multiply('T','N',1.0,f,temp,0.0);

        (*elestrain)(gp,0) = euler_almansi(0,0);
        (*elestrain)(gp,1) = euler_almansi(1,1);
        (*elestrain)(gp,2) = euler_almansi(2,2);
        (*elestrain)(gp,3) = euler_almansi(0,1);
        (*elestrain)(gp,4) = euler_almansi(1,2);
        (*elestrain)(gp,5) = euler_almansi(0,2);
      }
    }
    
    if (elestress)
    {
      if (!cauchy)
        for (int i = 0; i < 6; ++i) 
          (*elestress)(gp,i) = stress(i);
      else
        for (int i = 0; i < 6; ++i) 
          (*elestress)(gp,i) = cstress(i);
    }

    //-------------------------------------------- build operator Ypsilon
    LINALG::SerialDenseMatrix S(3,3);
    S(0,0) = stress(0);
    S(0,1) = stress(3);
    S(0,2) = stress(5);
    S(1,0) = stress(3);
    S(1,1) = stress(1);
    S(1,2) = stress(4);
    S(2,0) = stress(5);
    S(2,1) = stress(4);
    S(2,2) = stress(2);
    LINALG::SerialDenseMatrix Ypsilon(6,9);
    BuildYpsilon(Ypsilon,F,S);

    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * (*weights)(gp);

    //------------------------------------------- assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T','N',intfac,B,cstress,1.0);
    }

    //*******************************************************************
    if (stiffmatrix)
    {
      LINALG::SerialDenseMatrix sum(6,9);
      LINALG::SerialDenseMatrix tmp2(1,9);

      //================================== 1.) B^T * cstress * F^T * N_x
      sum.Multiply('N','T',1.0,cstress,FTvec,0.0);
      
      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      sum.Multiply('N','N',-j,Ypsilon,Lambda,1.0);
      
      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
      LINALG::SerialDenseMatrix sixnine(6,9);
      LINALG::SerialDenseMatrix sixsix(6,6);
      sixnine.Multiply('N','N',1.0,Theta,Lambda,0.0);
      sixsix.Multiply('N','N',1.0,IF,cmat,0.0);
      sum.Multiply('N','N',-j,sixsix,sixnine,1.0);

      //================ put everything together: K = intfac * B*T * sum * N_x
      LINALG::SerialDenseMatrix tmp3(6,NUMDOF_WEG6);
      tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
      (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);

    } // if (stiffmatrix)
    //*******************************************************************

    //*******************************************************************
    // Strictly, inverse design analysis is stationary and should not have
    // a mass term. Loosely, if no Dirichlet-BCs are present a small
    // mass term might be used. Note that we use density by unit deformed
    // volume here!
    if (massmatrix)
    { // evaluate mass matrix +++++++++++++++++++++++++
      // integrate concistent mass matrix
      for (int inod=0; inod<NUMNOD_WEG6; ++inod) 
      {
        for (int jnod=0; jnod<NUMNOD_WEG6; ++jnod) 
        {
          double massfactor = (*shapefct)(inod,gp) * density * (*shapefct)(jnod,gp)
                            * detj * (*weights)(gp);     // intermediate factor
          (*massmatrix)(NUMDIM_WEG6*inod+0,NUMDIM_WEG6*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_WEG6*inod+1,NUMDIM_WEG6*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_WEG6*inod+2,NUMDIM_WEG6*jnod+2) += massfactor;
        }
      }
    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  } // for (int gp=0; gp<NUMGPT_WEG6; ++gp) 
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  return;
} // DRT::ELEMENTS::InvDesign::sow6_nlnstiffmass


/*----------------------------------------------------------------------*
 |  compute and store material configuration                  (private) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::sow6_StoreMaterialConfiguration(
                                       DRT::ELEMENTS::So_weg6* ele,
                                       const vector<double>& disp)
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
** ============================================================================*/
/* pointer to (static) shape function array
 * for each node, evaluated at each gp*/
  Epetra_SerialDenseMatrix* shapefct; //[NUMNOD_WEG6][NUMGPT_WEG6]
/* pointer to (static) shape function derivatives array
 * for each node wrt to each direction, evaluated at each gp*/
  Epetra_SerialDenseMatrix* deriv;    //[NUMGPT_WEG6*NUMDIM][NUMNOD_WEG6]
/* pointer to (static) weight factors at each gp */
  Epetra_SerialDenseVector* weights;  //[NUMGPT_WEG6]
  ele->sow6_shapederiv(&shapefct,&deriv,&weights);   // call to evaluate
/* ============================================================================*/

  LINALG::SerialDenseMatrix xrefe(NUMNOD_WEG6,3);  // material coord. of element
  LINALG::SerialDenseMatrix xcurr(NUMNOD_WEG6,3);  // current  coord. of element
  for (int i=0; i<NUMNOD_WEG6; ++i)
  {
    xcurr(i,0) = ele->Nodes()[i]->X()[0];
    xcurr(i,1) = ele->Nodes()[i]->X()[1];
    xcurr(i,2) = ele->Nodes()[i]->X()[2];

    xrefe(i,0) = xcurr(i,0) + disp[i*NODDOF_WEG6+0];
    xrefe(i,1) = xcurr(i,1) + disp[i*NODDOF_WEG6+1];
    xrefe(i,2) = xcurr(i,2) + disp[i*NODDOF_WEG6+2];
  }
  
  LINALG::SerialDenseMatrix invJ(3,3);
  LINALG::SerialDenseMatrix N_XYZ(3,NUMNOD_WEG6);
  LINALG::SerialDenseMatrix F(3,3);
  for (int gp=0; gp<NUMGPT_WEG6; ++gp) 
  {
    // get submatrix of deriv at actual gp
    LINALG::SerialDenseMatrix deriv_gp(NUMDIM_WEG6,NUMGPT_WEG6);
    for (int m=0; m<NUMDIM_WEG6; ++m) 
      for (int n=0; n<NUMGPT_WEG6; ++n) 
        deriv_gp(m,n)=(*deriv)(NUMDIM_WEG6*gp+m,n);

    // compute invJ and detJ
    invJ.Multiply('N','N',1.0,deriv_gp,xrefe,0.0);
    detJ_[gp] = LINALG::NonsymInverse3x3(invJ);
    
    // compute F
    N_XYZ.Multiply('N','N',1.0,invJ,deriv_gp,0.0);
    F.Multiply('T','T',1.0,xcurr,N_XYZ,0.0);
    
    // put stuff into storage
    MatrixtoStorage(gp,invJ,JHistory());
    MatrixtoStorage(gp,F,FHistory());
  }  
  return;
}




/*----------------------------------------------------------------------*
 |  tet4 integration method                                    (public) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::so_tet4_nlnstiffmass(
      DRT::ELEMENTS::So_tet4*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    ///< element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     ///< element mass matrix
      Epetra_SerialDenseVector* force,          ///< element internal force vector
      Epetra_SerialDenseMatrix* elestress,      ///< stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      ///< strains at GP
      struct _MATERIAL*         material,       ///< element material data
      const bool                cauchy,         ///< stress output options
      const bool                ea)
{
/* =============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
** =============================================================================*/
  const static DRT::ELEMENTS::Integrator_tet4_1point tet4_dis;
  double density;

  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::SerialDenseMatrix xdisp(NUMNOD_SOTET4,NUMDIM_SOTET4);
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xdisp(i,0) = disp[i*NODDOF_SOTET4+0];
    xdisp(i,1) = disp[i*NODDOF_SOTET4+1];
    xdisp(i,2) = disp[i*NODDOF_SOTET4+2];
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  const double detj  = ele->V_;
  for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  {
    //----------------------------------- get inverse of Jacobian mapping
    //------------------------------ compute derivs wrt to spatial coords
    Epetra_SerialDenseMatrix& n_xyz = ele->nxyz_[gp];

    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::SerialDenseMatrix f(NUMDIM_SOTET4,NUMDIM_SOTET4);
    f.Multiply('T','N',1.0,xdisp,n_xyz,0.0);
    f(0,0)+=1;
    f(1,1)+=1;
    f(2,2)+=1;
    
    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::SerialDenseMatrix F(f);
    const double j = LINALG::NonsymInverse3x3(F);

    //------------------------------------ build F^T as vector 9x1
    LINALG::SerialDenseVector FTvec(9);
    FTvec(0) = F(0,0);
    FTvec(1) = F(0,1);
    FTvec(2) = F(0,2);
    FTvec(3) = F(1,0);
    FTvec(4) = F(1,1);
    FTvec(5) = F(1,2);
    FTvec(6) = F(2,0);
    FTvec(7) = F(2,1);
    FTvec(8) = F(2,2);
    
    //--------------------------------------------- build operator Lambda
    LINALG::SerialDenseMatrix Lambda(9,9);
    BuildLambda(Lambda,F);
    
    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::SerialDenseMatrix IF(6,6);
    BuildIF(IF,F);
    
    //---------------------------------------------- build operator Theta
    LINALG::SerialDenseMatrix Theta(6,9);
    BuildTheta(Theta,F);
    
    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOTET4,NUMDIM_SOTET4);
    cauchygreen.Multiply('T','N',1.0,F,F,0.0);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain(NUMSTR_SOTET4);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);
    
    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::SerialDenseMatrix B(NUMSTR_SOTET4,NUMDOF_SOTET4);
    for (int i=0; i<NUMNOD_SOTET4; ++i) 
    {
      B(0,NODDOF_SOTET4*i+0) = n_xyz(i,0);
      B(0,NODDOF_SOTET4*i+1) = 0.0;
      B(0,NODDOF_SOTET4*i+2) = 0.0;

      B(1,NODDOF_SOTET4*i+0) = 0.0;
      B(1,NODDOF_SOTET4*i+1) = n_xyz(i,1);
      B(1,NODDOF_SOTET4*i+2) = 0.0;

      B(2,NODDOF_SOTET4*i+0) = 0.0;
      B(2,NODDOF_SOTET4*i+1) = 0.0;
      B(2,NODDOF_SOTET4*i+2) = n_xyz(i,2);

      B(3,NODDOF_SOTET4*i+0) = n_xyz(i,1);
      B(3,NODDOF_SOTET4*i+1) = n_xyz(i,0);
      B(3,NODDOF_SOTET4*i+2) = 0.0;

      B(4,NODDOF_SOTET4*i+0) = 0.0;
      B(4,NODDOF_SOTET4*i+1) = n_xyz(i,2);
      B(4,NODDOF_SOTET4*i+2) = n_xyz(i,1);

      B(5,NODDOF_SOTET4*i+0) = n_xyz(i,2);
      B(5,NODDOF_SOTET4*i+1) = 0.0;
      B(5,NODDOF_SOTET4*i+2) = n_xyz(i,0);
    }
    
    //--------------------------- build N_x operator (wrt spatial config)
    Epetra_SerialDenseMatrix N_x(9,NUMDOF_SOTET4);
    for (int i=0; i<NUMNOD_SOTET4; ++i) 
    {
      N_x(0,3*i+0) = n_xyz(i,0);
      N_x(1,3*i+1) = n_xyz(i,0);
      N_x(2,3*i+2) = n_xyz(i,0);
      
      N_x(3,3*i+0) = n_xyz(i,1);
      N_x(4,3*i+1) = n_xyz(i,1);
      N_x(5,3*i+2) = n_xyz(i,1);
      
      N_x(6,3*i+0) = n_xyz(i,2);
      N_x(7,3*i+1) = n_xyz(i,2);
      N_x(8,3*i+2) = n_xyz(i,2);
    }

    //------------------------------------------------- call material law
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOTET4,NUMSTR_SOTET4);
    Epetra_SerialDenseVector stress(NUMSTR_SOTET4);
    ele->so_tet4_mat_sel(&stress,&cmat,&density,&glstrain, &F, gp);

    //------------------------------------------- compute cauchy stresses
    Epetra_SerialDenseVector cstress(NUMSTR_SOTET4);
    cstress.Multiply('N','N',j,IF,stress,0.0);

    //--------------------------------------- output strains and stresses
    if (elestrain)
    {
      if (!ea) // output Green-Lagrange strains
      {
        for (int i = 0; i < 3; ++i)
          (*elestrain)(gp,i) = glstrain(i);
        for (int i = 3; i < 6; ++i)
          (*elestrain)(gp,i) = 0.5 * glstrain(i);
      }
      else
      {
        // rewriting Green-Lagrange strains in matrix format
        LINALG::SerialDenseMatrix gl(NUMDIM_SOTET4,NUMDIM_SOTET4);
        gl(0,0) = glstrain(0);
        gl(0,1) = 0.5*glstrain(3);
        gl(0,2) = 0.5*glstrain(5);
        gl(1,0) = gl(0,1);
        gl(1,1) = glstrain(1);
        gl(1,2) = 0.5*glstrain(4);
        gl(2,0) = gl(0,2);
        gl(2,1) = gl(1,2);
        gl(2,2) = glstrain(2);

        LINALG::SerialDenseMatrix temp(NUMDIM_SOTET4,NUMDIM_SOTET4);
        LINALG::SerialDenseMatrix euler_almansi(NUMDIM_SOTET4,NUMDIM_SOTET4);
        temp.Multiply('N','N',1.0,gl,f,0.0);
        euler_almansi.Multiply('T','N',1.0,f,temp,0.0);

        (*elestrain)(gp,0) = euler_almansi(0,0);
        (*elestrain)(gp,1) = euler_almansi(1,1);
        (*elestrain)(gp,2) = euler_almansi(2,2);
        (*elestrain)(gp,3) = euler_almansi(0,1);
        (*elestrain)(gp,4) = euler_almansi(1,2);
        (*elestrain)(gp,5) = euler_almansi(0,2);
      }
    }

    if (elestress)
    { // return 2nd Piola-Kirchhoff stresses
      if (!cauchy)
        for (int i = 0; i < NUMSTR_SOTET4; ++i)
          (*elestress)(gp,i) = stress(i);
      else
        for (int i = 0; i < NUMSTR_SOTET4; ++i)
          (*elestress)(gp,i) = cstress(i);
    }

    //-------------------------------------------- build operator Ypsilon
    LINALG::SerialDenseMatrix S(3,3);
    S(0,0) = stress(0);
    S(0,1) = stress(3);
    S(0,2) = stress(5);
    S(1,0) = stress(3);
    S(1,1) = stress(1);
    S(1,2) = stress(4);
    S(2,0) = stress(5);
    S(2,1) = stress(4);
    S(2,2) = stress(2);
    LINALG::SerialDenseMatrix Ypsilon(6,9);
    BuildYpsilon(Ypsilon,F,S);

    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * (tet4_dis.weights)(gp);

    //------------------------------------------- assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T','N',intfac,B,cstress,1.0);
    }

    //*******************************************************************
    if (stiffmatrix)
    {
      LINALG::SerialDenseMatrix sum(6,9);
      LINALG::SerialDenseMatrix tmp2(1,9);

      //================================== 1.) B^T * cstress * F^T * N_x
      sum.Multiply('N','T',1.0,cstress,FTvec,0.0);
      
      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      sum.Multiply('N','N',-j,Ypsilon,Lambda,1.0);
      
      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
      LINALG::SerialDenseMatrix sixnine(6,9);
      LINALG::SerialDenseMatrix sixsix(6,6);
      sixnine.Multiply('N','N',1.0,Theta,Lambda,0.0);
      sixsix.Multiply('N','N',1.0,IF,cmat,0.0);
      sum.Multiply('N','N',-j,sixsix,sixnine,1.0);

      //================ put everything together: K = intfac * B*T * sum * N_x
      LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOTET4);
      tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
      (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
    } // if (stiffmatrix)
    //*******************************************************************

  } // for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  // static integrator created in any case to safe "if-case"
  const static DRT::ELEMENTS::Integrator_tet4_4point tet4_mass;
  // evaluate mass matrix
  if (massmatrix != NULL)
  {
    //consistent mass matrix evaluated using a 4-point rule
    for (int gp=0; gp<tet4_mass.num_gp; gp++)
    {
      for (int inod=0; inod<NUMNOD_SOTET4; ++inod)
      {
        for (int jnod=0; jnod<NUMNOD_SOTET4; ++jnod)
        {
          double massfactor = (tet4_mass.shapefct_gp[gp])(inod) * density *
                              (tet4_mass.shapefct_gp[gp])(jnod) * detj *
                              (tet4_mass.weights)(gp);
          (*massmatrix)(NUMDIM_SOTET4*inod+0,NUMDIM_SOTET4*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+1,NUMDIM_SOTET4*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4*inod+2,NUMDIM_SOTET4*jnod+2) += massfactor;
        }
      }
    }
  }// end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
  

  return;
} // DRT::ELEMENTS::InvDesign::so_tet4_nlnstiffmass


/*----------------------------------------------------------------------*
 |  compute and store material configuration                  (private) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::sot4_StoreMaterialConfiguration(
                                       DRT::ELEMENTS::So_tet4* ele,
                                       const vector<double>& disp)
{
  const static DRT::ELEMENTS::Integrator_tet4_1point tet4_dis;
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOTET4,NUMDIM_SOTET4);
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOTET4,NUMDIM_SOTET4);
  for (int i=0; i<NUMNOD_SOTET4; ++i)
  {
    xcurr(i,0) = ele->Nodes()[i]->X()[0];
    xcurr(i,1) = ele->Nodes()[i]->X()[1];
    xcurr(i,2) = ele->Nodes()[i]->X()[2];

    xrefe(i,0) = xcurr(i,0) + disp[i*NODDOF_SOTET4+0];
    xrefe(i,1) = xcurr(i,1) + disp[i*NODDOF_SOTET4+1];
    xrefe(i,2) = xcurr(i,2) + disp[i*NODDOF_SOTET4+2];
  }

  /* get the matrix of the coordinates of nodes needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  **         J = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  **		 [ Z_1  Z_2  Z_3  Z_4 ]
  */
  LINALG::SerialDenseMatrix J(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
  for (int i=0; i<4; i++)  J(0,i)=1;
  for (int row=0;row<3;row++)
    for (int col=0;col<4;col++)
      J(row+1,col)= xrefe(col,row);
  // volume of the element
  detJ_[0] = LINALG::DeterminantLU(J)/6.0;

  for (int gp=0; gp<NUMGPT_SOTET4; ++gp)
  {
    LINALG::SerialDenseMatrix jac(NUMCOORD_SOTET4,NUMCOORD_SOTET4);
    LINALG::SerialDenseMatrix N_XYZ(NUMNOD_SOTET4,NUMDIM_SOTET4);
    LINALG::SerialDenseMatrix F(3,3);
  
    {
      LINALG::SerialDenseMatrix tmp(NUMCOORD_SOTET4-1,NUMCOORD_SOTET4);
      tmp.Multiply('T','N',1.0,xrefe,tet4_dis.deriv_gp[gp],0.0);
      for (int i=0; i<4; i++) jac(0,i)=1;
      for (int row=0;row<3;row++)
        for (int col=0;col<4;col++)
          jac(row+1,col)=tmp(row,col);
    }
    // size is 4x3
    Epetra_SerialDenseMatrix  I_aug(NUMCOORD_SOTET4,NUMDIM_SOTET4);
    // size is 4x3
    Epetra_SerialDenseMatrix partials(NUMCOORD_SOTET4,NUMDIM_SOTET4);
    I_aug(1,0)=1;
    I_aug(2,1)=1;
    I_aug(3,2)=1;

    Epetra_SerialDenseSolver solve_for_inverseJac;  // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);            // set A=jac
    solve_for_inverseJac.SetVectors(partials,I_aug);// set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();         // partials = jac^-1.I_aug
    if ((err != 0) && (err2!=0))
    	dserror("Inversion of Jacobian failed");

    N_XYZ.Multiply('N','N',1.0,tet4_dis.deriv_gp[gp],partials,0.0);
    F.Multiply('T','N',1.0,xcurr,N_XYZ,0.0);

    // put stuff into storage
    MatrixtoStorage(gp,N_XYZ,JHistory());
    MatrixtoStorage(gp,F,FHistory());

  } // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  return;
}


#endif  // #if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#endif  // #ifdef CCADISCRET
