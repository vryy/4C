/*----------------------------------------------------------------------*/
/*! \file
\brief special element adaptions for inverse design
\maintainer Christoph Meier
\level 2

*----------------------------------------------------------------------*/

#include "inversedesign.H"
#include "../drt_mat/material.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "so_tet4.H"
#include "so_weg6.H"
#include "../drt_lib/drt_node.H"


/*----------------------------------------------------------------------*
 |  Lambda tensor                                             (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildLambda(
    LINALG::Matrix<9, 9>& L, const LINALG::Matrix<3, 3>& F) const

{
#if 0  // this is for theory only
  double Lambda4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int p=0; p<3; ++p)
      for (int q=0; q<3; ++q)
        for (int m=0; m<3; ++m)
          Lambda4[k][m][p][q] = F(k,p)*F(q,m);
#endif


#if 1
  L(0, 0) = F(0, 0) * F(0, 0);
  L(0, 1) = F(0, 1) * F(0, 0);
  L(0, 2) = F(0, 2) * F(0, 0);
  L(0, 3) = F(0, 0) * F(1, 0);
  L(0, 4) = F(0, 1) * F(1, 0);
  L(0, 5) = F(0, 2) * F(1, 0);
  L(0, 6) = F(0, 0) * F(2, 0);
  L(0, 7) = F(0, 1) * F(2, 0);
  L(0, 8) = F(0, 2) * F(2, 0);

  L(1, 0) = F(1, 0) * F(0, 0);
  L(1, 1) = F(1, 1) * F(0, 0);
  L(1, 2) = F(1, 2) * F(0, 0);
  L(1, 3) = F(1, 0) * F(1, 0);
  L(1, 4) = F(1, 1) * F(1, 0);
  L(1, 5) = F(1, 2) * F(1, 0);
  L(1, 6) = F(1, 0) * F(2, 0);
  L(1, 7) = F(1, 1) * F(2, 0);
  L(1, 8) = F(1, 2) * F(2, 0);

  L(2, 0) = F(2, 0) * F(0, 0);
  L(2, 1) = F(2, 1) * F(0, 0);
  L(2, 2) = F(2, 2) * F(0, 0);
  L(2, 3) = F(2, 0) * F(1, 0);
  L(2, 4) = F(2, 1) * F(1, 0);
  L(2, 5) = F(2, 2) * F(1, 0);
  L(2, 6) = F(2, 0) * F(2, 0);
  L(2, 7) = F(2, 1) * F(2, 0);
  L(2, 8) = F(2, 2) * F(2, 0);

  L(3, 0) = F(0, 0) * F(0, 1);
  L(3, 1) = F(0, 1) * F(0, 1);
  L(3, 2) = F(0, 2) * F(0, 1);
  L(3, 3) = F(0, 0) * F(1, 1);
  L(3, 4) = F(0, 1) * F(1, 1);
  L(3, 5) = F(0, 2) * F(1, 1);
  L(3, 6) = F(0, 0) * F(2, 1);
  L(3, 7) = F(0, 1) * F(2, 1);
  L(3, 8) = F(0, 2) * F(2, 1);

  L(4, 0) = F(1, 0) * F(0, 1);
  L(4, 1) = F(1, 1) * F(0, 1);
  L(4, 2) = F(1, 2) * F(0, 1);
  L(4, 3) = F(1, 0) * F(1, 1);
  L(4, 4) = F(1, 1) * F(1, 1);
  L(4, 5) = F(1, 2) * F(1, 1);
  L(4, 6) = F(1, 0) * F(2, 1);
  L(4, 7) = F(1, 1) * F(2, 1);
  L(4, 8) = F(1, 2) * F(2, 1);

  L(5, 0) = F(2, 0) * F(0, 1);
  L(5, 1) = F(2, 1) * F(0, 1);
  L(5, 2) = F(2, 2) * F(0, 1);
  L(5, 3) = F(2, 0) * F(1, 1);
  L(5, 4) = F(2, 1) * F(1, 1);
  L(5, 5) = F(2, 2) * F(1, 1);
  L(5, 6) = F(2, 0) * F(2, 1);
  L(5, 7) = F(2, 1) * F(2, 1);
  L(5, 8) = F(2, 2) * F(2, 1);

  L(6, 0) = F(0, 0) * F(0, 2);
  L(6, 1) = F(0, 1) * F(0, 2);
  L(6, 2) = F(0, 2) * F(0, 2);
  L(6, 3) = F(0, 0) * F(1, 2);
  L(6, 4) = F(0, 1) * F(1, 2);
  L(6, 5) = F(0, 2) * F(1, 2);
  L(6, 6) = F(0, 0) * F(2, 2);
  L(6, 7) = F(0, 1) * F(2, 2);
  L(6, 8) = F(0, 2) * F(2, 2);

  L(7, 0) = F(1, 0) * F(0, 2);
  L(7, 1) = F(1, 1) * F(0, 2);
  L(7, 2) = F(1, 2) * F(0, 2);
  L(7, 3) = F(1, 0) * F(1, 2);
  L(7, 4) = F(1, 1) * F(1, 2);
  L(7, 5) = F(1, 2) * F(1, 2);
  L(7, 6) = F(1, 0) * F(2, 2);
  L(7, 7) = F(1, 1) * F(2, 2);
  L(7, 8) = F(1, 2) * F(2, 2);

  L(8, 0) = F(2, 0) * F(0, 2);
  L(8, 1) = F(2, 1) * F(0, 2);
  L(8, 2) = F(2, 2) * F(0, 2);
  L(8, 3) = F(2, 0) * F(1, 2);
  L(8, 4) = F(2, 1) * F(1, 2);
  L(8, 5) = F(2, 2) * F(1, 2);
  L(8, 6) = F(2, 0) * F(2, 2);
  L(8, 7) = F(2, 1) * F(2, 2);
  L(8, 8) = F(2, 2) * F(2, 2);
#else
  //               k  m  p  q
  L(0, 0) = Lambda4[0][0][0][0];
  L(0, 1) = Lambda4[0][0][1][0];
  L(0, 2) = Lambda4[0][0][2][0];
  L(0, 3) = Lambda4[0][0][0][1];
  L(0, 4) = Lambda4[0][0][1][1];
  L(0, 5) = Lambda4[0][0][2][1];
  L(0, 6) = Lambda4[0][0][0][2];
  L(0, 7) = Lambda4[0][0][1][2];
  L(0, 8) = Lambda4[0][0][2][2];

  L(1, 0) = Lambda4[1][0][0][0];
  L(1, 1) = Lambda4[1][0][1][0];
  L(1, 2) = Lambda4[1][0][2][0];
  L(1, 3) = Lambda4[1][0][0][1];
  L(1, 4) = Lambda4[1][0][1][1];
  L(1, 5) = Lambda4[1][0][2][1];
  L(1, 6) = Lambda4[1][0][0][2];
  L(1, 7) = Lambda4[1][0][1][2];
  L(1, 8) = Lambda4[1][0][2][2];

  L(2, 0) = Lambda4[2][0][0][0];
  L(2, 1) = Lambda4[2][0][1][0];
  L(2, 2) = Lambda4[2][0][2][0];
  L(2, 3) = Lambda4[2][0][0][1];
  L(2, 4) = Lambda4[2][0][1][1];
  L(2, 5) = Lambda4[2][0][2][1];
  L(2, 6) = Lambda4[2][0][0][2];
  L(2, 7) = Lambda4[2][0][1][2];
  L(2, 8) = Lambda4[2][0][2][2];

  L(3, 0) = Lambda4[0][1][0][0];
  L(3, 1) = Lambda4[0][1][1][0];
  L(3, 2) = Lambda4[0][1][2][0];
  L(3, 3) = Lambda4[0][1][0][1];
  L(3, 4) = Lambda4[0][1][1][1];
  L(3, 5) = Lambda4[0][1][2][1];
  L(3, 6) = Lambda4[0][1][0][2];
  L(3, 7) = Lambda4[0][1][1][2];
  L(3, 8) = Lambda4[0][1][2][2];

  L(4, 0) = Lambda4[1][1][0][0];
  L(4, 1) = Lambda4[1][1][1][0];
  L(4, 2) = Lambda4[1][1][2][0];
  L(4, 3) = Lambda4[1][1][0][1];
  L(4, 4) = Lambda4[1][1][1][1];
  L(4, 5) = Lambda4[1][1][2][1];
  L(4, 6) = Lambda4[1][1][0][2];
  L(4, 7) = Lambda4[1][1][1][2];
  L(4, 8) = Lambda4[1][1][2][2];

  L(5, 0) = Lambda4[2][1][0][0];
  L(5, 1) = Lambda4[2][1][1][0];
  L(5, 2) = Lambda4[2][1][2][0];
  L(5, 3) = Lambda4[2][1][0][1];
  L(5, 4) = Lambda4[2][1][1][1];
  L(5, 5) = Lambda4[2][1][2][1];
  L(5, 6) = Lambda4[2][1][0][2];
  L(5, 7) = Lambda4[2][1][1][2];
  L(5, 8) = Lambda4[2][1][2][2];

  L(6, 0) = Lambda4[0][2][0][0];
  L(6, 1) = Lambda4[0][2][1][0];
  L(6, 2) = Lambda4[0][2][2][0];
  L(6, 3) = Lambda4[0][2][0][1];
  L(6, 4) = Lambda4[0][2][1][1];
  L(6, 5) = Lambda4[0][2][2][1];
  L(6, 6) = Lambda4[0][2][0][2];
  L(6, 7) = Lambda4[0][2][1][2];
  L(6, 8) = Lambda4[0][2][2][2];

  L(7, 0) = Lambda4[1][2][0][0];
  L(7, 1) = Lambda4[1][2][1][0];
  L(7, 2) = Lambda4[1][2][2][0];
  L(7, 3) = Lambda4[1][2][0][1];
  L(7, 4) = Lambda4[1][2][1][1];
  L(7, 5) = Lambda4[1][2][2][1];
  L(7, 6) = Lambda4[1][2][0][2];
  L(7, 7) = Lambda4[1][2][1][2];
  L(7, 8) = Lambda4[1][2][2][2];

  L(8, 0) = Lambda4[2][2][0][0];
  L(8, 1) = Lambda4[2][2][1][0];
  L(8, 2) = Lambda4[2][2][2][0];
  L(8, 3) = Lambda4[2][2][0][1];
  L(8, 4) = Lambda4[2][2][1][1];
  L(8, 5) = Lambda4[2][2][2][1];
  L(8, 6) = Lambda4[2][2][0][2];
  L(8, 7) = Lambda4[2][2][1][2];
  L(8, 8) = Lambda4[2][2][2][2];
#endif

  return;
}  // DRT::ELEMENTS::InvDesign::BuildLambda

/*----------------------------------------------------------------------*
 |  IF tensor                                                 (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildIF(
    LINALG::Matrix<6, 6>& IF, const LINALG::Matrix<3, 3>& F) const

{
#if 0  // this is for theory reference only
  double IF4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int l=0; l<3; ++l)
      for (int m=0; m<3; ++m)
        for (int n=0; n<3; ++n)
          IF4[k][l][m][n] = 0.5 * ( F(k,m)*F(l,n) + F(k,n)*F(l,m) );
#endif

#if 1
  IF(0, 0) = (F(0, 0) * F(0, 0) + F(0, 0) * F(0, 0));
  IF(0, 1) = (F(0, 1) * F(0, 1) + F(0, 1) * F(0, 1));
  IF(0, 2) = (F(0, 2) * F(0, 2) + F(0, 2) * F(0, 2));
  IF(0, 3) = (F(0, 0) * F(0, 1) + F(0, 1) * F(0, 0)) * 2.0;
  IF(0, 4) = (F(0, 1) * F(0, 2) + F(0, 2) * F(0, 1)) * 2.0;
  IF(0, 5) = (F(0, 2) * F(0, 0) + F(0, 0) * F(0, 2)) * 2.0;

  IF(1, 0) = (F(1, 0) * F(1, 0) + F(1, 0) * F(1, 0));
  IF(1, 1) = (F(1, 1) * F(1, 1) + F(1, 1) * F(1, 1));
  IF(1, 2) = (F(1, 2) * F(1, 2) + F(1, 2) * F(1, 2));
  IF(1, 3) = (F(1, 0) * F(1, 1) + F(1, 1) * F(1, 0)) * 2.0;
  IF(1, 4) = (F(1, 1) * F(1, 2) + F(1, 2) * F(1, 1)) * 2.0;
  IF(1, 5) = (F(1, 2) * F(1, 0) + F(1, 0) * F(1, 2)) * 2.0;

  IF(2, 0) = (F(2, 0) * F(2, 0) + F(2, 0) * F(2, 0));
  IF(2, 1) = (F(2, 1) * F(2, 1) + F(2, 1) * F(2, 1));
  IF(2, 2) = (F(2, 2) * F(2, 2) + F(2, 2) * F(2, 2));
  IF(2, 3) = (F(2, 0) * F(2, 1) + F(2, 1) * F(2, 0)) * 2.0;
  IF(2, 4) = (F(2, 1) * F(2, 2) + F(2, 2) * F(2, 1)) * 2.0;
  IF(2, 5) = (F(2, 2) * F(2, 0) + F(2, 0) * F(2, 2)) * 2.0;

  IF(3, 0) = (F(0, 0) * F(1, 0) + F(0, 0) * F(1, 0));
  IF(3, 1) = (F(0, 1) * F(1, 1) + F(0, 1) * F(1, 1));
  IF(3, 2) = (F(0, 2) * F(1, 2) + F(0, 2) * F(1, 2));
  IF(3, 3) = (F(0, 0) * F(1, 1) + F(0, 1) * F(1, 0)) * 2.0;
  IF(3, 4) = (F(0, 1) * F(1, 2) + F(0, 2) * F(1, 1)) * 2.0;
  IF(3, 5) = (F(0, 2) * F(1, 0) + F(0, 0) * F(1, 2)) * 2.0;

  IF(4, 0) = (F(1, 0) * F(2, 0) + F(1, 0) * F(2, 0));
  IF(4, 1) = (F(1, 1) * F(2, 1) + F(1, 1) * F(2, 1));
  IF(4, 2) = (F(1, 2) * F(2, 2) + F(1, 2) * F(2, 2));
  IF(4, 3) = (F(1, 0) * F(2, 1) + F(1, 1) * F(2, 0)) * 2.0;
  IF(4, 4) = (F(1, 1) * F(2, 2) + F(1, 2) * F(2, 1)) * 2.0;
  IF(4, 5) = (F(1, 2) * F(2, 0) + F(1, 0) * F(2, 2)) * 2.0;

  IF(5, 0) = (F(2, 0) * F(0, 0) + F(2, 0) * F(0, 0));
  IF(5, 1) = (F(2, 1) * F(0, 1) + F(2, 1) * F(0, 1));
  IF(5, 2) = (F(2, 2) * F(0, 2) + F(2, 2) * F(0, 2));
  IF(5, 3) = (F(2, 0) * F(0, 1) + F(2, 1) * F(0, 0)) * 2.0;
  IF(5, 4) = (F(2, 1) * F(0, 2) + F(2, 2) * F(0, 1)) * 2.0;
  IF(5, 5) = (F(2, 2) * F(0, 0) + F(2, 0) * F(0, 2)) * 2.0;

  IF.Scale(0.5);

#else
  //            k  l  m  n
  IF(0, 0) = IF4[0][0][0][0];
  IF(0, 1) = IF4[0][0][1][1];
  IF(0, 2) = IF4[0][0][2][2];
  IF(0, 3) = IF4[0][0][0][1] * 2.0;
  IF(0, 4) = IF4[0][0][1][2] * 2.0;
  IF(0, 5) = IF4[0][0][2][0] * 2.0;

  IF(1, 0) = IF4[1][1][0][0];
  IF(1, 1) = IF4[1][1][1][1];
  IF(1, 2) = IF4[1][1][2][2];
  IF(1, 3) = IF4[1][1][0][1] * 2.0;
  IF(1, 4) = IF4[1][1][1][2] * 2.0;
  IF(1, 5) = IF4[1][1][2][0] * 2.0;

  IF(2, 0) = IF4[2][2][0][0];
  IF(2, 1) = IF4[2][2][1][1];
  IF(2, 2) = IF4[2][2][2][2];
  IF(2, 3) = IF4[2][2][0][1] * 2.0;
  IF(2, 4) = IF4[2][2][1][2] * 2.0;
  IF(2, 5) = IF4[2][2][2][0] * 2.0;

  IF(3, 0) = IF4[0][1][0][0];
  IF(3, 1) = IF4[0][1][1][1];
  IF(3, 2) = IF4[0][1][2][2];
  IF(3, 3) = IF4[0][1][0][1] * 2.0;
  IF(3, 4) = IF4[0][1][1][2] * 2.0;
  IF(3, 5) = IF4[0][1][2][0] * 2.0;

  IF(4, 0) = IF4[1][2][0][0];
  IF(4, 1) = IF4[1][2][1][1];
  IF(4, 2) = IF4[1][2][2][2];
  IF(4, 3) = IF4[1][2][0][1] * 2.0;
  IF(4, 4) = IF4[1][2][1][2] * 2.0;
  IF(4, 5) = IF4[1][2][2][0] * 2.0;

  IF(5, 0) = IF4[2][0][0][0];
  IF(5, 1) = IF4[2][0][1][1];
  IF(5, 2) = IF4[2][0][2][2];
  IF(5, 3) = IF4[2][0][0][1] * 2.0;
  IF(5, 4) = IF4[2][0][1][2] * 2.0;
  IF(5, 5) = IF4[2][0][2][0] * 2.0;
#endif

  return;
}  // DRT::ELEMENTS::InvDesign::BuildIF

/*----------------------------------------------------------------------*
 |  Theta tensor                                              (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildTheta(
    LINALG::Matrix<6, 9>& Theta, const LINALG::Matrix<3, 3>& F) const

{
#if 0  // for theory reference only
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
  Theta(0, 0) = (F(0, 0));
  Theta(0, 1) = (F(1, 0));
  Theta(0, 2) = (F(2, 0));
  Theta(0, 3) = 0.0;
  Theta(0, 4) = 0.0;
  Theta(0, 5) = 0.0;
  Theta(0, 6) = 0.0;
  Theta(0, 7) = 0.0;
  Theta(0, 8) = 0.0;

  Theta(1, 0) = 0.0;
  Theta(1, 1) = 0.0;
  Theta(1, 2) = 0.0;
  Theta(1, 3) = (F(0, 1));
  Theta(1, 4) = (F(1, 1));
  Theta(1, 5) = (F(2, 1));
  Theta(1, 6) = 0.0;
  Theta(1, 7) = 0.0;
  Theta(1, 8) = 0.0;

  Theta(2, 0) = 0.0;
  Theta(2, 1) = 0.0;
  Theta(2, 2) = 0.0;
  Theta(2, 3) = 0.0;
  Theta(2, 4) = 0.0;
  Theta(2, 5) = 0.0;
  Theta(2, 6) = (F(0, 2));
  Theta(2, 7) = (F(1, 2));
  Theta(2, 8) = (F(2, 2));

  Theta(3, 0) = (F(0, 1));
  Theta(3, 1) = (F(1, 1));
  Theta(3, 2) = (F(2, 1));
  Theta(3, 3) = (F(0, 0));
  Theta(3, 4) = (F(1, 0));
  Theta(3, 5) = (F(2, 0));
  Theta(3, 6) = 0.0;
  Theta(3, 7) = 0.0;
  Theta(3, 8) = 0.0;

  Theta(4, 0) = 0.0;
  Theta(4, 1) = 0.0;
  Theta(4, 2) = 0.0;
  Theta(4, 3) = (F(0, 2));
  Theta(4, 4) = (F(1, 2));
  Theta(4, 5) = (F(2, 2));
  Theta(4, 6) = (F(0, 1));
  Theta(4, 7) = (F(1, 1));
  Theta(4, 8) = (F(2, 1));

  Theta(5, 0) = (F(0, 2));
  Theta(5, 1) = (F(1, 2));
  Theta(5, 2) = (F(2, 2));
  Theta(5, 3) = 0.0;
  Theta(5, 4) = 0.0;
  Theta(5, 5) = 0.0;
  Theta(5, 6) = (F(0, 0));
  Theta(5, 7) = (F(1, 0));
  Theta(5, 8) = (F(2, 0));
#else
  //                  i  j  k  l
  Theta(0, 0) = Theta4[0][0][0][0];
  Theta(0, 1) = Theta4[0][0][1][0];
  Theta(0, 2) = Theta4[0][0][2][0];
  Theta(0, 3) = Theta4[0][0][0][1];
  Theta(0, 4) = Theta4[0][0][1][1];
  Theta(0, 5) = Theta4[0][0][2][1];
  Theta(0, 6) = Theta4[0][0][0][2];
  Theta(0, 7) = Theta4[0][0][1][2];
  Theta(0, 8) = Theta4[0][0][2][2];

  Theta(1, 0) = Theta4[1][1][0][0];
  Theta(1, 1) = Theta4[1][1][1][0];
  Theta(1, 2) = Theta4[1][1][2][0];
  Theta(1, 3) = Theta4[1][1][0][1];
  Theta(1, 4) = Theta4[1][1][1][1];
  Theta(1, 5) = Theta4[1][1][2][1];
  Theta(1, 6) = Theta4[1][1][0][2];
  Theta(1, 7) = Theta4[1][1][1][2];
  Theta(1, 8) = Theta4[1][1][2][2];

  Theta(2, 0) = Theta4[2][2][0][0];
  Theta(2, 1) = Theta4[2][2][1][0];
  Theta(2, 2) = Theta4[2][2][2][0];
  Theta(2, 3) = Theta4[2][2][0][1];
  Theta(2, 4) = Theta4[2][2][1][1];
  Theta(2, 5) = Theta4[2][2][2][1];
  Theta(2, 6) = Theta4[2][2][0][2];
  Theta(2, 7) = Theta4[2][2][1][2];
  Theta(2, 8) = Theta4[2][2][2][2];

  Theta(3, 0) = Theta4[0][1][0][0] * 2.0;
  Theta(3, 1) = Theta4[0][1][1][0] * 2.0;
  Theta(3, 2) = Theta4[0][1][2][0] * 2.0;
  Theta(3, 3) = Theta4[0][1][0][1] * 2.0;
  Theta(3, 4) = Theta4[0][1][1][1] * 2.0;
  Theta(3, 5) = Theta4[0][1][2][1] * 2.0;
  Theta(3, 6) = Theta4[0][1][0][2] * 2.0;
  Theta(3, 7) = Theta4[0][1][1][2] * 2.0;
  Theta(3, 8) = Theta4[0][1][2][2] * 2.0;

  Theta(4, 0) = Theta4[1][2][0][0] * 2.0;
  Theta(4, 1) = Theta4[1][2][1][0] * 2.0;
  Theta(4, 2) = Theta4[1][2][2][0] * 2.0;
  Theta(4, 3) = Theta4[1][2][0][1] * 2.0;
  Theta(4, 4) = Theta4[1][2][1][1] * 2.0;
  Theta(4, 5) = Theta4[1][2][2][1] * 2.0;
  Theta(4, 6) = Theta4[1][2][0][2] * 2.0;
  Theta(4, 7) = Theta4[1][2][1][2] * 2.0;
  Theta(4, 8) = Theta4[1][2][2][2] * 2.0;

  Theta(5, 0) = Theta4[2][0][0][0] * 2.0;
  Theta(5, 1) = Theta4[2][0][1][0] * 2.0;
  Theta(5, 2) = Theta4[2][0][2][0] * 2.0;
  Theta(5, 3) = Theta4[2][0][0][1] * 2.0;
  Theta(5, 4) = Theta4[2][0][1][1] * 2.0;
  Theta(5, 5) = Theta4[2][0][2][1] * 2.0;
  Theta(5, 6) = Theta4[2][0][0][2] * 2.0;
  Theta(5, 7) = Theta4[2][0][1][2] * 2.0;
  Theta(5, 8) = Theta4[2][0][2][2] * 2.0;
#endif

  return;
}  // DRT::ELEMENTS::InvDesign::BuildTheta


/*----------------------------------------------------------------------*
 |  Ypsilon tensor                                            (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildYpsilon(
    LINALG::Matrix<6, 9>& Y, const LINALG::Matrix<3, 3>& F, const LINALG::Matrix<3, 3>& S) const

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
  Y.Clear();

  Y(0, 0) = (4.0 * F(0, 0)) * S(0, 0) + (2.0 * F(0, 1)) * S(0, 1) + (2.0 * F(0, 2)) * S(0, 2) +
            (2.0 * F(0, 1)) * S(1, 0) + (2.0 * F(0, 2)) * S(2, 0);
  Y(0, 3) = +(2.0 * F(0, 0)) * S(0, 1) + (2.0 * F(0, 0)) * S(1, 0) + (4.0 * F(0, 1)) * S(1, 1) +
            (2.0 * F(0, 2)) * S(1, 2) + (2.0 * F(0, 2)) * S(2, 1);
  Y(0, 6) = +(2.0 * F(0, 0)) * S(0, 2) + (2.0 * F(0, 1)) * S(1, 2) + (2.0 * F(0, 0)) * S(2, 0) +
            (2.0 * F(0, 1)) * S(2, 1) + (4.0 * F(0, 2)) * S(2, 2);
  Y(1, 1) = (4.0 * F(1, 0)) * S(0, 0) + (2.0 * F(1, 1)) * S(0, 1) + (2.0 * F(1, 2)) * S(0, 2) +
            (2.0 * F(1, 1)) * S(1, 0) + (2.0 * F(1, 2)) * S(2, 0);
  Y(1, 4) = +(2.0 * F(1, 0)) * S(0, 1) + (2.0 * F(1, 0)) * S(1, 0) + (4.0 * F(1, 1)) * S(1, 1) +
            (2.0 * F(1, 2)) * S(1, 2) + (2.0 * F(1, 2)) * S(2, 1);
  Y(1, 7) = +(2.0 * F(1, 0)) * S(0, 2) + (2.0 * F(1, 1)) * S(1, 2) + (2.0 * F(1, 0)) * S(2, 0) +
            (2.0 * F(1, 1)) * S(2, 1) + (4.0 * F(1, 2)) * S(2, 2);
  Y(2, 2) = (4.0 * F(2, 0)) * S(0, 0) + (2.0 * F(2, 1)) * S(0, 1) + (2.0 * F(2, 2)) * S(0, 2) +
            (2.0 * F(2, 1)) * S(1, 0) + (2.0 * F(2, 2)) * S(2, 0);
  Y(2, 5) = +(2.0 * F(2, 0)) * S(0, 1) + (2.0 * F(2, 0)) * S(1, 0) + (4.0 * F(2, 1)) * S(1, 1) +
            (2.0 * F(2, 2)) * S(1, 2) + (2.0 * F(2, 2)) * S(2, 1);
  Y(2, 8) = +(2.0 * F(2, 0)) * S(0, 2) + (2.0 * F(2, 1)) * S(1, 2) + (2.0 * F(2, 0)) * S(2, 0) +
            (2.0 * F(2, 1)) * S(2, 1) + (4.0 * F(2, 2)) * S(2, 2);
  Y(3, 0) = (2.0 * F(1, 0)) * S(0, 0) + F(1, 1) * S(0, 1) + F(1, 2) * S(0, 2) + F(1, 1) * S(1, 0) +
            F(1, 2) * S(2, 0);
  Y(3, 1) = (2.0 * F(0, 0)) * S(0, 0) + F(0, 1) * S(0, 1) + F(0, 2) * S(0, 2) + F(0, 1) * S(1, 0) +
            F(0, 2) * S(2, 0);
  Y(3, 3) = +F(1, 0) * S(0, 1) + F(1, 0) * S(1, 0) + (2.0 * F(1, 1)) * S(1, 1) + F(1, 2) * S(1, 2) +
            F(1, 2) * S(2, 1);
  Y(3, 4) = +F(0, 0) * S(0, 1) + F(0, 0) * S(1, 0) + (2.0 * F(0, 1)) * S(1, 1) + F(0, 2) * S(1, 2) +
            F(0, 2) * S(2, 1);
  Y(3, 6) = +F(1, 0) * S(0, 2) + F(1, 1) * S(1, 2) + F(1, 0) * S(2, 0) + F(1, 1) * S(2, 1) +
            (2.0 * F(1, 2)) * S(2, 2);
  Y(3, 7) = +F(0, 0) * S(0, 2) + F(0, 1) * S(1, 2) + F(0, 0) * S(2, 0) + F(0, 1) * S(2, 1) +
            (2.0 * F(0, 2)) * S(2, 2);
  Y(4, 1) = (2.0 * F(2, 0)) * S(0, 0) + F(2, 1) * S(0, 1) + F(2, 2) * S(0, 2) + F(2, 1) * S(1, 0) +
            F(2, 2) * S(2, 0);
  Y(4, 2) = (2.0 * F(1, 0)) * S(0, 0) + F(1, 1) * S(0, 1) + F(1, 2) * S(0, 2) + F(1, 1) * S(1, 0) +
            F(1, 2) * S(2, 0);
  Y(4, 4) = +F(2, 0) * S(0, 1) + F(2, 0) * S(1, 0) + (2.0 * F(2, 1)) * S(1, 1) + F(2, 2) * S(1, 2) +
            F(2, 2) * S(2, 1);
  Y(4, 5) = +F(1, 0) * S(0, 1) + F(1, 0) * S(1, 0) + (2.0 * F(1, 1)) * S(1, 1) + F(1, 2) * S(1, 2) +
            F(1, 2) * S(2, 1);
  Y(4, 7) = +F(2, 0) * S(0, 2) + F(2, 1) * S(1, 2) + F(2, 0) * S(2, 0) + F(2, 1) * S(2, 1) +
            (2.0 * F(2, 2)) * S(2, 2);
  Y(4, 8) = +F(1, 0) * S(0, 2) + F(1, 1) * S(1, 2) + F(1, 0) * S(2, 0) + F(1, 1) * S(2, 1) +
            (2.0 * F(1, 2)) * S(2, 2);
  Y(5, 0) = (2.0 * F(2, 0)) * S(0, 0) + F(2, 1) * S(0, 1) + F(2, 2) * S(0, 2) + F(2, 1) * S(1, 0) +
            F(2, 2) * S(2, 0);
  Y(5, 2) = (2.0 * F(0, 0)) * S(0, 0) + F(0, 1) * S(0, 1) + F(0, 2) * S(0, 2) + F(0, 1) * S(1, 0) +
            F(0, 2) * S(2, 0);
  Y(5, 3) = +F(2, 0) * S(0, 1) + F(2, 0) * S(1, 0) + (2.0 * F(2, 1)) * S(1, 1) + F(2, 2) * S(1, 2) +
            F(2, 2) * S(2, 1);
  Y(5, 5) = +F(0, 0) * S(0, 1) + F(0, 0) * S(1, 0) + (2.0 * F(0, 1)) * S(1, 1) + F(0, 2) * S(1, 2) +
            F(0, 2) * S(2, 1);
  Y(5, 6) = +F(2, 0) * S(0, 2) + F(2, 1) * S(1, 2) + F(2, 0) * S(2, 0) + F(2, 1) * S(2, 1) +
            (2.0 * F(2, 2)) * S(2, 2);
  Y(5, 8) = +F(0, 0) * S(0, 2) + F(0, 1) * S(1, 2) + F(0, 0) * S(2, 0) + F(0, 1) * S(2, 1) +
            (2.0 * F(0, 2)) * S(2, 2);
  Y.Scale(0.5);

#else
  //          k  l  p  q
  Y(0, 0) = Y4[0][0][0][0];
  Y(0, 1) = Y4[0][0][1][0];
  Y(0, 2) = Y4[0][0][2][0];
  Y(0, 3) = Y4[0][0][0][1];
  Y(0, 4) = Y4[0][0][1][1];
  Y(0, 5) = Y4[0][0][2][1];
  Y(0, 6) = Y4[0][0][0][2];
  Y(0, 7) = Y4[0][0][1][2];
  Y(0, 8) = Y4[0][0][2][2];

  Y(1, 0) = Y4[1][1][0][0];
  Y(1, 1) = Y4[1][1][1][0];
  Y(1, 2) = Y4[1][1][2][0];
  Y(1, 3) = Y4[1][1][0][1];
  Y(1, 4) = Y4[1][1][1][1];
  Y(1, 5) = Y4[1][1][2][1];
  Y(1, 6) = Y4[1][1][0][2];
  Y(1, 7) = Y4[1][1][1][2];
  Y(1, 8) = Y4[1][1][2][2];

  Y(2, 0) = Y4[2][2][0][0];
  Y(2, 1) = Y4[2][2][1][0];
  Y(2, 2) = Y4[2][2][2][0];
  Y(2, 3) = Y4[2][2][0][1];
  Y(2, 4) = Y4[2][2][1][1];
  Y(2, 5) = Y4[2][2][2][1];
  Y(2, 6) = Y4[2][2][0][2];
  Y(2, 7) = Y4[2][2][1][2];
  Y(2, 8) = Y4[2][2][2][2];

  Y(3, 0) = Y4[0][1][0][0];
  Y(3, 1) = Y4[0][1][1][0];
  Y(3, 2) = Y4[0][1][2][0];
  Y(3, 3) = Y4[0][1][0][1];
  Y(3, 4) = Y4[0][1][1][1];
  Y(3, 5) = Y4[0][1][2][1];
  Y(3, 6) = Y4[0][1][0][2];
  Y(3, 7) = Y4[0][1][1][2];
  Y(3, 8) = Y4[0][1][2][2];

  Y(4, 0) = Y4[1][2][0][0];
  Y(4, 1) = Y4[1][2][1][0];
  Y(4, 2) = Y4[1][2][2][0];
  Y(4, 3) = Y4[1][2][0][1];
  Y(4, 4) = Y4[1][2][1][1];
  Y(4, 5) = Y4[1][2][2][1];
  Y(4, 6) = Y4[1][2][0][2];
  Y(4, 7) = Y4[1][2][1][2];
  Y(4, 8) = Y4[1][2][2][2];

  Y(5, 0) = Y4[2][0][0][0];
  Y(5, 1) = Y4[2][0][1][0];
  Y(5, 2) = Y4[2][0][2][0];
  Y(5, 3) = Y4[2][0][0][1];
  Y(5, 4) = Y4[2][0][1][1];
  Y(5, 5) = Y4[2][0][2][1];
  Y(5, 6) = Y4[2][0][0][2];
  Y(5, 7) = Y4[2][0][1][2];
  Y(5, 8) = Y4[2][0][2][2];
#endif

  return;
}  // DRT::ELEMENTS::InvDesign::BuildYpsilon


/*----------------------------------------------------------------------*
 |  hex8 integration method                                    (public) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::soh8_nlnstiffmass(DRT::ELEMENTS::So_hex8* ele,  ///< this element
    std::vector<int>& lm,                                                      ///< location matrix
    std::vector<double>& disp,                                   ///< current displacements
    std::vector<double>& residual,                               ///< current residuum
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* stiffmatrix,       ///< element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH8, NUMDOF_SOH8>* massmatrix,        ///< element mass matrix
    LINALG::Matrix<NUMDOF_SOH8, 1>* force,                       ///< element internal force vector
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_SOH8, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  ///< stress output option
    const INPAR::STR::StrainType iostrain)  ///< strain output option
{
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8, 1>> shapefcts = ele->soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = ele->soh8_derivs();
  const static std::vector<double> gpweights = ele->soh8_weights();

  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8, NUMDIM_SOH8> xcurr;  // current  coord. of element
  DRT::Node** nodes = ele->Nodes();
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xcurr(i, 0) = x[0];
    xcurr(i, 1) = x[1];
    xcurr(i, 2) = x[2];

    xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_SOH8 + 0];
    xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_SOH8 + 1];
    xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_SOH8 + 2];
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    //----------------------------------- get inverse of Jacobian mapping
    const double detj = ele->detJ_[gp];
    LINALG::Matrix<3, 3>& invj = ele->invJ_[gp];

    //------------------------------ compute derivs wrt to spatial coords
    LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8> n_xyz;
    n_xyz.Multiply(invj, derivs[gp]);

    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> f;
    f.MultiplyTT(xrefe, n_xyz);

    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> F(false);
    const double j = F.Invert(f);

    //------------------------------------ build F^T as vector 9x1
    LINALG::Matrix<9, 1> FTvec;
    FTvec(0) = F(0, 0);
    FTvec(1) = F(0, 1);
    FTvec(2) = F(0, 2);
    FTvec(3) = F(1, 0);
    FTvec(4) = F(1, 1);
    FTvec(5) = F(1, 2);
    FTvec(6) = F(2, 0);
    FTvec(7) = F(2, 1);
    FTvec(8) = F(2, 2);

    //--------------------------------------------- build operator Lambda
    LINALG::Matrix<9, 9> Lambda;
    BuildLambda(Lambda, F);

    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::Matrix<6, 6> IF;
    BuildIF(IF, F);

    //---------------------------------------------- build operator Theta
    LINALG::Matrix<6, 9> Theta;
    BuildTheta(Theta, F);

    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> cauchygreen;
    cauchygreen.MultiplyTN(F, F);

    //-- Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOH8> B;
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      B(0, NODDOF_SOH8 * i + 0) = n_xyz(0, i);
      B(0, NODDOF_SOH8 * i + 1) = 0.0;
      B(0, NODDOF_SOH8 * i + 2) = 0.0;

      B(1, NODDOF_SOH8 * i + 0) = 0.0;
      B(1, NODDOF_SOH8 * i + 1) = n_xyz(1, i);
      B(1, NODDOF_SOH8 * i + 2) = 0.0;

      B(2, NODDOF_SOH8 * i + 0) = 0.0;
      B(2, NODDOF_SOH8 * i + 1) = 0.0;
      B(2, NODDOF_SOH8 * i + 2) = n_xyz(2, i);

      B(3, NODDOF_SOH8 * i + 0) = n_xyz(1, i);
      B(3, NODDOF_SOH8 * i + 1) = n_xyz(0, i);
      B(3, NODDOF_SOH8 * i + 2) = 0.0;

      B(4, NODDOF_SOH8 * i + 0) = 0.0;
      B(4, NODDOF_SOH8 * i + 1) = n_xyz(2, i);
      B(4, NODDOF_SOH8 * i + 2) = n_xyz(1, i);

      B(5, NODDOF_SOH8 * i + 0) = n_xyz(2, i);
      B(5, NODDOF_SOH8 * i + 1) = 0.0;
      B(5, NODDOF_SOH8 * i + 2) = n_xyz(0, i);
    }

    //--------------------------- build N_x operator (wrt spatial config)
    LINALG::Matrix<9, NUMDOF_SOH8> N_x(true);  // set to zero
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      N_x(0, 3 * i + 0) = n_xyz(0, i);
      N_x(1, 3 * i + 1) = n_xyz(0, i);
      N_x(2, 3 * i + 2) = n_xyz(0, i);

      N_x(3, 3 * i + 0) = n_xyz(1, i);
      N_x(4, 3 * i + 1) = n_xyz(1, i);
      N_x(5, 3 * i + 2) = n_xyz(1, i);

      N_x(6, 3 * i + 0) = n_xyz(2, i);
      N_x(7, 3 * i + 1) = n_xyz(2, i);
      N_x(8, 3 * i + 2) = n_xyz(2, i);
    }

    //------------------------------------------------- call material law
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
    params.set<int>("gp", gp);
    Teuchos::RCP<MAT::So3Material> so3mat =
        Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
    so3mat->Evaluate(&F, &glstrain, params, &stress, &cmat, gp, ele->Id());

    //------------------------------------------- compute cauchy stresses
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> cstress;
    cstress.Multiply(j, IF, stress);

    //--------------------------------------- output strains and stresses
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == NULL) dserror("strain data not available");
        // rewriting Green-Lagrange strains in matrix format
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> gl;
        gl(0, 0) = glstrain(0);
        gl(0, 1) = 0.5 * glstrain(3);
        gl(0, 2) = 0.5 * glstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = glstrain(1);
        gl(1, 2) = 0.5 * glstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = glstrain(2);

        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> temp;
        LINALG::Matrix<NUMDIM_SOH8, NUMDIM_SOH8> euler_almansi;
        temp.Multiply(gl, f);
        euler_almansi.MultiplyTN(f, temp);

        (*elestrain)(gp, 0) = euler_almansi(0, 0);
        (*elestrain)(gp, 1) = euler_almansi(1, 1);
        (*elestrain)(gp, 2) = euler_almansi(2, 2);
        (*elestrain)(gp, 3) = euler_almansi(0, 1);
        (*elestrain)(gp, 4) = euler_almansi(1, 2);
        (*elestrain)(gp, 5) = euler_almansi(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
    }

    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < 6; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < 6; ++i) (*elestress)(gp, i) = cstress(i);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
    }

    //-------------------------------------------- build operator Ypsilon
    LINALG::Matrix<3, 3> S;
    S(0, 0) = stress(0);
    S(0, 1) = stress(3);
    S(0, 2) = stress(5);
    S(1, 0) = stress(3);
    S(1, 1) = stress(1);
    S(1, 2) = stress(4);
    S(2, 0) = stress(5);
    S(2, 1) = stress(4);
    S(2, 2) = stress(2);
    LINALG::Matrix<6, 9> Ypsilon;
    BuildYpsilon(Ypsilon, F, S);

    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * gpweights[gp];

    //------------------------------------------- assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).MultiplyTN(intfac, B, cstress, 1.0);
    }

    //*******************************************************************
    if (stiffmatrix)
    {
      LINALG::Matrix<6, 9> sum;

      //================================== 1.) B^T * cstress * F^T * N_x
      sum.MultiplyNT(cstress, FTvec);

      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      sum.Multiply(-j, Ypsilon, Lambda, 1.0);

      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
      LINALG::Matrix<6, 9> sixnine;
      LINALG::Matrix<6, 6> sixsix;
      sixnine.Multiply(Theta, Lambda);
      sixsix.Multiply(IF, cmat);
      sum.Multiply(-j, sixsix, sixnine, 1.0);

      //================ put everything together: K = intfac * B*T * sum * N_x
      LINALG::Matrix<6, NUMDOF_SOH8> tmp3;
      tmp3.Multiply(sum, N_x);
      (*stiffmatrix).MultiplyTN(intfac, B, tmp3, 1.0);

    }  // if (stiffmatrix)
    //*******************************************************************


    //*******************************************************************
    // Strictly, inverse design analysis is stationary and should not have
    // a mass term. Loosely, if no Dirichlet-BCs are present a small
    // mass term might be used. Note that we use density per unit deformed
    // volume here!
    if (massmatrix)
    {
      double density = ele->Material()->Density(gp);
      const double fac = density * detj * gpweights[gp];
      for (int inod = 0; inod < NUMNOD_SOH8; ++inod)
        for (int jnod = 0; jnod < NUMNOD_SOH8; ++jnod)
        {
          const double massfactor = (shapefcts[gp])(inod) * (shapefcts[gp])(jnod)*fac;
          (*massmatrix)(NUMDIM_SOH8 * inod + 0, NUMDIM_SOH8 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 1, NUMDIM_SOH8 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8 * inod + 2, NUMDIM_SOH8 * jnod + 2) += massfactor;
        }
    }  // if (massmatrix)
    //*******************************************************************

  }  // for (unsigned gp=0; gp<NUMGPT_SOH8; ++gp)
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------


}  // DRT::ELEMENTS::InvDesign::soh8_nlnstiffmass



/*----------------------------------------------------------------------*
 |  compute and store material configuration                  (private) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::soh8_StoreMaterialConfiguration(
    DRT::ELEMENTS::So_hex8* ele, const std::vector<double>& disp)
{
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8, NUMNOD_SOH8>> derivs = ele->soh8_derivs();
  LINALG::Matrix<NUMNOD_SOH8, 3> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8, 3> xcurr;  // current  coord. of element
  DRT::Node** nodes = ele->Nodes();
  for (int i = 0; i < NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xcurr(i, 0) = x[0];
    xcurr(i, 1) = x[1];
    xcurr(i, 2) = x[2];

    xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_SOH8 + 0];
    xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_SOH8 + 1];
    xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_SOH8 + 2];
  }

  LINALG::Matrix<3, 3> invJ;
  LINALG::Matrix<3, NUMNOD_SOH8> N_XYZ;
  LINALG::Matrix<3, 3> F;
  for (unsigned gp = 0; gp < NUMGPT_SOH8; ++gp)
  {
    // compute invJ and detJ
    invJ.Multiply(derivs[gp], xrefe);
    detJ_[gp] = invJ.Invert();

    // compute F
    N_XYZ.Multiply(invJ, derivs[gp]);
    F.MultiplyTT(xcurr, N_XYZ);

    // put stuff into storage
    MatrixtoStorage(gp, invJ, JHistory());
    MatrixtoStorage(gp, F, FHistory());
  }

  return;
}


/*----------------------------------------------------------------------*
 |  weg6 integration method                                    (public) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::sow6_nlnstiffmass(DRT::ELEMENTS::So_weg6* ele,  ///< this element
    std::vector<int>& lm,                                                      ///< location matrix
    std::vector<double>& disp,                                   ///< current displacements
    std::vector<double>& residual,                               ///< current residuum
    LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* stiffmatrix,       ///< element stiffness matrix
    LINALG::Matrix<NUMDOF_WEG6, NUMDOF_WEG6>* massmatrix,        ///< element mass matrix
    LINALG::Matrix<NUMDOF_WEG6, 1>* force,                       ///< element internal force vector
    LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_WEG6, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    Teuchos::ParameterList& params,         ///< algorithmic parameters e.g. time
    const INPAR::STR::StressType iostress,  ///< stress output option
    const INPAR::STR::StrainType iostrain)  ///< strain output option
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
  ** ============================================================================*/
  /* pointer to (static) shape function array
   * for each node, evaluated at each gp*/
  LINALG::Matrix<NUMNOD_WEG6, NUMGPT_WEG6>* shapefct;
  /* pointer to (static) shape function derivatives array
   * for each node wrt to each direction, evaluated at each gp*/
  LINALG::Matrix<NUMGPT_WEG6 * NUMDIM_WEG6, NUMNOD_WEG6>* deriv;
  /* pointer to (static) weight factors at each gp */
  LINALG::Matrix<NUMGPT_WEG6, 1>* weights;
  ele->sow6_shapederiv(&shapefct, &deriv, &weights);  // call to evaluate
  /* ============================================================================*/

  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::Matrix<NUMNOD_WEG6, NUMDIM_WEG6> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_WEG6, NUMDIM_WEG6> xcurr;  // current  coord. of element
  DRT::Node** nodes = ele->Nodes();
  for (int i = 0; i < NUMNOD_WEG6; ++i)
  {
    const double* x = nodes[i]->X();
    xcurr(i, 0) = x[0];
    xcurr(i, 1) = x[1];
    xcurr(i, 2) = x[2];

    xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_WEG6 + 0];
    xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_WEG6 + 1];
    xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_WEG6 + 2];
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  for (int gp = 0; gp < NUMGPT_WEG6; ++gp)
  {
    // get submatrix of deriv at actual gp
    LINALG::Matrix<NUMDIM_WEG6, NUMGPT_WEG6> deriv_gp;
    for (int m = 0; m < NUMDIM_WEG6; ++m)
      for (int n = 0; n < NUMGPT_WEG6; ++n) deriv_gp(m, n) = (*deriv)(NUMDIM_WEG6 * gp + m, n);

    const double detj = ele->detJ_[gp];
    LINALG::Matrix<3, 3>& invj = ele->invJ_[gp];

    //------------------------------ compute derivs wrt to spatial coords
    LINALG::Matrix<NUMDIM_WEG6, NUMNOD_WEG6> n_xyz;
    n_xyz.Multiply(invj, deriv_gp);

    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> f;
    f.MultiplyTT(xrefe, n_xyz);

    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::Matrix<3, 3> F(false);
    const double j = F.Invert(f);

    //------------------------------------ build F^T as vector 9x1
    LINALG::Matrix<9, 1> FTvec;
    FTvec(0) = F(0, 0);
    FTvec(1) = F(0, 1);
    FTvec(2) = F(0, 2);
    FTvec(3) = F(1, 0);
    FTvec(4) = F(1, 1);
    FTvec(5) = F(1, 2);
    FTvec(6) = F(2, 0);
    FTvec(7) = F(2, 1);
    FTvec(8) = F(2, 2);

    //--------------------------------------------- build operator Lambda
    LINALG::Matrix<9, 9> Lambda;
    BuildLambda(Lambda, F);

    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::Matrix<6, 6> IF;
    BuildIF(IF, F);

    //---------------------------------------------- build operator Theta
    LINALG::Matrix<6, 9> Theta;
    BuildTheta(Theta, F);

    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> cauchygreen;
    cauchygreen.MultiplyTN(F, F);

    //-- Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_WEG6> B;
    for (int i = 0; i < NUMNOD_WEG6; ++i)
    {
      B(0, NODDOF_WEG6 * i + 0) = n_xyz(0, i);
      B(0, NODDOF_WEG6 * i + 1) = 0.0;
      B(0, NODDOF_WEG6 * i + 2) = 0.0;

      B(1, NODDOF_WEG6 * i + 0) = 0.0;
      B(1, NODDOF_WEG6 * i + 1) = n_xyz(1, i);
      B(1, NODDOF_WEG6 * i + 2) = 0.0;

      B(2, NODDOF_WEG6 * i + 0) = 0.0;
      B(2, NODDOF_WEG6 * i + 1) = 0.0;
      B(2, NODDOF_WEG6 * i + 2) = n_xyz(2, i);

      B(3, NODDOF_WEG6 * i + 0) = n_xyz(1, i);
      B(3, NODDOF_WEG6 * i + 1) = n_xyz(0, i);
      B(3, NODDOF_WEG6 * i + 2) = 0.0;

      B(4, NODDOF_WEG6 * i + 0) = 0.0;
      B(4, NODDOF_WEG6 * i + 1) = n_xyz(2, i);
      B(4, NODDOF_WEG6 * i + 2) = n_xyz(1, i);

      B(5, NODDOF_WEG6 * i + 0) = n_xyz(2, i);
      B(5, NODDOF_WEG6 * i + 1) = 0.0;
      B(5, NODDOF_WEG6 * i + 2) = n_xyz(0, i);
    }

    //--------------------------- build N_x operator (wrt spatial config)
    LINALG::Matrix<9, NUMDOF_WEG6> N_x(true);  // set to zero
    for (int i = 0; i < NUMNOD_WEG6; ++i)
    {
      N_x(0, 3 * i + 0) = n_xyz(0, i);
      N_x(1, 3 * i + 1) = n_xyz(0, i);
      N_x(2, 3 * i + 2) = n_xyz(0, i);

      N_x(3, 3 * i + 0) = n_xyz(1, i);
      N_x(4, 3 * i + 1) = n_xyz(1, i);
      N_x(5, 3 * i + 2) = n_xyz(1, i);

      N_x(6, 3 * i + 0) = n_xyz(2, i);
      N_x(7, 3 * i + 1) = n_xyz(2, i);
      N_x(8, 3 * i + 2) = n_xyz(2, i);
    }

    //------------------------------------------------- call material law
    LINALG::Matrix<6, 6> cmat(true);
    LINALG::Matrix<6, 1> stress(true);
    params.set<int>("gp", gp);
    Teuchos::RCP<MAT::So3Material> so3mat =
        Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
    so3mat->Evaluate(&F, &glstrain, params, &stress, &cmat, gp, ele->Id());

    //------------------------------------------- compute cauchy stresses
    LINALG::Matrix<6, 1> cstress;
    cstress.Multiply(j, IF, stress);

    //--------------------------------------- output strains and stresses
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == NULL) dserror("strain data not available");
        // rewriting Green-Lagrange strains in matrix format
        LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> gl;
        gl(0, 0) = glstrain(0);
        gl(0, 1) = 0.5 * glstrain(3);
        gl(0, 2) = 0.5 * glstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = glstrain(1);
        gl(1, 2) = 0.5 * glstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = glstrain(2);

        LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> temp;
        LINALG::Matrix<NUMDIM_WEG6, NUMDIM_WEG6> euler_almansi;
        temp.Multiply(gl, f);
        euler_almansi.MultiplyTN(f, temp);

        (*elestrain)(gp, 0) = euler_almansi(0, 0);
        (*elestrain)(gp, 1) = euler_almansi(1, 1);
        (*elestrain)(gp, 2) = euler_almansi(2, 2);
        (*elestrain)(gp, 3) = euler_almansi(0, 1);
        (*elestrain)(gp, 4) = euler_almansi(1, 2);
        (*elestrain)(gp, 5) = euler_almansi(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
    }

    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < 6; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < 6; ++i) (*elestress)(gp, i) = cstress(i);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
    }

    //-------------------------------------------- build operator Ypsilon
    LINALG::Matrix<3, 3> S;
    S(0, 0) = stress(0);
    S(0, 1) = stress(3);
    S(0, 2) = stress(5);
    S(1, 0) = stress(3);
    S(1, 1) = stress(1);
    S(1, 2) = stress(4);
    S(2, 0) = stress(5);
    S(2, 1) = stress(4);
    S(2, 2) = stress(2);
    LINALG::Matrix<6, 9> Ypsilon;
    BuildYpsilon(Ypsilon, F, S);

    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * (*weights)(gp);

    //------------------------------------------- assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).MultiplyTN(intfac, B, cstress, 1.0);
    }

    //*******************************************************************
    if (stiffmatrix)
    {
      LINALG::Matrix<6, 9> sum;

      //================================== 1.) B^T * cstress * F^T * N_x
      sum.MultiplyNT(cstress, FTvec);

      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      sum.Multiply(-j, Ypsilon, Lambda, 1.0);

      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
      LINALG::Matrix<6, 9> sixnine;
      LINALG::Matrix<6, 6> sixsix;
      sixnine.Multiply(Theta, Lambda);
      sixsix.Multiply(IF, cmat);
      sum.Multiply(-j, sixsix, sixnine, 1.0);

      //================ put everything together: K = intfac * B*T * sum * N_x
      LINALG::Matrix<6, NUMDOF_WEG6> tmp3;
      tmp3.Multiply(sum, N_x);
      (*stiffmatrix).MultiplyTN(intfac, B, tmp3, 1.0);

    }  // if (stiffmatrix)
    //*******************************************************************

    //*******************************************************************
    // Strictly, inverse design analysis is stationary and should not have
    // a mass term. Loosely, if no Dirichlet-BCs are present a small
    // mass term might be used. Note that we use density by unit deformed
    // volume here!
    if (massmatrix)
    {  // evaluate mass matrix +++++++++++++++++++++++++
      // integrate concistent mass matrix
      double density = ele->Material()->Density(gp);

      for (int inod = 0; inod < NUMNOD_WEG6; ++inod)
      {
        for (int jnod = 0; jnod < NUMNOD_WEG6; ++jnod)
        {
          double massfactor = (*shapefct)(inod, gp) * density * (*shapefct)(jnod, gp) * detj *
                              (*weights)(gp);  // intermediate factor
          (*massmatrix)(NUMDIM_WEG6 * inod + 0, NUMDIM_WEG6 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_WEG6 * inod + 1, NUMDIM_WEG6 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_WEG6 * inod + 2, NUMDIM_WEG6 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  }  // for (int gp=0; gp<NUMGPT_WEG6; ++gp)
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  return;
}  // DRT::ELEMENTS::InvDesign::sow6_nlnstiffmass


/*----------------------------------------------------------------------*
 |  compute and store material configuration                  (private) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::sow6_StoreMaterialConfiguration(
    DRT::ELEMENTS::So_weg6* ele, const std::vector<double>& disp)
{
  /* ============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for Wedge_6 with 6 GAUSS POINTS*
  ** ============================================================================*/
  /* pointer to (static) shape function array
   * for each node, evaluated at each gp*/
  LINALG::Matrix<NUMNOD_WEG6, NUMGPT_WEG6>* shapefct;
  /* pointer to (static) shape function derivatives array
   * for each node wrt to each direction, evaluated at each gp*/
  LINALG::Matrix<NUMGPT_WEG6 * NUMDIM_WEG6, NUMNOD_WEG6>* deriv;
  /* pointer to (static) weight factors at each gp */
  LINALG::Matrix<NUMGPT_WEG6, 1>* weights;
  ele->sow6_shapederiv(&shapefct, &deriv, &weights);  // call to evaluate
  /* ============================================================================*/

  LINALG::Matrix<NUMNOD_WEG6, 3> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_WEG6, 3> xcurr;  // current  coord. of element
  DRT::Node** nodes = ele->Nodes();
  for (int i = 0; i < NUMNOD_WEG6; ++i)
  {
    const double* x = nodes[i]->X();
    xcurr(i, 0) = x[0];
    xcurr(i, 1) = x[1];
    xcurr(i, 2) = x[2];

    xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_WEG6 + 0];
    xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_WEG6 + 1];
    xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_WEG6 + 2];
  }

  LINALG::Matrix<3, 3> invJ;
  LINALG::Matrix<3, NUMNOD_WEG6> N_XYZ;
  LINALG::Matrix<3, 3> F;
  for (int gp = 0; gp < NUMGPT_WEG6; ++gp)
  {
    // get submatrix of deriv at actual gp
    LINALG::Matrix<NUMDIM_WEG6, NUMGPT_WEG6> deriv_gp;
    for (int m = 0; m < NUMDIM_WEG6; ++m)
      for (int n = 0; n < NUMGPT_WEG6; ++n) deriv_gp(m, n) = (*deriv)(NUMDIM_WEG6 * gp + m, n);

    // compute invJ and detJ
    invJ.Multiply(deriv_gp, xrefe);
    detJ_[gp] = invJ.Invert();

    // compute F
    N_XYZ.Multiply(invJ, deriv_gp);
    F.MultiplyTT(xcurr, N_XYZ);

    // put stuff into storage
    MatrixtoStorage(gp, invJ, JHistory());
    MatrixtoStorage(gp, F, FHistory());
  }
  return;
}



/*----------------------------------------------------------------------*
 |  tet4 integration method                                    (public) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::so_tet4_nlnstiffmass(Teuchos::ParameterList& params,
    DRT::ELEMENTS::So_tet4* ele,                                ///< this element
    std::vector<int>& lm,                                       ///< location matrix
    std::vector<double>& disp,                                  ///< current displacements
    std::vector<double>& residual,                              ///< current residuum
    LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* stiffmatrix,  ///< element stiffness matrix
    LINALG::Matrix<NUMDOF_SOTET4, NUMDOF_SOTET4>* massmatrix,   ///< element mass matrix
    LINALG::Matrix<NUMDOF_SOTET4, 1>* force,                    ///< element internal force vector
    LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestress,  ///< stresses at GP
    LINALG::Matrix<NUMGPT_SOTET4, MAT::NUM_STRESS_3D>* elestrain,  ///< strains at GP
    const INPAR::STR::StressType iostress,                         ///< stress output options
    const INPAR::STR::StrainType iostrain)
{
  /* =============================================================================*
  ** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for TET_4  with 1 GAUSS POINTS*
  ** =============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
      ele->so_tet4_1gp_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
      ele->so_tet4_1gp_derivs();
  const static std::vector<double> gpweights = ele->so_tet4_1gp_weights();

  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xdisp;
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    xdisp(i, 0) = disp[i * NODDOF_SOTET4 + 0];
    xdisp(i, 1) = disp[i * NODDOF_SOTET4 + 1];
    xdisp(i, 2) = disp[i * NODDOF_SOTET4 + 2];
  }

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  const double detj = ele->V_;
  for (int gp = 0; gp < NUMGPT_SOTET4; gp++)
  {
    //----------------------------------- get inverse of Jacobian mapping
    //------------------------------ compute derivs wrt to spatial coords
    LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4>& n_xyz = ele->nxyz_;  // [gp];

    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> f;
    f.MultiplyTN(xdisp, n_xyz);
    f(0, 0) += 1;
    f(1, 1) += 1;
    f(2, 2) += 1;

    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> F(false);
    const double j = F.Invert(f);

    //------------------------------------ build F^T as vector 9x1
    LINALG::Matrix<9, 1> FTvec;
    FTvec(0) = F(0, 0);
    FTvec(1) = F(0, 1);
    FTvec(2) = F(0, 2);
    FTvec(3) = F(1, 0);
    FTvec(4) = F(1, 1);
    FTvec(5) = F(1, 2);
    FTvec(6) = F(2, 0);
    FTvec(7) = F(2, 1);
    FTvec(8) = F(2, 2);

    //--------------------------------------------- build operator Lambda
    LINALG::Matrix<9, 9> Lambda;
    BuildLambda(Lambda, F);

    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::Matrix<6, 6> IF;
    BuildIF(IF, F);

    //---------------------------------------------- build operator Theta
    LINALG::Matrix<6, 9> Theta;
    BuildTheta(Theta, F);

    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> cauchygreen;
    cauchygreen.MultiplyTN(F, F);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDOF_SOTET4> B;
    for (int i = 0; i < NUMNOD_SOTET4; ++i)
    {
      B(0, NODDOF_SOTET4 * i + 0) = n_xyz(i, 0);
      B(0, NODDOF_SOTET4 * i + 1) = 0.0;
      B(0, NODDOF_SOTET4 * i + 2) = 0.0;

      B(1, NODDOF_SOTET4 * i + 0) = 0.0;
      B(1, NODDOF_SOTET4 * i + 1) = n_xyz(i, 1);
      B(1, NODDOF_SOTET4 * i + 2) = 0.0;

      B(2, NODDOF_SOTET4 * i + 0) = 0.0;
      B(2, NODDOF_SOTET4 * i + 1) = 0.0;
      B(2, NODDOF_SOTET4 * i + 2) = n_xyz(i, 2);

      B(3, NODDOF_SOTET4 * i + 0) = n_xyz(i, 1);
      B(3, NODDOF_SOTET4 * i + 1) = n_xyz(i, 0);
      B(3, NODDOF_SOTET4 * i + 2) = 0.0;

      B(4, NODDOF_SOTET4 * i + 0) = 0.0;
      B(4, NODDOF_SOTET4 * i + 1) = n_xyz(i, 2);
      B(4, NODDOF_SOTET4 * i + 2) = n_xyz(i, 1);

      B(5, NODDOF_SOTET4 * i + 0) = n_xyz(i, 2);
      B(5, NODDOF_SOTET4 * i + 1) = 0.0;
      B(5, NODDOF_SOTET4 * i + 2) = n_xyz(i, 0);
    }

    //--------------------------- build N_x operator (wrt spatial config)
    LINALG::Matrix<9, NUMDOF_SOTET4> N_x(true);  // set to zero
    for (int i = 0; i < NUMNOD_SOTET4; ++i)
    {
      N_x(0, 3 * i + 0) = n_xyz(i, 0);
      N_x(1, 3 * i + 1) = n_xyz(i, 0);
      N_x(2, 3 * i + 2) = n_xyz(i, 0);

      N_x(3, 3 * i + 0) = n_xyz(i, 1);
      N_x(4, 3 * i + 1) = n_xyz(i, 1);
      N_x(5, 3 * i + 2) = n_xyz(i, 1);

      N_x(6, 3 * i + 0) = n_xyz(i, 2);
      N_x(7, 3 * i + 1) = n_xyz(i, 2);
      N_x(8, 3 * i + 2) = n_xyz(i, 2);
    }

    //------------------------------------------------- call material law
    LINALG::Matrix<6, 6> cmat(true);
    LINALG::Matrix<6, 1> stress(true);
    params.set<int>("gp", gp);
    Teuchos::RCP<MAT::So3Material> so3mat =
        Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
    so3mat->Evaluate(&F, &glstrain, params, &stress, &cmat, gp, ele->Id());

    //------------------------------------------- compute cauchy stresses
    LINALG::Matrix<6, 1> cstress;
    cstress.Multiply(j, IF, stress);

    //--------------------------------------- output strains and stresses
    switch (iostrain)
    {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("strain data not available");
        for (int i = 0; i < 3; ++i) (*elestrain)(gp, i) = glstrain(i);
        for (int i = 3; i < 6; ++i) (*elestrain)(gp, i) = 0.5 * glstrain(i);
      }
      break;
      case INPAR::STR::strain_ea:
      {
        if (elestrain == NULL) dserror("strain data not available");
        // rewriting Green-Lagrange strains in matrix format
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> gl;
        gl(0, 0) = glstrain(0);
        gl(0, 1) = 0.5 * glstrain(3);
        gl(0, 2) = 0.5 * glstrain(5);
        gl(1, 0) = gl(0, 1);
        gl(1, 1) = glstrain(1);
        gl(1, 2) = 0.5 * glstrain(4);
        gl(2, 0) = gl(0, 2);
        gl(2, 1) = gl(1, 2);
        gl(2, 2) = glstrain(2);

        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> temp;
        LINALG::Matrix<NUMDIM_SOTET4, NUMDIM_SOTET4> euler_almansi;
        temp.Multiply(gl, f);
        euler_almansi.MultiplyTN(f, temp);

        (*elestrain)(gp, 0) = euler_almansi(0, 0);
        (*elestrain)(gp, 1) = euler_almansi(1, 1);
        (*elestrain)(gp, 2) = euler_almansi(2, 2);
        (*elestrain)(gp, 3) = euler_almansi(0, 1);
        (*elestrain)(gp, 4) = euler_almansi(1, 2);
        (*elestrain)(gp, 5) = euler_almansi(0, 2);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      default:
        dserror("requested strain type not available");
    }

    switch (iostress)
    {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = stress(i);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("stress data not available");
        for (int i = 0; i < MAT::NUM_STRESS_3D; ++i) (*elestress)(gp, i) = cstress(i);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not available");
    }

    //-------------------------------------------- build operator Ypsilon
    LINALG::Matrix<3, 3> S;
    S(0, 0) = stress(0);
    S(0, 1) = stress(3);
    S(0, 2) = stress(5);
    S(1, 0) = stress(3);
    S(1, 1) = stress(1);
    S(1, 2) = stress(4);
    S(2, 0) = stress(5);
    S(2, 1) = stress(4);
    S(2, 2) = stress(2);
    LINALG::Matrix<6, 9> Ypsilon;
    BuildYpsilon(Ypsilon, F, S);

    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * (gpweights)[gp];

    //------------------------------------------- assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).MultiplyTN(intfac, B, cstress, 1.0);
    }

    //*******************************************************************
    if (stiffmatrix)
    {
      LINALG::Matrix<6, 9> sum;

      //================================== 1.) B^T * cstress * F^T * N_x
      sum.MultiplyNT(cstress, FTvec);

      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      sum.Multiply(-j, Ypsilon, Lambda, 1.0);

      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
      LINALG::Matrix<6, 9> sixnine;
      LINALG::Matrix<6, 6> sixsix;
      sixnine.Multiply(Theta, Lambda);
      sixsix.Multiply(IF, cmat);
      sum.Multiply(-j, sixsix, sixnine, 1.0);

      //================ put everything together: K = intfac * B*T * sum * N_x
      LINALG::Matrix<6, NUMDOF_SOTET4> tmp3;
      tmp3.Multiply(sum, N_x);
      (*stiffmatrix).MultiplyTN(intfac, B, tmp3, 1.0);
    }  // if (stiffmatrix)
    //*******************************************************************

  }  // for (int gp=0; gp<NUMGPT_SOTET4; gp++)
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  // static vectors created in any case to safe "if-case"
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts4gp =
      ele->so_tet4_4gp_shapefcts();
  const static std::vector<double> gpweights4gp = ele->so_tet4_4gp_weights();
  // evaluate mass matrix
  if (massmatrix != NULL)
  {
    double density =
        ele->Material()->Density(0);  // density at the only Gauss point the material has!

    // consistent mass matrix evaluated using a 4-point rule
    for (int gp = 0; gp < 4; gp++)
    {
      for (int inod = 0; inod < NUMNOD_SOTET4; ++inod)
      {
        for (int jnod = 0; jnod < NUMNOD_SOTET4; ++jnod)
        {
          double massfactor =
              (shapefcts4gp[gp])(inod)*density * (shapefcts4gp[gp])(jnod)*detj * (gpweights4gp)[gp];
          (*massmatrix)(NUMDIM_SOTET4 * inod + 0, NUMDIM_SOTET4 * jnod + 0) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4 * inod + 1, NUMDIM_SOTET4 * jnod + 1) += massfactor;
          (*massmatrix)(NUMDIM_SOTET4 * inod + 2, NUMDIM_SOTET4 * jnod + 2) += massfactor;
        }
      }
    }
  }  // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++


  return;
}  // DRT::ELEMENTS::InvDesign::so_tet4_nlnstiffmass


/*----------------------------------------------------------------------*
 |  compute and store material configuration                  (private) |
 |                                                             gee 09/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::sot4_StoreMaterialConfiguration(
    DRT::ELEMENTS::So_tet4* ele, const std::vector<double>& disp)
{
  const static std::vector<LINALG::Matrix<NUMNOD_SOTET4, 1>> shapefcts =
      ele->so_tet4_1gp_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOTET4 + 1, NUMNOD_SOTET4>> derivs =
      ele->so_tet4_1gp_derivs();
  const static std::vector<double> gpweights = ele->so_tet4_1gp_weights();
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xrefe;
  LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> xcurr;
  DRT::Node** nodes = ele->Nodes();
  for (int i = 0; i < NUMNOD_SOTET4; ++i)
  {
    const double* x = nodes[i]->X();
    xcurr(i, 0) = x[0];
    xcurr(i, 1) = x[1];
    xcurr(i, 2) = x[2];

    xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_SOTET4 + 0];
    xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_SOTET4 + 1];
    xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_SOTET4 + 2];
  }

  /* get the matrix of the coordinates of nodes needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  **         J = [ X_1  X_2  X_3  X_4 ]
  **             [ Y_1  Y_2  Y_3  Y_4 ]
  ** [ Z_1  Z_2  Z_3  Z_4 ]
  */
  LINALG::Matrix<NUMCOORD_SOTET4, NUMCOORD_SOTET4> J;
  for (int i = 0; i < 4; i++) J(0, i) = 1;
  for (int row = 0; row < 3; row++)
    for (int col = 0; col < 4; col++) J(row + 1, col) = xrefe(col, row);
  // volume of the element
  detJ_[0] = J.Determinant() / 6.0;

  for (int gp = 0; gp < NUMGPT_SOTET4; ++gp)
  {
    LINALG::Matrix<NUMCOORD_SOTET4, NUMCOORD_SOTET4> jac;
    LINALG::Matrix<NUMNOD_SOTET4, NUMDIM_SOTET4> N_XYZ;
    LINALG::Matrix<3, 3> F;

    {
      LINALG::Matrix<NUMCOORD_SOTET4 - 1, NUMCOORD_SOTET4> tmp;
      tmp.MultiplyTN(xrefe, derivs[gp]);
      for (int i = 0; i < 4; i++) jac(0, i) = 1;
      for (int row = 0; row < 3; row++)
        for (int col = 0; col < 4; col++) jac(row + 1, col) = tmp(row, col);
    }
    // size is 4x3
    LINALG::Matrix<NUMCOORD_SOTET4, NUMDIM_SOTET4> I_aug(true);  // set to zero
    // size is 4x3
    LINALG::Matrix<NUMCOORD_SOTET4, NUMDIM_SOTET4> partials(true);  // set to zero
    I_aug(1, 0) = 1;
    I_aug(2, 1) = 1;
    I_aug(3, 2) = 1;

    LINALG::FixedSizeSerialDenseSolver<NUMCOORD_SOTET4, NUMCOORD_SOTET4, NUMDIM_SOTET4>
        solve_for_inverseJac;                          // solve A.X=B
    solve_for_inverseJac.SetMatrix(jac);               // set A=jac
    solve_for_inverseJac.SetVectors(partials, I_aug);  // set X=partials, B=I_aug
    solve_for_inverseJac.FactorWithEquilibration(true);
    int err2 = solve_for_inverseJac.Factor();
    int err = solve_for_inverseJac.Solve();  // partials = jac^-1.I_aug
    if ((err != 0) && (err2 != 0)) dserror("Inversion of Jacobian failed");

    N_XYZ.Multiply(derivs[gp], partials);
    F.MultiplyTN(xcurr, N_XYZ);

    // put stuff into storage
    MatrixtoStorage(gp, N_XYZ, JHistory());
    MatrixtoStorage(gp, F, FHistory());

  }  // for (int gp=0; gp<NUMGPT_SOTET4; ++gp)

  return;
}
