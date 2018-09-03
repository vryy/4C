/*!----------------------------------------------------------------------
\file inversedesign_evaluate_fd.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "inversedesign.H"
#include "so_hex8.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  Lambda tensor (finite differenced)                        (private) |
 |  d f^{-1}_km / d f_pq                                                |
 |  which is the 9x9 nonsym Boeppel product                             |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FDLambda(
    LINALG::SerialDenseMatrix& L, const LINALG::SerialDenseMatrix& f) const

{
  // non-permuted original
  LINALG::SerialDenseMatrix F(f);
  LINALG::NonsymInverse3x3(F);

  double Lambda4[3][3][3][3];
  const double eps = 1.0e-8;
  for (int k = 0; k < 3; ++k)
    for (int m = 0; m < 3; ++m)
    {
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q)
        {
          LINALG::SerialDenseMatrix fperm(f);
          fperm(p, q) += eps;
          LINALG::NonsymInverse3x3(fperm);
          Lambda4[k][m][p][q] = (fperm(k, m) - F(k, m)) / eps;
        }
    }



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

  return;
}  // DRT::ELEMENTS::InvDesign::FDLambda


/*----------------------------------------------------------------------*
 |  Lambda tensor (finite differenced)                        (private) |
 |  d f^{-1}_km / d f_pq                                                |
 |  which is the 9x9 nonsym Boeppel product                             |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FDLambdaT(
    LINALG::SerialDenseMatrix& L, const LINALG::SerialDenseMatrix& f) const

{
  // non-permuted original function
  LINALG::SerialDenseMatrix F(f);
  LINALG::NonsymInverse3x3(F);
  LINALG::SerialDenseMatrix FT(3, 3);
  FT(0, 0) = F(0, 0);
  FT(0, 1) = F(1, 0);
  FT(0, 2) = F(2, 0);
  FT(1, 0) = F(0, 1);
  FT(1, 1) = F(1, 1);
  FT(1, 2) = F(2, 1);
  FT(2, 0) = F(0, 2);
  FT(2, 1) = F(1, 2);
  FT(2, 2) = F(2, 2);

  double Lambda4[3][3][3][3];
  const double eps = 1.0e-8;
  for (int k = 0; k < 3; ++k)
    for (int m = 0; m < 3; ++m)
    {
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q)
        {
          LINALG::SerialDenseMatrix fperm(f);
          fperm(p, q) += eps;
          // evaluate the function, which is inverse followed by transpose
          LINALG::NonsymInverse3x3(fperm);
          LINALG::SerialDenseMatrix fpermT(3, 3);
          fpermT(0, 0) = fperm(0, 0);
          fpermT(0, 1) = fperm(1, 0);
          fpermT(0, 2) = fperm(2, 0);
          fpermT(1, 0) = fperm(0, 1);
          fpermT(1, 1) = fperm(1, 1);
          fpermT(1, 2) = fperm(2, 1);
          fpermT(2, 0) = fperm(0, 2);
          fpermT(2, 1) = fperm(1, 2);
          fpermT(2, 2) = fperm(2, 2);
          Lambda4[k][m][p][q] = (fpermT(k, m) - FT(k, m)) / eps;
        }
    }



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

  return;
}  // DRT::ELEMENTS::InvDesign::FDLambdaT


/*----------------------------------------------------------------------*
 |  Theta tensor                                              (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FDTheta(double Theta4[][3][3][3], LINALG::SerialDenseMatrix& Theta,
    const LINALG::SerialDenseMatrix& F) const

{
  // non-permuted original function
  LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8, NUMDIM_SOH8);
  cauchygreen.Multiply('T', 'N', 1.0, F, F, 0.0);
  cauchygreen(0, 0) -= 1.0;
  cauchygreen(1, 1) -= 1.0;
  cauchygreen(2, 1) -= 1.0;
  LINALG::SerialDenseMatrix E(cauchygreen);
  E.Scale(0.5);

  const double eps = 1.0e-8;

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l)
        {
          LINALG::SerialDenseMatrix Fperm(F);
          Fperm(k, l) += eps;
          // evaluate Eperm
          LINALG::SerialDenseMatrix cauchygreenperm(NUMDIM_SOH8, NUMDIM_SOH8);
          cauchygreenperm.Multiply('T', 'N', 1.0, Fperm, Fperm, 0.0);
          cauchygreenperm(0, 0) -= 1.0;
          cauchygreenperm(1, 1) -= 1.0;
          cauchygreenperm(2, 1) -= 1.0;
          LINALG::SerialDenseMatrix Eperm(cauchygreenperm);
          Eperm.Scale(0.5);
          // do the differencing
          Theta4[i][j][k][l] = (Eperm(i, j) - E(i, j)) / eps;
        }
    }


  Theta(0, 0) = Theta4[0][0][0][0];  // ok
  Theta(0, 1) = Theta4[0][0][1][0];
  Theta(0, 2) = Theta4[0][0][2][0];  // nok
  Theta(0, 3) = Theta4[0][0][0][1];
  Theta(0, 4) = Theta4[0][0][1][1];
  Theta(0, 5) = Theta4[0][0][2][1];
  Theta(0, 6) = Theta4[0][0][0][2];  // nok
  Theta(0, 7) = Theta4[0][0][1][2];
  Theta(0, 8) = Theta4[0][0][2][2];

  Theta(1, 0) = Theta4[1][1][0][0];
  Theta(1, 1) = Theta4[1][1][1][0];
  Theta(1, 2) = Theta4[1][1][2][0];
  Theta(1, 3) = Theta4[1][1][0][1];
  Theta(1, 4) = Theta4[1][1][1][1];  // ok
  Theta(1, 5) = Theta4[1][1][2][1];
  Theta(1, 6) = Theta4[1][1][0][2];
  Theta(1, 7) = Theta4[1][1][1][2];
  Theta(1, 8) = Theta4[1][1][2][2];

  Theta(2, 0) = Theta4[2][2][0][0];
  Theta(2, 1) = Theta4[2][2][1][0];
  Theta(2, 2) = Theta4[2][2][2][0];  // nok
  Theta(2, 3) = Theta4[2][2][0][1];
  Theta(2, 4) = Theta4[2][2][1][1];
  Theta(2, 5) = Theta4[2][2][2][1];
  Theta(2, 6) = Theta4[2][2][0][2];  // nok
  Theta(2, 7) = Theta4[2][2][1][2];
  Theta(2, 8) = Theta4[2][2][2][2];  // ok

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

  return;
}  // DRT::ELEMENTS::InvDesign::FDTheta



/*----------------------------------------------------------------------*
 |  Ypsilon tensor                                            (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FDYpsilon(double Y4[][3][3][3], LINALG::SerialDenseMatrix& Y,
    const LINALG::SerialDenseMatrix& F, const LINALG::SerialDenseMatrix& S) const

{
  // non-permuted original function
  Epetra_SerialDenseMatrix IFS(3, 3);
  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
    {
      IFS(k, l) = 0.0;
      for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
          IFS(k, l) += 0.5 * (F(k, m) * F(l, n) + F(k, n) * F(l, m)) * S(m, n);
    }

  const double eps = 1.0e-8;

  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
    {
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q)
        {
          // permute F -> Fperm
          LINALG::SerialDenseMatrix Fperm(F);
          Fperm(p, q) += eps;
          // evaluate IFSperm
          Epetra_SerialDenseMatrix IFSperm(3, 3);
          for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
            {
              IFSperm(a, b) = 0.0;
              for (int c = 0; c < 3; ++c)
                for (int d = 0; d < 3; ++d)
                  IFSperm(a, b) +=
                      0.5 * (Fperm(a, c) * Fperm(b, d) + Fperm(a, d) * Fperm(b, c)) * S(c, d);
            }
          // do the differencing
          Y4[k][l][p][q] = (IFSperm(k, l) - IFS(k, l)) / eps;
        }  // pq
    }      // kl



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

  return;
}  // DRT::ELEMENTS::InvDesign::FDYpsilon


/*----------------------------------------------------------------------*
 |  FD whole stiffness                                        (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FDstiffmatrix(Epetra_SerialDenseMatrix& stiff,
    const std::vector<double>& disp, const int gp, DRT::ELEMENTS::So_hex8* ele,
    Teuchos::ParameterList& params) const

{
  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;

  //********************** non-permuted original function cstress(disp)
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8, NUMDIM_SOH8);
  Epetra_SerialDenseVector cstress(MAT::NUM_STRESS_3D);
  {
    LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8, NUMDIM_SOH8);
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      xcurr(i, 0) = ele->Nodes()[i]->X()[0];
      xcurr(i, 1) = ele->Nodes()[i]->X()[1];
      xcurr(i, 2) = ele->Nodes()[i]->X()[2];

      xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_SOH8 + 0];
      xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_SOH8 + 1];
      xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_SOH8 + 2];
    }

    // no gauss point loop, we just look at THE gauss point gp here
    const double detj = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];

    LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8, NUMNOD_SOH8);
    n_xyz.Multiply('N', 'N', 1.0, invj, int_hex8.deriv_gp[gp], 0.0);

    LINALG::SerialDenseMatrix f(NUMDIM_SOH8, NUMDIM_SOH8);
    f.Multiply('T', 'T', 1.0, xrefe, n_xyz, 0.0);

    LINALG::SerialDenseMatrix F(f);
    const double detf = LINALG::NonsymInverse3x3(F);
    const double detF = 1.0 / detf;

    LINALG::SerialDenseMatrix IF(6, 6);
    BuildIF(IF, F);

    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8, NUMDIM_SOH8);
    cauchygreen.Multiply('T', 'N', 1.0, F, F, 0.0);

    LINALG::SerialDenseVector glstrain(MAT::NUM_STRESS_3D);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    Epetra_SerialDenseMatrix cmat(MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D);
    Epetra_SerialDenseVector stress(MAT::NUM_STRESS_3D);
    params.set<int>("gp", gp);
    Teuchos::RCP<MAT::So3Material> so3mat =
        Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
    so3mat->Evaluate(&F, &glstrain, params, &stress, &cmat, ele->Id());

    cstress.Multiply('N', 'N', detf, IF, stress, 0.0);
    // std::cout << "unperm cstress\n" << cstress;
  }

  const double eps = 1.0e-8;

  for (int i = 0; i < 6; ++i)
  {
    std::vector<double> dispperm(NUMDOF_SOH8);
    for (int j = 0; j < NUMDOF_SOH8; ++j)
    {
      // permuted the displacement in direction j
      dispperm = disp;
      dispperm[j] += eps;

      // evaluated cstressperm
      LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8, NUMDIM_SOH8);
      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        xrefe(k, 0) = xcurr(k, 0) + dispperm[k * NODDOF_SOH8 + 0];
        xrefe(k, 1) = xcurr(k, 1) + dispperm[k * NODDOF_SOH8 + 1];
        xrefe(k, 2) = xcurr(k, 2) + dispperm[k * NODDOF_SOH8 + 2];
      }
      const double detj = ele->detJ_[gp];
      Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];
      LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8, NUMNOD_SOH8);
      n_xyz.Multiply('N', 'N', 1.0, invj, int_hex8.deriv_gp[gp], 0.0);

      LINALG::SerialDenseMatrix f(NUMDIM_SOH8, NUMDIM_SOH8);
      f.Multiply('T', 'T', 1.0, xrefe, n_xyz, 0.0);

      LINALG::SerialDenseMatrix F(f);
      const double detf = LINALG::NonsymInverse3x3(F);
      const double detF = 1.0 / detf;

      LINALG::SerialDenseMatrix IF(6, 6);
      BuildIF(IF, F);

      LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8, NUMDIM_SOH8);
      cauchygreen.Multiply('T', 'N', 1.0, F, F, 0.0);

      LINALG::SerialDenseVector glstrain(MAT::NUM_STRESS_3D);
      glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
      glstrain(3) = cauchygreen(0, 1);
      glstrain(4) = cauchygreen(1, 2);
      glstrain(5) = cauchygreen(2, 0);

      Epetra_SerialDenseMatrix cmat(MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D);
      Epetra_SerialDenseVector stress(MAT::NUM_STRESS_3D);
      params.set<int>("gp", gp);
      Teuchos::RCP<MAT::So3Material> so3mat =
          Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
      so3mat->Evaluate(&F, &glstrain, params, &stress, &cmat, ele->Id());

      Epetra_SerialDenseVector cstressperm(MAT::NUM_STRESS_3D);
      cstressperm.Multiply('N', 'N', detf, IF, stress, 0.0);
      // std::cout << "cstressperm (i,j)=(" << i << "," << j << ")\n" << cstressperm;


      // do the finite differencing
      stiff(i, j) = (cstressperm(i) - cstress(i)) / eps;

    }  // j
  }    // i



  return;
}  // DRT::ELEMENTS::InvDesign::FDstiffmatrix



/*----------------------------------------------------------------------*
 |  FD dj/dX                                                  (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FD_djdX(Epetra_SerialDenseMatrix& djdX,
    const std::vector<double>& disp, const int gp, DRT::ELEMENTS::So_hex8* ele) const

{
  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;

  //********************** non-permuted original function cstress(disp)
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8, NUMDIM_SOH8);
  double detf = 0.0;
  {
    LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8, NUMDIM_SOH8);
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      xcurr(i, 0) = ele->Nodes()[i]->X()[0];
      xcurr(i, 1) = ele->Nodes()[i]->X()[1];
      xcurr(i, 2) = ele->Nodes()[i]->X()[2];

      xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_SOH8 + 0];
      xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_SOH8 + 1];
      xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_SOH8 + 2];
    }

    // no gauss point loop, we just look at THE gauss point gp here
    const double detj = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];

    LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8, NUMNOD_SOH8);
    n_xyz.Multiply('N', 'N', 1.0, invj, int_hex8.deriv_gp[gp], 0.0);

    LINALG::SerialDenseMatrix f(NUMDIM_SOH8, NUMDIM_SOH8);
    f.Multiply('T', 'T', 1.0, xrefe, n_xyz, 0.0);

    LINALG::SerialDenseMatrix F(f);
    detf = LINALG::NonsymInverse3x3(F);
  }

  const double eps = 1.0e-8;

  std::vector<double> dispperm(NUMDOF_SOH8);
  for (int j = 0; j < NUMDOF_SOH8; ++j)
  {
    // permuted the displacement in direction j
    dispperm = disp;
    dispperm[j] += eps;

    // evaluated detjperm
    LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8, NUMDIM_SOH8);
    for (int k = 0; k < NUMNOD_SOH8; ++k)
    {
      xrefe(k, 0) = xcurr(k, 0) + dispperm[k * NODDOF_SOH8 + 0];
      xrefe(k, 1) = xcurr(k, 1) + dispperm[k * NODDOF_SOH8 + 1];
      xrefe(k, 2) = xcurr(k, 2) + dispperm[k * NODDOF_SOH8 + 2];
    }
    const double detj = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];
    LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8, NUMNOD_SOH8);
    n_xyz.Multiply('N', 'N', 1.0, invj, int_hex8.deriv_gp[gp], 0.0);

    LINALG::SerialDenseMatrix f(NUMDIM_SOH8, NUMDIM_SOH8);
    f.Multiply('T', 'T', 1.0, xrefe, n_xyz, 0.0);

    LINALG::SerialDenseMatrix F(f);
    const double detfperm = LINALG::NonsymInverse3x3(F);

    // do the finite differencing
    djdX(0, j) = (detfperm - detf) / eps;

  }  // j



  return;
}  // DRT::ELEMENTS::InvDesign::FD_djdX



/*----------------------------------------------------------------------*
 |  FD d IS / d X                                             (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FD_dISdX(Epetra_SerialDenseMatrix& stiff,
    const std::vector<double>& disp, const int gp, DRT::ELEMENTS::So_hex8* ele,
    Teuchos::ParameterList& params) const

{
  const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;

  //********************** non-permuted original function cstress(disp)
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8, NUMDIM_SOH8);
  Epetra_SerialDenseVector cstress(MAT::NUM_STRESS_3D);
  Epetra_SerialDenseVector stress(MAT::NUM_STRESS_3D);
  double detf;
  {
    LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8, NUMDIM_SOH8);
    for (int i = 0; i < NUMNOD_SOH8; ++i)
    {
      xcurr(i, 0) = ele->Nodes()[i]->X()[0];
      xcurr(i, 1) = ele->Nodes()[i]->X()[1];
      xcurr(i, 2) = ele->Nodes()[i]->X()[2];

      xrefe(i, 0) = xcurr(i, 0) + disp[i * NODDOF_SOH8 + 0];
      xrefe(i, 1) = xcurr(i, 1) + disp[i * NODDOF_SOH8 + 1];
      xrefe(i, 2) = xcurr(i, 2) + disp[i * NODDOF_SOH8 + 2];
    }

    // no gauss point loop, we just look at THE gauss point gp here
    const double detj = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];

    LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8, NUMNOD_SOH8);
    n_xyz.Multiply('N', 'N', 1.0, invj, int_hex8.deriv_gp[gp], 0.0);

    LINALG::SerialDenseMatrix f(NUMDIM_SOH8, NUMDIM_SOH8);
    f.Multiply('T', 'T', 1.0, xrefe, n_xyz, 0.0);

    LINALG::SerialDenseMatrix F(f);
    detf = LINALG::NonsymInverse3x3(F);
    // const double detF = 1.0/detf;

    LINALG::SerialDenseMatrix IF(6, 6);
    BuildIF(IF, F);

    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8, NUMDIM_SOH8);
    cauchygreen.Multiply('T', 'N', 1.0, F, F, 0.0);

    LINALG::SerialDenseVector glstrain(MAT::NUM_STRESS_3D);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    Epetra_SerialDenseMatrix cmat(MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D);
    params.set<int>("gp", gp);
    Teuchos::RCP<MAT::So3Material> so3mat =
        Teuchos::rcp_dynamic_cast<MAT::So3Material>(ele->Material());
    so3mat->Evaluate(&F, &glstrain, params, &stress, &cmat, ele->Id());

    cstress.Multiply('N', 'N', detf, IF, stress, 0.0);
    // std::cout << "unperm cstress\n" << cstress;
  }

  const double eps = 1.0e-8;

  for (int i = 0; i < 6; ++i)
  {
    std::vector<double> dispperm(NUMDOF_SOH8);
    for (int j = 0; j < NUMDOF_SOH8; ++j)
    {
      // permuted the displacement in direction j
      dispperm = disp;
      dispperm[j] += eps;

      // evaluated cstressperm
      LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8, NUMDIM_SOH8);
      for (int k = 0; k < NUMNOD_SOH8; ++k)
      {
        xrefe(k, 0) = xcurr(k, 0) + dispperm[k * NODDOF_SOH8 + 0];
        xrefe(k, 1) = xcurr(k, 1) + dispperm[k * NODDOF_SOH8 + 1];
        xrefe(k, 2) = xcurr(k, 2) + dispperm[k * NODDOF_SOH8 + 2];
      }
      const double detj = ele->detJ_[gp];
      Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];
      LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8, NUMNOD_SOH8);
      n_xyz.Multiply('N', 'N', 1.0, invj, int_hex8.deriv_gp[gp], 0.0);

      LINALG::SerialDenseMatrix f(NUMDIM_SOH8, NUMDIM_SOH8);
      f.Multiply('T', 'T', 1.0, xrefe, n_xyz, 0.0);

      LINALG::SerialDenseMatrix F(f);
      LINALG::NonsymInverse3x3(F);
      // const double detF = 1.0/detf;

      LINALG::SerialDenseMatrix IF(6, 6);
      BuildIF(IF, F);

      Epetra_SerialDenseVector cstressperm(MAT::NUM_STRESS_3D);
      cstressperm.Multiply('N', 'N', detf, IF, stress, 0.0);
      // std::cout << "cstressperm (i,j)=(" << i << "," << j << ")\n" << cstressperm;


      // do the finite differencing
      stiff(i, j) = (cstressperm(i) - cstress(i)) / eps;

    }  // j
  }    // i



  return;
}  // DRT::ELEMENTS::InvDesign::FD_dISdX



/*----------------------------------------------------------------------*
 |  FD d IS / d X                                             (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::FD_dISdf(Epetra_SerialDenseMatrix& dISdf,
    const Epetra_SerialDenseVector& stress, const LINALG::SerialDenseMatrix& f) const

{
  //********************** non-permuted original function cstress(disp)
  Epetra_SerialDenseVector cstress(MAT::NUM_STRESS_3D);
  double detf;
  {
    LINALG::SerialDenseMatrix F(f);
    detf = LINALG::NonsymInverse3x3(F);

    LINALG::SerialDenseMatrix IF(6, 6);
    BuildIF(IF, F);

    // cstress.Multiply('N','N',detf,IF,stress,0.0);
    cstress.Multiply('N', 'N', 1.0, IF, stress, 0.0);  // detf is done outside
  }

  const double eps = 1.0e-8;

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 9; ++j)
    {
      // permute f
      Epetra_SerialDenseVector fpermvec(9);
      fpermvec(0) = f(0, 0);
      fpermvec(1) = f(1, 0);
      fpermvec(2) = f(2, 0);
      fpermvec(3) = f(0, 1);
      fpermvec(4) = f(1, 1);
      fpermvec(5) = f(2, 1);
      fpermvec(6) = f(0, 2);
      fpermvec(7) = f(1, 2);
      fpermvec(8) = f(2, 2);

      fpermvec(j) += eps;

      LINALG::SerialDenseMatrix fperm(3, 3);
      fperm(0, 0) = fpermvec(0);
      fperm(1, 0) = fpermvec(1);
      fperm(2, 0) = fpermvec(2);
      fperm(0, 1) = fpermvec(3);
      fperm(1, 1) = fpermvec(4);
      fperm(2, 1) = fpermvec(5);
      fperm(0, 2) = fpermvec(6);
      fperm(1, 2) = fpermvec(7);
      fperm(2, 2) = fpermvec(8);

      LINALG::SerialDenseMatrix Fperm(fperm);
      LINALG::NonsymInverse3x3(Fperm);

      LINALG::SerialDenseMatrix IFperm(6, 6);
      BuildIF(IFperm, Fperm);

      Epetra_SerialDenseVector cstressperm(MAT::NUM_STRESS_3D);
      // cstressperm.Multiply('N','N',detf,IFperm,stress,0.0);
      cstressperm.Multiply('N', 'N', 1.0, IFperm, stress, 0.0);  // detf is done outside


      // do the finite differencing
      dISdf(i, j) = (cstressperm(i) - cstress(i)) / eps;
    }  // j


  return;
}  // DRT::ELEMENTS::InvDesign::FD_dISdf


/*----------------------------------------------------------------------*
 |  Ypsilon tensor                                            (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::TensorMultiply(
    Epetra_SerialDenseMatrix& sum, double fac, double Y4[][3][3][3], double L4[][3][3][3])

{
  double result[3][3][3][3];

  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l)
    {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
          result[k][l][i][j] = 0.0;
          for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q) result[k][l][i][j] += Y4[k][l][p][q] * L4[p][q][i][j];
        }
    }


  sum(0, 0) = result[0][0][0][0];
  sum(0, 1) = result[0][0][1][0];
  sum(0, 2) = result[0][0][2][0];
  sum(0, 3) = result[0][0][0][1];
  sum(0, 4) = result[0][0][1][1];
  sum(0, 5) = result[0][0][2][1];
  sum(0, 6) = result[0][0][0][2];
  sum(0, 7) = result[0][0][1][2];
  sum(0, 8) = result[0][0][2][2];

  sum(1, 0) = result[1][1][0][0];
  sum(1, 1) = result[1][1][1][0];
  sum(1, 2) = result[1][1][2][0];
  sum(1, 3) = result[1][1][0][1];
  sum(1, 4) = result[1][1][1][1];
  sum(1, 5) = result[1][1][2][1];
  sum(1, 6) = result[1][1][0][2];
  sum(1, 7) = result[1][1][1][2];
  sum(1, 8) = result[1][1][2][2];

  sum(2, 0) = result[2][2][0][0];
  sum(2, 1) = result[2][2][1][0];
  sum(2, 2) = result[2][2][2][0];
  sum(2, 3) = result[2][2][0][1];
  sum(2, 4) = result[2][2][1][1];
  sum(2, 5) = result[2][2][2][1];
  sum(2, 6) = result[2][2][0][2];
  sum(2, 7) = result[2][2][1][2];
  sum(2, 8) = result[2][2][2][2];

  sum(3, 0) = result[0][1][0][0];
  sum(3, 1) = result[0][1][1][0];
  sum(3, 2) = result[0][1][2][0];
  sum(3, 3) = result[0][1][0][1];
  sum(3, 4) = result[0][1][1][1];
  sum(3, 5) = result[0][1][2][1];
  sum(3, 6) = result[0][1][0][2];
  sum(3, 7) = result[0][1][1][2];
  sum(3, 8) = result[0][1][2][2];

  sum(4, 0) = result[1][2][0][0];
  sum(4, 1) = result[1][2][1][0];
  sum(4, 2) = result[1][2][2][0];
  sum(4, 3) = result[1][2][0][1];
  sum(4, 4) = result[1][2][1][1];
  sum(4, 5) = result[1][2][2][1];
  sum(4, 6) = result[1][2][0][2];
  sum(4, 7) = result[1][2][1][2];
  sum(4, 8) = result[1][2][2][2];

  sum(5, 0) = result[2][0][0][0];
  sum(5, 1) = result[2][0][1][0];
  sum(5, 2) = result[2][0][2][0];
  sum(5, 3) = result[2][0][0][1];
  sum(5, 4) = result[2][0][1][1];
  sum(5, 5) = result[2][0][2][1];
  sum(5, 6) = result[2][0][0][2];
  sum(5, 7) = result[2][0][1][2];
  sum(5, 8) = result[2][0][2][2];

  sum.Scale(fac);

  return;
}  // DRT::ELEMENTS::InvDesign::BuildYpsilon
