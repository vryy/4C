/*---------------------------------------------------------------------*/
/*! \file

\brief Material conversion

\maintainer Martin Kronbichler

\level 2

*/
/*---------------------------------------------------------------------*/


#include "microstatic.H"

void STRUMULTI::MicroStatic::ConvertMat(const Epetra_MultiVector& cmatpf,
    const LINALG::Matrix<3, 3>& F_inv, const LINALG::Matrix<6, 1>& S, LINALG::Matrix<6, 6>& cmat)
{
  // For details concerning the conversion of the constitutive tensor
  // relating first Piola Kirchhoff stresses and deformation gradient
  // into a constitutive tensor relating second Piola Kirchhoff
  // stresses and Green-Lagrange strains cf.
  //
  // Marsden and Hughes, Mathematical Foundations of Elasticity,
  // Dover, pg. 215.

  int StressIndex[24][6] = {{0, 3, 5, 3, 5, 0}, {3, 1, 4, 1, 4, 3}, {5, 4, 2, 4, 2, 5},
      {3, 1, 4, 1, 4, 3}, {5, 4, 2, 4, 2, 5}, {0, 3, 5, 3, 5, 0},

      {0, 3, 5, 0, 3, 5}, {3, 1, 4, 3, 1, 4}, {5, 4, 2, 5, 4, 2}, {3, 1, 4, 3, 1, 4},
      {5, 4, 2, 5, 4, 2}, {0, 3, 5, 0, 3, 5},

      {0, 3, 5, 3, 5, 0}, {3, 1, 4, 1, 4, 3}, {5, 4, 2, 4, 2, 5}, {0, 3, 5, 3, 5, 0},
      {3, 1, 4, 1, 4, 3}, {5, 4, 2, 4, 2, 5},

      {0, 3, 5, 0, 3, 5}, {3, 1, 4, 3, 1, 4}, {5, 4, 2, 5, 4, 2}, {0, 3, 5, 0, 3, 5},
      {3, 1, 4, 3, 1, 4}, {5, 4, 2, 5, 4, 2}};


  int cmatpfIndex1[24][3] = {
      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},

      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},

      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},

      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
  };

  int cmatpfIndex2[24][3] = {
      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},

      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},

      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},

      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
      {0, 7, 5},
      {3, 1, 8},
      {6, 4, 2},
  };

  int FinvIndex1[4][6] = {
      {0, 1, 2, 0, 1, 2}, {0, 1, 2, 0, 1, 2}, {0, 1, 2, 1, 2, 0}, {0, 1, 2, 1, 2, 0}};

  int FinvIndex2[4][6] = {
      {0, 1, 2, 0, 1, 2}, {0, 1, 2, 1, 2, 0}, {0, 1, 2, 0, 1, 2}, {0, 1, 2, 1, 2, 0}};

  cmat.PutScalar(0.0);

  for (int i = 0; i < 6; ++i)
  {
    for (int j = 0; j < 6; ++j)
    {
      for (int k = 0; k < 4; ++k)
      {
        cmat(i, j) += ((*cmatpf(cmatpfIndex2[j + k * 6][0]))[cmatpfIndex1[i + k * 6][0]] -
                          S(StressIndex[i + k * 6][j])) *
                          F_inv(FinvIndex1[k][i], 0) * F_inv(FinvIndex2[k][j], 0) +
                      (*cmatpf(cmatpfIndex2[j + k * 6][1]))[cmatpfIndex1[i + k * 6][0]] *
                          F_inv(FinvIndex1[k][i], 0) * F_inv(FinvIndex2[k][j], 1) +
                      (*cmatpf(cmatpfIndex2[j + k * 6][2]))[cmatpfIndex1[i + k * 6][0]] *
                          F_inv(FinvIndex1[k][i], 0) * F_inv(FinvIndex2[k][j], 2) +
                      (*cmatpf(cmatpfIndex2[j + k * 6][0]))[cmatpfIndex1[i + k * 6][1]] *
                          F_inv(FinvIndex1[k][i], 1) * F_inv(FinvIndex2[k][j], 0) +
                      ((*cmatpf(cmatpfIndex2[j + k * 6][1]))[cmatpfIndex1[i + k * 6][1]] -
                          S(StressIndex[i + k * 6][j])) *
                          F_inv(FinvIndex1[k][i], 1) * F_inv(FinvIndex2[k][j], 1) +
                      (*cmatpf(cmatpfIndex2[j + k * 6][2]))[cmatpfIndex1[i + k * 6][1]] *
                          F_inv(FinvIndex1[k][i], 1) * F_inv(FinvIndex2[k][j], 2) +
                      (*cmatpf(cmatpfIndex2[j + k * 6][0]))[cmatpfIndex1[i + k * 6][2]] *
                          F_inv(FinvIndex1[k][i], 2) * F_inv(FinvIndex2[k][j], 0) +
                      (*cmatpf(cmatpfIndex2[j + k * 6][1]))[cmatpfIndex1[i + k * 6][2]] *
                          F_inv(FinvIndex1[k][i], 2) * F_inv(FinvIndex2[k][j], 1) +
                      ((*cmatpf(cmatpfIndex2[j + k * 6][2]))[cmatpfIndex1[i + k * 6][2]] -
                          S(StressIndex[i + k * 6][j])) *
                          F_inv(FinvIndex1[k][i], 2) * F_inv(FinvIndex2[k][j], 2);
      }
    }
  }

  cmat.Scale(0.25);

  return;
}
