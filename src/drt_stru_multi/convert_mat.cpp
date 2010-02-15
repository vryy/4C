/*!----------------------------------------------------------------------
\file convert_mat.cpp

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "microstatic.H"

void STRUMULTI::MicroStatic::ConvertMat(const Epetra_MultiVector& cmatpf,
                                        const LINALG::Matrix<3,3>& F_inv,
                                        const LINALG::Matrix<6,1>& S,
                                        LINALG::Matrix<6,6>& cmat)
{
  // For details concerning the conversion of the constitutive tensor
  // relating first Piola Kirchhoff stresses and deformation gradient
  // into a constitutive tensor relating second Piola Kirchhoff
  // stresses and Green-Lagrange strains cf.
  //
  // Marsden and Hughes, Mathematical Foundations of Elasticity,
  // Dover, pg. 215.

  int StressIndex[6][6] =
    {
      {0,3,5,3,5,0},
      {3,1,4,1,4,3},
      {5,4,2,4,2,5},
      {3,1,4,1,4,3},
      {5,4,2,4,2,5},
      {0,3,5,3,5,0}
    };

  int cmatpfIndex[6][6] =
    {
      {0,7,5},
      {3,1,8},
      {6,4,2},
      {3,1,8},
      {6,4,2},
      {0,7,5}
    };

  int FinvIndex[6] = {0,1,2,0,1,2};


  for (int i=0; i<6; ++i)
  {
    for (int j=0; j<6; ++j)
    {
      cmat(i,j) = ((*cmatpf(cmatpfIndex[j][0]))[cmatpfIndex[i][0]]-S(StressIndex[i][j]))*F_inv(FinvIndex[i],0)*F_inv(FinvIndex[j],0)
                  +(*cmatpf(cmatpfIndex[j][1]))[cmatpfIndex[i][0]]                      *F_inv(FinvIndex[i],0)*F_inv(FinvIndex[j],1)
                  +(*cmatpf(cmatpfIndex[j][2]))[cmatpfIndex[i][0]]                      *F_inv(FinvIndex[i],0)*F_inv(FinvIndex[j],2)
                  +(*cmatpf(cmatpfIndex[j][0]))[cmatpfIndex[i][1]]                      *F_inv(FinvIndex[i],1)*F_inv(FinvIndex[j],0)
                 +((*cmatpf(cmatpfIndex[j][1]))[cmatpfIndex[i][1]]-S(StressIndex[i][j]))*F_inv(FinvIndex[i],1)*F_inv(FinvIndex[j],1)
                  +(*cmatpf(cmatpfIndex[j][2]))[cmatpfIndex[i][1]]                      *F_inv(FinvIndex[i],1)*F_inv(FinvIndex[j],2)
                  +(*cmatpf(cmatpfIndex[j][0]))[cmatpfIndex[i][2]]                      *F_inv(FinvIndex[i],2)*F_inv(FinvIndex[j],0)
                  +(*cmatpf(cmatpfIndex[j][1]))[cmatpfIndex[i][2]]                      *F_inv(FinvIndex[i],2)*F_inv(FinvIndex[j],1)
                 +((*cmatpf(cmatpfIndex[j][2]))[cmatpfIndex[i][2]]-S(StressIndex[i][j]))*F_inv(FinvIndex[i],2)*F_inv(FinvIndex[j],2);
    }
  }

  return;
}

#endif
