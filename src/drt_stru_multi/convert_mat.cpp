#ifdef CCADISCRET

#include "microstatic.H"

void MicroStatic::ConvertMat(const Epetra_MultiVector& cmatpf,
                             const Epetra_SerialDenseMatrix& F_inv,
                             const Epetra_SerialDenseVector& S,
                             Epetra_SerialDenseMatrix& cmat)
{
  // For details concerning the conversion of the constitutive tensor
  // relating first Piola Kirchhoff stresses and deformation gradient
  // into a constitutive tensor relating second Piola Kirchhoff
  // stresses and Green-Lagrange strains cf.
  //
  // Marsden and Hughes, Mathematical Foundations of Elasticity,
  // Dover, pg. 215.

  cmat(0,0) = ((*cmatpf(0))[0]-S[0])*F_inv(0,0)*F_inv(0,0)
                 +(*cmatpf(7))[0]*F_inv(0,0)*F_inv(0,1)
                 +(*cmatpf(5))[0]*F_inv(0,0)*F_inv(0,2)
                 +(*cmatpf(0))[7]*F_inv(0,1)*F_inv(0,0)
                +((*cmatpf(7))[7]-S[0])*F_inv(0,1)*F_inv(0,1)
                 +(*cmatpf(5))[7]*F_inv(0,1)*F_inv(0,2)
                 +(*cmatpf(0))[5]*F_inv(0,2)*F_inv(0,0)
                 +(*cmatpf(7))[5]*F_inv(0,2)*F_inv(0,1)
                +((*cmatpf(5))[5]-S[0])*F_inv(0,2)*F_inv(0,2);

  cmat(0,1) = ((*cmatpf(3))[0]-S[3])*F_inv(0,0)*F_inv(1,0)
                 +(*cmatpf(1))[0]*F_inv(0,0)*F_inv(1,1)
                 +(*cmatpf(8))[0]*F_inv(0,0)*F_inv(1,2)
                 +(*cmatpf(3))[7]*F_inv(0,1)*F_inv(1,0)
                +((*cmatpf(1))[7]-S[3])*F_inv(0,1)*F_inv(1,1)
                 +(*cmatpf(8))[7]*F_inv(0,1)*F_inv(1,2)
                 +(*cmatpf(3))[5]*F_inv(0,2)*F_inv(1,0)
                 +(*cmatpf(1))[5]*F_inv(0,2)*F_inv(1,1)
                +((*cmatpf(8))[5]-S[3])*F_inv(0,2)*F_inv(1,2);

  cmat(0,2) = ((*cmatpf(6))[0]-S[5])*F_inv(0,0)*F_inv(2,0)
                 +(*cmatpf(4))[0]*F_inv(0,0)*F_inv(2,1)
                 +(*cmatpf(2))[0]*F_inv(0,0)*F_inv(2,2)
                 +(*cmatpf(6))[7]*F_inv(0,1)*F_inv(2,0)
                +((*cmatpf(4))[7]-S[5])*F_inv(0,1)*F_inv(2,1)
                 +(*cmatpf(2))[7]*F_inv(0,1)*F_inv(2,2)
                 +(*cmatpf(6))[5]*F_inv(0,2)*F_inv(2,0)
                 +(*cmatpf(4))[5]*F_inv(0,2)*F_inv(2,1)
                +((*cmatpf(2))[5]-S[5])*F_inv(0,2)*F_inv(2,2);

  cmat(0,3) = ((*cmatpf(3))[0]-S[3])*F_inv(0,0)*F_inv(0,0)
                 +(*cmatpf(1))[0]*F_inv(0,0)*F_inv(0,1)
                 +(*cmatpf(8))[0]*F_inv(0,0)*F_inv(0,2)
                 +(*cmatpf(3))[7]*F_inv(0,1)*F_inv(0,0)
                +((*cmatpf(1))[7]-S[3])*F_inv(0,1)*F_inv(0,1)
                 +(*cmatpf(8))[7]*F_inv(0,1)*F_inv(0,2)
                 +(*cmatpf(3))[5]*F_inv(0,2)*F_inv(0,0)
                 +(*cmatpf(1))[5]*F_inv(0,2)*F_inv(0,1)
                +((*cmatpf(8))[5]-S[3])*F_inv(0,2)*F_inv(0,2);

  cmat(0,4) = ((*cmatpf(6))[0]-S[5])*F_inv(0,0)*F_inv(1,0)
                 +(*cmatpf(4))[0]*F_inv(0,0)*F_inv(1,1)
                 +(*cmatpf(2))[0]*F_inv(0,0)*F_inv(1,2)
                 +(*cmatpf(6))[7]*F_inv(0,1)*F_inv(1,0)
                +((*cmatpf(4))[7]-S[5])*F_inv(0,1)*F_inv(1,1)
                 +(*cmatpf(2))[7]*F_inv(0,1)*F_inv(1,2)
                 +(*cmatpf(6))[5]*F_inv(0,2)*F_inv(1,0)
                 +(*cmatpf(4))[5]*F_inv(0,2)*F_inv(1,1)
                +((*cmatpf(2))[5]-S[5])*F_inv(0,2)*F_inv(1,2);

  cmat(0,5) = ((*cmatpf(0))[0]-S[0])*F_inv(0,0)*F_inv(2,0)
                 +(*cmatpf(7))[0]*F_inv(0,0)*F_inv(2,1)
                 +(*cmatpf(5))[0]*F_inv(0,0)*F_inv(2,2)
                 +(*cmatpf(0))[7]*F_inv(0,1)*F_inv(2,0)
                +((*cmatpf(7))[7]-S[0])*F_inv(0,1)*F_inv(2,1)
                 +(*cmatpf(5))[7]*F_inv(0,1)*F_inv(2,2)
                 +(*cmatpf(0))[5]*F_inv(0,2)*F_inv(2,0)
                 +(*cmatpf(7))[5]*F_inv(0,2)*F_inv(2,1)
                +((*cmatpf(5))[5]-S[0])*F_inv(0,2)*F_inv(2,2);

  cmat(1,1) = ((*cmatpf(3))[3]-S[1])*F_inv(1,0)*F_inv(1,0)
                 +(*cmatpf(1))[3]*F_inv(1,0)*F_inv(1,1)
                 +(*cmatpf(8))[3]*F_inv(1,0)*F_inv(1,2)
                 +(*cmatpf(3))[1]*F_inv(1,1)*F_inv(1,0)
                +((*cmatpf(1))[1]-S[1])*F_inv(1,1)*F_inv(1,1)
                 +(*cmatpf(8))[1]*F_inv(1,1)*F_inv(1,2)
                 +(*cmatpf(3))[8]*F_inv(1,2)*F_inv(1,0)
                 +(*cmatpf(1))[8]*F_inv(1,2)*F_inv(1,1)
                +((*cmatpf(8))[8]-S[1])*F_inv(1,2)*F_inv(1,2);

  cmat(1,2) = ((*cmatpf(6))[3]-S[4])*F_inv(1,0)*F_inv(2,0)
                 +(*cmatpf(4))[3]*F_inv(1,0)*F_inv(2,1)
                 +(*cmatpf(2))[3]*F_inv(1,0)*F_inv(2,2)
                 +(*cmatpf(6))[1]*F_inv(1,1)*F_inv(2,0)
                +((*cmatpf(4))[1]-S[4])*F_inv(1,1)*F_inv(2,1)
                 +(*cmatpf(2))[1]*F_inv(1,1)*F_inv(2,2)
                 +(*cmatpf(6))[8]*F_inv(1,2)*F_inv(2,0)
                 +(*cmatpf(4))[8]*F_inv(1,2)*F_inv(2,1)
                +((*cmatpf(2))[8]-S[4])*F_inv(1,2)*F_inv(2,2);

  cmat(1,3) = ((*cmatpf(3))[3]-S[1])*F_inv(1,0)*F_inv(0,0)
                 +(*cmatpf(1))[3]*F_inv(1,0)*F_inv(0,1)
                 +(*cmatpf(8))[3]*F_inv(1,0)*F_inv(0,2)
                 +(*cmatpf(3))[1]*F_inv(1,1)*F_inv(0,0)
                +((*cmatpf(1))[1]-S[1])*F_inv(1,1)*F_inv(0,1)
                 +(*cmatpf(8))[1]*F_inv(1,1)*F_inv(0,2)
                 +(*cmatpf(3))[8]*F_inv(1,2)*F_inv(0,0)
                 +(*cmatpf(1))[8]*F_inv(1,2)*F_inv(0,1)
                +((*cmatpf(8))[8]-S[1])*F_inv(1,2)*F_inv(0,2);

  cmat(1,4) = ((*cmatpf(6))[3]-S[4])*F_inv(1,0)*F_inv(1,0)
                 +(*cmatpf(4))[3]*F_inv(1,0)*F_inv(1,1)
                 +(*cmatpf(2))[3]*F_inv(1,0)*F_inv(1,2)
                 +(*cmatpf(6))[1]*F_inv(1,1)*F_inv(1,0)
                +((*cmatpf(4))[1]-S[4])*F_inv(1,1)*F_inv(1,1)
                 +(*cmatpf(2))[1]*F_inv(1,1)*F_inv(1,2)
                 +(*cmatpf(6))[8]*F_inv(1,2)*F_inv(1,0)
                 +(*cmatpf(4))[8]*F_inv(1,2)*F_inv(1,1)
                +((*cmatpf(2))[8]-S[4])*F_inv(1,2)*F_inv(1,2);

  cmat(1,5) = ((*cmatpf(0))[3]-S[3])*F_inv(1,0)*F_inv(2,0)
                 +(*cmatpf(7))[3]*F_inv(1,0)*F_inv(2,1)
                 +(*cmatpf(5))[3]*F_inv(1,0)*F_inv(2,2)
                 +(*cmatpf(0))[1]*F_inv(1,1)*F_inv(2,0)
                +((*cmatpf(7))[1]-S[3])*F_inv(1,1)*F_inv(2,1)
                 +(*cmatpf(5))[1]*F_inv(1,1)*F_inv(2,2)
                 +(*cmatpf(0))[8]*F_inv(1,2)*F_inv(2,0)
                 +(*cmatpf(7))[8]*F_inv(1,2)*F_inv(2,1)
                +((*cmatpf(5))[8]-S[3])*F_inv(1,2)*F_inv(2,2);

  cmat(2,2) = ((*cmatpf(6))[6]-S[2])*F_inv(2,0)*F_inv(2,0)
                 +(*cmatpf(4))[6]*F_inv(2,0)*F_inv(2,1)
                 +(*cmatpf(2))[6]*F_inv(2,0)*F_inv(2,2)
                 +(*cmatpf(6))[4]*F_inv(2,1)*F_inv(2,0)
                +((*cmatpf(4))[4]-S[2])*F_inv(2,1)*F_inv(2,1)
                 +(*cmatpf(2))[4]*F_inv(2,1)*F_inv(2,2)
                 +(*cmatpf(6))[2]*F_inv(2,2)*F_inv(2,0)
                 +(*cmatpf(4))[2]*F_inv(2,2)*F_inv(2,1)
                +((*cmatpf(2))[2]-S[2])*F_inv(2,2)*F_inv(2,2);

  cmat(2,3) = ((*cmatpf(3))[6]-S[4])*F_inv(2,0)*F_inv(0,0)
                 +(*cmatpf(1))[6]*F_inv(2,0)*F_inv(0,1)
                 +(*cmatpf(8))[6]*F_inv(2,0)*F_inv(0,2)
                 +(*cmatpf(3))[4]*F_inv(2,1)*F_inv(0,0)
                +((*cmatpf(1))[4]-S[4])*F_inv(2,1)*F_inv(0,1)
                 +(*cmatpf(8))[4]*F_inv(2,1)*F_inv(0,2)
                 +(*cmatpf(3))[2]*F_inv(2,2)*F_inv(0,0)
                 +(*cmatpf(1))[2]*F_inv(2,2)*F_inv(0,1)
                +((*cmatpf(8))[2]-S[4])*F_inv(2,2)*F_inv(0,2);

  cmat(2,4) = ((*cmatpf(6))[6]-S[2])*F_inv(2,0)*F_inv(1,0)
                 +(*cmatpf(4))[6]*F_inv(2,0)*F_inv(1,1)
                 +(*cmatpf(2))[6]*F_inv(2,0)*F_inv(1,2)
                 +(*cmatpf(6))[4]*F_inv(2,1)*F_inv(1,0)
                +((*cmatpf(4))[4]-S[2])*F_inv(2,1)*F_inv(1,1)
                 +(*cmatpf(2))[4]*F_inv(2,1)*F_inv(1,2)
                 +(*cmatpf(6))[2]*F_inv(2,2)*F_inv(1,0)
                 +(*cmatpf(4))[2]*F_inv(2,2)*F_inv(1,1)
                +((*cmatpf(2))[2]-S[2])*F_inv(2,2)*F_inv(1,2);

  cmat(2,5) = ((*cmatpf(0))[6]-S[5])*F_inv(2,0)*F_inv(2,0)
                 +(*cmatpf(7))[6]*F_inv(2,0)*F_inv(2,1)
                 +(*cmatpf(5))[6]*F_inv(2,0)*F_inv(2,2)
                 +(*cmatpf(0))[4]*F_inv(2,1)*F_inv(2,0)
                +((*cmatpf(7))[4]-S[5])*F_inv(2,1)*F_inv(2,1)
                 +(*cmatpf(5))[4]*F_inv(2,1)*F_inv(2,2)
                 +(*cmatpf(0))[2]*F_inv(2,2)*F_inv(2,0)
                 +(*cmatpf(7))[2]*F_inv(2,2)*F_inv(2,1)
                +((*cmatpf(5))[2]-S[5])*F_inv(2,2)*F_inv(2,2);

  cmat(3,3) = ((*cmatpf(3))[3]-S[1])*F_inv(0,0)*F_inv(0,0)
                 +(*cmatpf(1))[3]*F_inv(0,0)*F_inv(0,1)
                 +(*cmatpf(8))[3]*F_inv(0,0)*F_inv(0,2)
                 +(*cmatpf(3))[1]*F_inv(0,1)*F_inv(0,0)
                +((*cmatpf(1))[1]-S[1])*F_inv(0,1)*F_inv(0,1)
                 +(*cmatpf(8))[1]*F_inv(0,1)*F_inv(0,2)
                 +(*cmatpf(3))[8]*F_inv(0,2)*F_inv(0,0)
                 +(*cmatpf(1))[8]*F_inv(0,2)*F_inv(0,1)
                +((*cmatpf(8))[8]-S[1])*F_inv(0,2)*F_inv(0,2);

  cmat(3,4) = ((*cmatpf(6))[3]-S[4])*F_inv(0,0)*F_inv(1,0)
                 +(*cmatpf(4))[3]*F_inv(0,0)*F_inv(1,1)
                 +(*cmatpf(2))[3]*F_inv(0,0)*F_inv(1,2)
                 +(*cmatpf(6))[1]*F_inv(0,1)*F_inv(1,0)
                +((*cmatpf(4))[1]-S[4])*F_inv(0,1)*F_inv(1,1)
                 +(*cmatpf(2))[1]*F_inv(0,1)*F_inv(1,2)
                 +(*cmatpf(6))[8]*F_inv(0,2)*F_inv(1,0)
                 +(*cmatpf(4))[8]*F_inv(0,2)*F_inv(1,1)
                +((*cmatpf(2))[8]-S[4])*F_inv(0,2)*F_inv(1,2);

  cmat(3,5) = ((*cmatpf(0))[3]-S[3])*F_inv(0,0)*F_inv(2,0)
                 +(*cmatpf(7))[3]*F_inv(0,0)*F_inv(2,1)
                 +(*cmatpf(5))[3]*F_inv(0,0)*F_inv(2,2)
                 +(*cmatpf(0))[1]*F_inv(0,1)*F_inv(2,0)
                +((*cmatpf(7))[1]-S[3])*F_inv(0,1)*F_inv(2,1)
                 +(*cmatpf(5))[1]*F_inv(0,1)*F_inv(2,2)
                 +(*cmatpf(0))[8]*F_inv(0,2)*F_inv(2,0)
                 +(*cmatpf(7))[8]*F_inv(0,2)*F_inv(2,1)
                +((*cmatpf(5))[8]-S[3])*F_inv(0,2)*F_inv(2,2);

  cmat(4,4) = ((*cmatpf(6))[6]-S[2])*F_inv(1,0)*F_inv(1,0)
                 +(*cmatpf(4))[6]*F_inv(1,0)*F_inv(1,1)
                 +(*cmatpf(2))[6]*F_inv(1,0)*F_inv(1,2)
                 +(*cmatpf(6))[4]*F_inv(1,1)*F_inv(1,0)
                +((*cmatpf(4))[4]-S[2])*F_inv(1,1)*F_inv(1,1)
                 +(*cmatpf(2))[4]*F_inv(1,1)*F_inv(1,2)
                 +(*cmatpf(6))[2]*F_inv(1,2)*F_inv(1,0)
                 +(*cmatpf(4))[2]*F_inv(1,2)*F_inv(1,1)
                +((*cmatpf(2))[2]-S[2])*F_inv(1,2)*F_inv(1,2);

  cmat(4,5) = ((*cmatpf(0))[6]-S[5])*F_inv(1,0)*F_inv(2,0)
                 +(*cmatpf(7))[6]*F_inv(1,0)*F_inv(2,1)
                 +(*cmatpf(5))[6]*F_inv(1,0)*F_inv(2,2)
                 +(*cmatpf(0))[4]*F_inv(1,1)*F_inv(2,0)
                +((*cmatpf(7))[4]-S[5])*F_inv(1,1)*F_inv(2,1)
                 +(*cmatpf(5))[4]*F_inv(1,1)*F_inv(2,2)
                 +(*cmatpf(0))[2]*F_inv(1,2)*F_inv(2,0)
                 +(*cmatpf(7))[2]*F_inv(1,2)*F_inv(2,1)
                +((*cmatpf(5))[2]-S[5])*F_inv(1,2)*F_inv(2,2);

  cmat(5,5) = ((*cmatpf(0))[0]-S[0])*F_inv(2,0)*F_inv(2,0)
                 +(*cmatpf(7))[0]*F_inv(2,0)*F_inv(2,1)
                 +(*cmatpf(5))[0]*F_inv(2,0)*F_inv(2,2)
                 +(*cmatpf(0))[7]*F_inv(2,1)*F_inv(2,0)
                +((*cmatpf(7))[7]-S[0])*F_inv(2,1)*F_inv(2,1)
                 +(*cmatpf(5))[7]*F_inv(2,1)*F_inv(2,2)
                 +(*cmatpf(0))[5]*F_inv(2,2)*F_inv(2,0)
                 +(*cmatpf(7))[5]*F_inv(2,2)*F_inv(2,1)
                +((*cmatpf(5))[5]-S[0])*F_inv(2,2)*F_inv(2,2);

  // now assign all the symmetric components

  cmat(1,0) = cmat(0,1);
  cmat(2,0) = cmat(0,2);
  cmat(2,1) = cmat(1,2);
  cmat(3,0) = cmat(0,3);
  cmat(3,1) = cmat(1,3);
  cmat(3,2) = cmat(2,3);
  cmat(4,0) = cmat(0,4);
  cmat(4,1) = cmat(1,4);
  cmat(4,2) = cmat(2,4);
  cmat(4,3) = cmat(3,4);
  cmat(5,0) = cmat(0,5);
  cmat(5,1) = cmat(1,5);
  cmat(5,2) = cmat(2,5);
  cmat(5,3) = cmat(3,5);
  cmat(5,4) = cmat(4,5);

  cout << "cmat:\n" << cmat << "\n";
}

#endif
