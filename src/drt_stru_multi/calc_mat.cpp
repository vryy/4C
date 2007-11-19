#ifdef CCADISCRET

#include "microstrugenalpha.H"

void MicroStruGenAlpha::calc_cmat(const Epetra_MultiVector& K_M,
                                  const Epetra_SerialDenseMatrix& inv_defgrd,
                                  const Epetra_SerialDenseVector& S,
                                  Epetra_SerialDenseMatrix* cmat,
                                  const Epetra_SerialDenseMatrix* F)
{
  // The calculation of the consistent macroscopic constitutive tensor
  // follows
  //
  // C. Miehe, Computational micro-to-macro transitions for
  // discretized micro-structures of heterogeneous materials at finite
  // strains based on a minimization of averaged incremental energy.
  // Computer Methods in Applied Mechanics and Engineering 192: 559-591, 2003.

  Epetra_Map Dmap(9, 0, Epetra_SerialComm());
  Epetra_MultiVector D(Dmap, np_);
  for (int n=0;n<np_/3;++n)
  {
    Epetra_Vector* temp1 = D(3*n);
    (*temp1)[0] = (*Xp_)[3*n];
    (*temp1)[3] = (*Xp_)[3*n+1];
    (*temp1)[6] = (*Xp_)[3*n+2];
    Epetra_Vector* temp2 = D(3*n+1);
    (*temp2)[1] = (*Xp_)[3*n+1];
    (*temp2)[4] = (*Xp_)[3*n+2];
    (*temp2)[7] = (*Xp_)[3*n];
    Epetra_Vector* temp3 = D(3*n+2);
    (*temp3)[2] = (*Xp_)[3*n+2];
    (*temp3)[5] = (*Xp_)[3*n];
    (*temp3)[8] = (*Xp_)[3*n+1];
  }

  Epetra_MultiVector cmatpf_temp(Dmap, np_);
  cmatpf_temp.Multiply('N', 'N', 1.0, D, K_M, 0.0);
  Epetra_MultiVector cmatpf(Dmap, 9);
  cmatpf.Multiply('N', 'T', 1.0, cmatpf_temp, D, 0.0);

  for (int i=0;i<9;++i)
  {
    for (int j=0;j<9;++j)
    {
      (*cmatpf(j))[i] /= V0_;
    }
  }

  // We now have to transform the calculated constitutive tensor
  // relating first Piola-Kirchhoff stresses to the deformation
  // gradient into a constitutive tensor relating second
  // Piola-Kirchhoff stresses to Green-Lagrange strains. For
  // conversion details cf.
  //
  // Marsden and Hughes, Mathematical Foundations of Elasticity,
  // Dover, pg. 215.

  (*cmat)(0,0) = ((*cmatpf(0))[0]-S[0])*inv_defgrd(0,0)*inv_defgrd(0,0)
                 +(*cmatpf(7))[0]*inv_defgrd(0,0)*inv_defgrd(0,1)
                 +(*cmatpf(5))[0]*inv_defgrd(0,0)*inv_defgrd(0,2)
                 +(*cmatpf(0))[7]*inv_defgrd(0,1)*inv_defgrd(0,0)
                +((*cmatpf(7))[7]-S[0])*inv_defgrd(0,1)*inv_defgrd(0,1)
                 +(*cmatpf(5))[7]*inv_defgrd(0,1)*inv_defgrd(0,2)
                 +(*cmatpf(0))[5]*inv_defgrd(0,2)*inv_defgrd(0,0)
                 +(*cmatpf(7))[5]*inv_defgrd(0,2)*inv_defgrd(0,1)
                +((*cmatpf(5))[5]-S[0])*inv_defgrd(0,2)*inv_defgrd(0,2);

  (*cmat)(0,1) = ((*cmatpf(3))[0]-S[3])*inv_defgrd(0,0)*inv_defgrd(1,0)
                 +(*cmatpf(1))[0]*inv_defgrd(0,0)*inv_defgrd(1,1)
                 +(*cmatpf(8))[0]*inv_defgrd(0,0)*inv_defgrd(1,2)
                 +(*cmatpf(3))[7]*inv_defgrd(0,1)*inv_defgrd(1,0)
                +((*cmatpf(1))[7]-S[3])*inv_defgrd(0,1)*inv_defgrd(1,1)
                 +(*cmatpf(8))[7]*inv_defgrd(0,1)*inv_defgrd(1,2)
                 +(*cmatpf(3))[5]*inv_defgrd(0,2)*inv_defgrd(1,0)
                 +(*cmatpf(1))[5]*inv_defgrd(0,2)*inv_defgrd(1,1)
                +((*cmatpf(8))[5]-S[3])*inv_defgrd(0,2)*inv_defgrd(1,2);

  (*cmat)(0,2) = ((*cmatpf(6))[0]-S[5])*inv_defgrd(0,0)*inv_defgrd(2,0)
                 +(*cmatpf(4))[0]*inv_defgrd(0,0)*inv_defgrd(2,1)
                 +(*cmatpf(2))[0]*inv_defgrd(0,0)*inv_defgrd(2,2)
                 +(*cmatpf(6))[7]*inv_defgrd(0,1)*inv_defgrd(2,0)
                +((*cmatpf(4))[7]-S[5])*inv_defgrd(0,1)*inv_defgrd(2,1)
                 +(*cmatpf(2))[7]*inv_defgrd(0,1)*inv_defgrd(2,2)
                 +(*cmatpf(6))[5]*inv_defgrd(0,2)*inv_defgrd(2,0)
                 +(*cmatpf(4))[5]*inv_defgrd(0,2)*inv_defgrd(2,1)
                +((*cmatpf(2))[5]-S[5])*inv_defgrd(0,2)*inv_defgrd(2,2);

  (*cmat)(0,3) = ((*cmatpf(3))[0]-S[3])*inv_defgrd(0,0)*inv_defgrd(0,0)
                 +(*cmatpf(1))[0]*inv_defgrd(0,0)*inv_defgrd(0,1)
                 +(*cmatpf(8))[0]*inv_defgrd(0,0)*inv_defgrd(0,2)
                 +(*cmatpf(3))[7]*inv_defgrd(0,1)*inv_defgrd(0,0)
                +((*cmatpf(1))[7]-S[3])*inv_defgrd(0,1)*inv_defgrd(0,1)
                 +(*cmatpf(8))[7]*inv_defgrd(0,1)*inv_defgrd(0,2)
                 +(*cmatpf(3))[5]*inv_defgrd(0,2)*inv_defgrd(0,0)
                 +(*cmatpf(1))[5]*inv_defgrd(0,2)*inv_defgrd(0,1)
                +((*cmatpf(8))[5]-S[3])*inv_defgrd(0,2)*inv_defgrd(0,2);

  (*cmat)(0,4) = ((*cmatpf(6))[0]-S[5])*inv_defgrd(0,0)*inv_defgrd(1,0)
                 +(*cmatpf(4))[0]*inv_defgrd(0,0)*inv_defgrd(1,1)
                 +(*cmatpf(2))[0]*inv_defgrd(0,0)*inv_defgrd(1,2)
                 +(*cmatpf(6))[7]*inv_defgrd(0,1)*inv_defgrd(1,0)
                +((*cmatpf(4))[7]-S[5])*inv_defgrd(0,1)*inv_defgrd(1,1)
                 +(*cmatpf(2))[7]*inv_defgrd(0,1)*inv_defgrd(1,2)
                 +(*cmatpf(6))[5]*inv_defgrd(0,2)*inv_defgrd(1,0)
                 +(*cmatpf(4))[5]*inv_defgrd(0,2)*inv_defgrd(1,1)
                +((*cmatpf(2))[5]-S[5])*inv_defgrd(0,2)*inv_defgrd(1,2);

  (*cmat)(0,5) = ((*cmatpf(0))[0]-S[0])*inv_defgrd(0,0)*inv_defgrd(2,0)
                 +(*cmatpf(7))[0]*inv_defgrd(0,0)*inv_defgrd(2,1)
                 +(*cmatpf(5))[0]*inv_defgrd(0,0)*inv_defgrd(2,2)
                 +(*cmatpf(0))[7]*inv_defgrd(0,1)*inv_defgrd(2,0)
                +((*cmatpf(7))[7]-S[0])*inv_defgrd(0,1)*inv_defgrd(2,1)
                 +(*cmatpf(5))[7]*inv_defgrd(0,1)*inv_defgrd(2,2)
                 +(*cmatpf(0))[5]*inv_defgrd(0,2)*inv_defgrd(2,0)
                 +(*cmatpf(7))[5]*inv_defgrd(0,2)*inv_defgrd(2,1)
                +((*cmatpf(5))[5]-S[0])*inv_defgrd(0,2)*inv_defgrd(2,2);

  (*cmat)(1,1) = ((*cmatpf(3))[3]-S[1])*inv_defgrd(1,0)*inv_defgrd(1,0)
                 +(*cmatpf(1))[3]*inv_defgrd(1,0)*inv_defgrd(1,1)
                 +(*cmatpf(8))[3]*inv_defgrd(1,0)*inv_defgrd(1,2)
                 +(*cmatpf(3))[1]*inv_defgrd(1,1)*inv_defgrd(1,0)
                +((*cmatpf(1))[1]-S[1])*inv_defgrd(1,1)*inv_defgrd(1,1)
                 +(*cmatpf(8))[1]*inv_defgrd(1,1)*inv_defgrd(1,2)
                 +(*cmatpf(3))[8]*inv_defgrd(1,2)*inv_defgrd(1,0)
                 +(*cmatpf(1))[8]*inv_defgrd(1,2)*inv_defgrd(1,1)
                +((*cmatpf(8))[8]-S[1])*inv_defgrd(1,2)*inv_defgrd(1,2);

  (*cmat)(1,2) = ((*cmatpf(6))[3]-S[4])*inv_defgrd(1,0)*inv_defgrd(2,0)
                 +(*cmatpf(4))[3]*inv_defgrd(1,0)*inv_defgrd(2,1)
                 +(*cmatpf(2))[3]*inv_defgrd(1,0)*inv_defgrd(2,2)
                 +(*cmatpf(6))[1]*inv_defgrd(1,1)*inv_defgrd(2,0)
                +((*cmatpf(4))[1]-S[4])*inv_defgrd(1,1)*inv_defgrd(2,1)
                 +(*cmatpf(2))[1]*inv_defgrd(1,1)*inv_defgrd(2,2)
                 +(*cmatpf(6))[8]*inv_defgrd(1,2)*inv_defgrd(2,0)
                 +(*cmatpf(4))[8]*inv_defgrd(1,2)*inv_defgrd(2,1)
                +((*cmatpf(2))[8]-S[4])*inv_defgrd(1,2)*inv_defgrd(2,2);

  (*cmat)(1,3) = ((*cmatpf(3))[3]-S[1])*inv_defgrd(1,0)*inv_defgrd(0,0)
                 +(*cmatpf(1))[3]*inv_defgrd(1,0)*inv_defgrd(0,1)
                 +(*cmatpf(8))[3]*inv_defgrd(1,0)*inv_defgrd(0,2)
                 +(*cmatpf(3))[1]*inv_defgrd(1,1)*inv_defgrd(0,0)
                +((*cmatpf(1))[1]-S[1])*inv_defgrd(1,1)*inv_defgrd(0,1)
                 +(*cmatpf(8))[1]*inv_defgrd(1,1)*inv_defgrd(0,2)
                 +(*cmatpf(3))[8]*inv_defgrd(1,2)*inv_defgrd(0,0)
                 +(*cmatpf(1))[8]*inv_defgrd(1,2)*inv_defgrd(0,1)
                +((*cmatpf(8))[8]-S[1])*inv_defgrd(1,2)*inv_defgrd(0,2);

  (*cmat)(1,4) = ((*cmatpf(6))[3]-S[4])*inv_defgrd(1,0)*inv_defgrd(1,0)
                 +(*cmatpf(4))[3]*inv_defgrd(1,0)*inv_defgrd(1,1)
                 +(*cmatpf(2))[3]*inv_defgrd(1,0)*inv_defgrd(1,2)
                 +(*cmatpf(6))[1]*inv_defgrd(1,1)*inv_defgrd(1,0)
                +((*cmatpf(4))[1]-S[4])*inv_defgrd(1,1)*inv_defgrd(1,1)
                 +(*cmatpf(2))[1]*inv_defgrd(1,1)*inv_defgrd(1,2)
                 +(*cmatpf(6))[8]*inv_defgrd(1,2)*inv_defgrd(1,0)
                 +(*cmatpf(4))[8]*inv_defgrd(1,2)*inv_defgrd(1,1)
                +((*cmatpf(2))[8]-S[4])*inv_defgrd(1,2)*inv_defgrd(1,2);

  (*cmat)(1,5) = ((*cmatpf(0))[3]-S[3])*inv_defgrd(1,0)*inv_defgrd(2,0)
                 +(*cmatpf(7))[3]*inv_defgrd(1,0)*inv_defgrd(2,1)
                 +(*cmatpf(5))[3]*inv_defgrd(1,0)*inv_defgrd(2,2)
                 +(*cmatpf(0))[1]*inv_defgrd(1,1)*inv_defgrd(2,0)
                +((*cmatpf(7))[1]-S[3])*inv_defgrd(1,1)*inv_defgrd(2,1)
                 +(*cmatpf(5))[1]*inv_defgrd(1,1)*inv_defgrd(2,2)
                 +(*cmatpf(0))[8]*inv_defgrd(1,2)*inv_defgrd(2,0)
                 +(*cmatpf(7))[8]*inv_defgrd(1,2)*inv_defgrd(2,1)
                +((*cmatpf(5))[8]-S[3])*inv_defgrd(1,2)*inv_defgrd(2,2);

  (*cmat)(2,2) = ((*cmatpf(6))[6]-S[2])*inv_defgrd(2,0)*inv_defgrd(2,0)
                 +(*cmatpf(4))[6]*inv_defgrd(2,0)*inv_defgrd(2,1)
                 +(*cmatpf(2))[6]*inv_defgrd(2,0)*inv_defgrd(2,2)
                 +(*cmatpf(6))[4]*inv_defgrd(2,1)*inv_defgrd(2,0)
                +((*cmatpf(4))[4]-S[2])*inv_defgrd(2,1)*inv_defgrd(2,1)
                 +(*cmatpf(2))[4]*inv_defgrd(2,1)*inv_defgrd(2,2)
                 +(*cmatpf(6))[2]*inv_defgrd(2,2)*inv_defgrd(2,0)
                 +(*cmatpf(4))[2]*inv_defgrd(2,2)*inv_defgrd(2,1)
                +((*cmatpf(2))[2]-S[2])*inv_defgrd(2,2)*inv_defgrd(2,2);

  (*cmat)(2,3) = ((*cmatpf(3))[6]-S[4])*inv_defgrd(2,0)*inv_defgrd(0,0)
                 +(*cmatpf(1))[6]*inv_defgrd(2,0)*inv_defgrd(0,1)
                 +(*cmatpf(8))[6]*inv_defgrd(2,0)*inv_defgrd(0,2)
                 +(*cmatpf(3))[4]*inv_defgrd(2,1)*inv_defgrd(0,0)
                +((*cmatpf(1))[4]-S[4])*inv_defgrd(2,1)*inv_defgrd(0,1)
                 +(*cmatpf(8))[4]*inv_defgrd(2,1)*inv_defgrd(0,2)
                 +(*cmatpf(3))[2]*inv_defgrd(2,2)*inv_defgrd(0,0)
                 +(*cmatpf(1))[2]*inv_defgrd(2,2)*inv_defgrd(0,1)
                +((*cmatpf(8))[2]-S[4])*inv_defgrd(2,2)*inv_defgrd(0,2);

  (*cmat)(2,4) = ((*cmatpf(6))[6]-S[2])*inv_defgrd(2,0)*inv_defgrd(1,0)
                 +(*cmatpf(4))[6]*inv_defgrd(2,0)*inv_defgrd(1,1)
                 +(*cmatpf(2))[6]*inv_defgrd(2,0)*inv_defgrd(1,2)
                 +(*cmatpf(6))[4]*inv_defgrd(2,1)*inv_defgrd(1,0)
                +((*cmatpf(4))[4]-S[2])*inv_defgrd(2,1)*inv_defgrd(1,1)
                 +(*cmatpf(2))[4]*inv_defgrd(2,1)*inv_defgrd(1,2)
                 +(*cmatpf(6))[2]*inv_defgrd(2,2)*inv_defgrd(1,0)
                 +(*cmatpf(4))[2]*inv_defgrd(2,2)*inv_defgrd(1,1)
                +((*cmatpf(2))[2]-S[2])*inv_defgrd(2,2)*inv_defgrd(1,2);

  (*cmat)(2,5) = ((*cmatpf(0))[6]-S[5])*inv_defgrd(2,0)*inv_defgrd(2,0)
                 +(*cmatpf(7))[6]*inv_defgrd(2,0)*inv_defgrd(2,1)
                 +(*cmatpf(5))[6]*inv_defgrd(2,0)*inv_defgrd(2,2)
                 +(*cmatpf(0))[4]*inv_defgrd(2,1)*inv_defgrd(2,0)
                +((*cmatpf(7))[4]-S[5])*inv_defgrd(2,1)*inv_defgrd(2,1)
                 +(*cmatpf(5))[4]*inv_defgrd(2,1)*inv_defgrd(2,2)
                 +(*cmatpf(0))[2]*inv_defgrd(2,2)*inv_defgrd(2,0)
                 +(*cmatpf(7))[2]*inv_defgrd(2,2)*inv_defgrd(2,1)
                +((*cmatpf(5))[2]-S[5])*inv_defgrd(2,2)*inv_defgrd(2,2);

  (*cmat)(3,3) = ((*cmatpf(3))[3]-S[1])*inv_defgrd(0,0)*inv_defgrd(0,0)
                 +(*cmatpf(1))[3]*inv_defgrd(0,0)*inv_defgrd(0,1)
                 +(*cmatpf(8))[3]*inv_defgrd(0,0)*inv_defgrd(0,2)
                 +(*cmatpf(3))[1]*inv_defgrd(0,1)*inv_defgrd(0,0)
                +((*cmatpf(1))[1]-S[1])*inv_defgrd(0,1)*inv_defgrd(0,1)
                 +(*cmatpf(8))[1]*inv_defgrd(0,1)*inv_defgrd(0,2)
                 +(*cmatpf(3))[8]*inv_defgrd(0,2)*inv_defgrd(0,0)
                 +(*cmatpf(1))[8]*inv_defgrd(0,2)*inv_defgrd(0,1)
                +((*cmatpf(8))[8]-S[1])*inv_defgrd(0,2)*inv_defgrd(0,2);

  (*cmat)(3,4) = ((*cmatpf(6))[3]-S[4])*inv_defgrd(0,0)*inv_defgrd(1,0)
                 +(*cmatpf(4))[3]*inv_defgrd(0,0)*inv_defgrd(1,1)
                 +(*cmatpf(2))[3]*inv_defgrd(0,0)*inv_defgrd(1,2)
                 +(*cmatpf(6))[1]*inv_defgrd(0,1)*inv_defgrd(1,0)
                +((*cmatpf(4))[1]-S[4])*inv_defgrd(0,1)*inv_defgrd(1,1)
                 +(*cmatpf(2))[1]*inv_defgrd(0,1)*inv_defgrd(1,2)
                 +(*cmatpf(6))[8]*inv_defgrd(0,2)*inv_defgrd(1,0)
                 +(*cmatpf(4))[8]*inv_defgrd(0,2)*inv_defgrd(1,1)
                +((*cmatpf(2))[8]-S[4])*inv_defgrd(0,2)*inv_defgrd(1,2);

  (*cmat)(3,5) = ((*cmatpf(0))[3]-S[3])*inv_defgrd(0,0)*inv_defgrd(2,0)
                 +(*cmatpf(7))[3]*inv_defgrd(0,0)*inv_defgrd(2,1)
                 +(*cmatpf(5))[3]*inv_defgrd(0,0)*inv_defgrd(2,2)
                 +(*cmatpf(0))[1]*inv_defgrd(0,1)*inv_defgrd(2,0)
                +((*cmatpf(7))[1]-S[3])*inv_defgrd(0,1)*inv_defgrd(2,1)
                 +(*cmatpf(5))[1]*inv_defgrd(0,1)*inv_defgrd(2,2)
                 +(*cmatpf(0))[8]*inv_defgrd(0,2)*inv_defgrd(2,0)
                 +(*cmatpf(7))[8]*inv_defgrd(0,2)*inv_defgrd(2,1)
                +((*cmatpf(5))[8]-S[3])*inv_defgrd(0,2)*inv_defgrd(2,2);

  (*cmat)(4,4) = ((*cmatpf(6))[6]-S[2])*inv_defgrd(1,0)*inv_defgrd(1,0)
                 +(*cmatpf(4))[6]*inv_defgrd(1,0)*inv_defgrd(1,1)
                 +(*cmatpf(2))[6]*inv_defgrd(1,0)*inv_defgrd(1,2)
                 +(*cmatpf(6))[4]*inv_defgrd(1,1)*inv_defgrd(1,0)
                +((*cmatpf(4))[4]-S[2])*inv_defgrd(1,1)*inv_defgrd(1,1)
                 +(*cmatpf(2))[4]*inv_defgrd(1,1)*inv_defgrd(1,2)
                 +(*cmatpf(6))[2]*inv_defgrd(1,2)*inv_defgrd(1,0)
                 +(*cmatpf(4))[2]*inv_defgrd(1,2)*inv_defgrd(1,1)
                +((*cmatpf(2))[2]-S[2])*inv_defgrd(1,2)*inv_defgrd(1,2);

  (*cmat)(4,5) = ((*cmatpf(0))[6]-S[5])*inv_defgrd(1,0)*inv_defgrd(2,0)
                 +(*cmatpf(7))[6]*inv_defgrd(1,0)*inv_defgrd(2,1)
                 +(*cmatpf(5))[6]*inv_defgrd(1,0)*inv_defgrd(2,2)
                 +(*cmatpf(0))[4]*inv_defgrd(1,1)*inv_defgrd(2,0)
                +((*cmatpf(7))[4]-S[5])*inv_defgrd(1,1)*inv_defgrd(2,1)
                 +(*cmatpf(5))[4]*inv_defgrd(1,1)*inv_defgrd(2,2)
                 +(*cmatpf(0))[2]*inv_defgrd(1,2)*inv_defgrd(2,0)
                 +(*cmatpf(7))[2]*inv_defgrd(1,2)*inv_defgrd(2,1)
                +((*cmatpf(5))[2]-S[5])*inv_defgrd(1,2)*inv_defgrd(2,2);

  (*cmat)(5,5) = ((*cmatpf(0))[0]-S[0])*inv_defgrd(2,0)*inv_defgrd(2,0)
                 +(*cmatpf(7))[0]*inv_defgrd(2,0)*inv_defgrd(2,1)
                 +(*cmatpf(5))[0]*inv_defgrd(2,0)*inv_defgrd(2,2)
                 +(*cmatpf(0))[7]*inv_defgrd(2,1)*inv_defgrd(2,0)
                +((*cmatpf(7))[7]-S[0])*inv_defgrd(2,1)*inv_defgrd(2,1)
                 +(*cmatpf(5))[7]*inv_defgrd(2,1)*inv_defgrd(2,2)
                 +(*cmatpf(0))[5]*inv_defgrd(2,2)*inv_defgrd(2,0)
                 +(*cmatpf(7))[5]*inv_defgrd(2,2)*inv_defgrd(2,1)
                +((*cmatpf(5))[5]-S[0])*inv_defgrd(2,2)*inv_defgrd(2,2);

  // now assign all the symmetric components

  (*cmat)(1,0) = (*cmat)(0,1);
  (*cmat)(2,0) = (*cmat)(0,2);
  (*cmat)(2,1) = (*cmat)(1,2);
  (*cmat)(3,0) = (*cmat)(0,3);
  (*cmat)(3,1) = (*cmat)(1,3);
  (*cmat)(3,2) = (*cmat)(2,3);
  (*cmat)(4,0) = (*cmat)(0,4);
  (*cmat)(4,1) = (*cmat)(1,4);
  (*cmat)(4,2) = (*cmat)(2,4);
  (*cmat)(4,3) = (*cmat)(3,4);
  (*cmat)(5,0) = (*cmat)(0,5);
  (*cmat)(5,1) = (*cmat)(1,5);
  (*cmat)(5,2) = (*cmat)(2,5);
  (*cmat)(5,3) = (*cmat)(3,5);
  (*cmat)(5,4) = (*cmat)(4,5);
}

#endif
