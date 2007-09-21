#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "microstrugenalpha.H"

void MicroStruGenAlpha::calc_cmat(const Epetra_MultiVector& K_M,
                                  const Epetra_SerialDenseMatrix& inv_defgrd,
                                  const Epetra_SerialDenseVector& S,
                                  Epetra_SerialDenseMatrix* cmat)
{
  DOUBLE cmat_pf[3][3][3][3];

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        for (int l=0;l<3;++l)
        {
          cmat_pf[i][j][k][l] = 0.;

          for (int m=0;m<np_/3;++m)
          {
            for (int n=0;n<np_/3;++n)
            {
              cmat_pf[i][j][k][l] += (*Xp_)[3*m+i]*K_M[3*m+j][3*n+k]*(*Xp_)[3*n+l];
            }
          }

        }
      }
    }
  }

  // definition of Kronecker delta

  Epetra_SerialDenseMatrix delta(3,3);

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      delta(i,j) = 0.;
    }
  }
  delta(0,0) = 1.;
  delta(1,1) = 1.;
  delta(2,2) = 1.;


  for (int i=0;i<6;++i)
  {
    for (int j=0;j<6;++j)
    {
      (*cmat)(i,j) = 0.;
    }
  }

  // take advantage of symmetry of the constitutive tensor relating
  // Green-Lagrange strains and second Piola-Kirchhoff stresses

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      (*cmat)(0,0) += (cmat_pf[i][0][j][0]-S[0]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(0,j);
      (*cmat)(0,1) += (cmat_pf[i][0][j][1]-S[3]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(1,j);
      (*cmat)(0,2) += (cmat_pf[i][0][j][2]-S[5]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(2,j);
      (*cmat)(0,3) += (cmat_pf[i][0][j][1]-S[3]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(0,j);
      (*cmat)(0,4) += (cmat_pf[i][0][j][2]-S[5]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(1,j);
      (*cmat)(0,5) += (cmat_pf[i][0][j][0]-S[0]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(2,j);

      (*cmat)(1,1) += (cmat_pf[i][1][j][1]-S[1]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(1,j);
      (*cmat)(1,2) += (cmat_pf[i][1][j][2]-S[4]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(2,j);
      (*cmat)(1,3) += (cmat_pf[i][1][j][1]-S[1]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(0,j);
      (*cmat)(1,4) += (cmat_pf[i][1][j][2]-S[4]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(1,j);
      (*cmat)(1,5) += (cmat_pf[i][1][j][0]-S[3]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(2,j);

      (*cmat)(2,2) += (cmat_pf[i][2][j][2]-S[2]*delta(i,j))*inv_defgrd(2,i)*inv_defgrd(2,j);
      (*cmat)(2,3) += (cmat_pf[i][2][j][1]-S[4]*delta(i,j))*inv_defgrd(2,i)*inv_defgrd(0,j);
      (*cmat)(2,4) += (cmat_pf[i][2][j][2]-S[2]*delta(i,j))*inv_defgrd(2,i)*inv_defgrd(1,j);
      (*cmat)(2,5) += (cmat_pf[i][2][j][0]-S[5]*delta(i,j))*inv_defgrd(2,i)*inv_defgrd(2,j);

      (*cmat)(3,3) += (cmat_pf[i][1][j][1]-S[1]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(0,j);
      (*cmat)(3,4) += (cmat_pf[i][1][j][2]-S[4]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(1,j);
      (*cmat)(3,5) += (cmat_pf[i][1][j][0]-S[3]*delta(i,j))*inv_defgrd(0,i)*inv_defgrd(2,j);

      (*cmat)(4,4) += (cmat_pf[i][2][j][2]-S[2]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(1,j);
      (*cmat)(4,5) += (cmat_pf[i][2][j][0]-S[5]*delta(i,j))*inv_defgrd(1,i)*inv_defgrd(2,j);

      (*cmat)(5,5) += (cmat_pf[i][0][j][0]-S[0]*delta(i,j))*inv_defgrd(2,i)*inv_defgrd(2,j);
    }
  }

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
#endif
