#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "microstrugenalpha.H"

void MicroStruGenAlpha::calc_cmat(const Epetra_MultiVector& K_M,
                                  const Epetra_SerialDenseMatrix& inv_defgrd,
                                  const Epetra_SerialDenseVector& S,
                                  Epetra_SerialDenseMatrix* cmat,
                                  const Epetra_SerialDenseMatrix* F)
{
  DOUBLE cmat_pf_temp[3][3][3][3];
  DOUBLE cmat_pf_temp2[3][3][3][3];
  DOUBLE cmat_pf[3][3][3][3];

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        for (int l=0;l<3;++l)
        {
          cmat_pf_temp[i][j][k][l] = 0.;

          for (int m=0;m<np_/3;++m)
          {
            for (int n=0;n<np_/3;++n)
            {
              const Epetra_Vector* temp=K_M(3*n+k);  // this is a pointer to
                                                     // the corresponding
                                                     // vector (-> column!)

              cmat_pf_temp[i][j][k][l]+=(*Xp_)[3*m+j]*(*temp)[3*m+i]*(*Xp_)[3*n+l];
            }
          }
        }
      }
    }
  }

  // left conjugation of constitutive tensor (cf. Kouznetsova, thesis)
  // (T^ijkl)LC = T^jikl (change in first two indices)

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        for (int l=0;l<3;++l)
        {
          cmat_pf_temp2[i][j][k][l] = 1/V0_*cmat_pf_temp[j][i][k][l];
        }
      }
    }
  }

  // note that deformation gradient for Kouznetsova's approach is
  // conjugated, therefore the latter two indices need to be changed,
  // too (-> check out Miehe's approach)

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<3;++j)
    {
      for (int k=0;k<3;++k)
      {
        for (int l=0;l<3;++l)
        {
          cmat_pf[i][j][k][l] = cmat_pf_temp2[i][j][l][k];
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

//   double Emod = 100.;
//   double nu = 0.2;
//   double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
//   /* write non-zero components */
//   (*cmat)(0,0) = mfac*(1.0-nu);
//   (*cmat)(0,1) = mfac*nu;
//   (*cmat)(0,2) = mfac*nu;
//   (*cmat)(1,0) = mfac*nu;
//   (*cmat)(1,1) = mfac*(1.0-nu);
//   (*cmat)(1,2) = mfac*nu;
//   (*cmat)(2,0) = mfac*nu;
//   (*cmat)(2,1) = mfac*nu;
//   (*cmat)(2,2) = mfac*(1.0-nu);
//   /* ~~~ */
//   (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
//   (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
//   (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);

//   cout << "cmat:\n" << *cmat << endl;
//   exit(0);

//   FILE* output_fileptr;
//   output_fileptr=fopen("shear_cmat_rho1.0.txt","a");
//   fprintf(output_fileptr, "stress:\n");
//   for (int i=0;i<6;++i)
//     fprintf(output_fileptr, "%lf\n", S[i]);
//   fprintf(output_fileptr, "\n");
//   fprintf(output_fileptr, "cmat:\n");
//   for (int i=0;i<6;++i)
//   {
//     for (int j=0;j<6;++j)
//       fprintf(output_fileptr, "%lf\t", (*cmat)(i, j));
//     fprintf(output_fileptr, "\n");
//   }
//   fprintf(output_fileptr, "\n\n");
//   fclose(output_fileptr);
}

#endif
#endif
