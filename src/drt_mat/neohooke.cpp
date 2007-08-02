/*!----------------------------------------------------------------------
\file neohooke.cpp
\brief contains the functions to establish local material law
       stress-strain law for isotropic material for a 3D hex element
       following compressible Neo Hookean material law
       see Holzapfel, Nonlinear Solid Mechanics, pp. 243
       example input line:
       MAT 1 MAT_Struct_NeoHooke  YOUNG 100.0 NUE 0.3 DENS 1.0
<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
\param  Epetra_SerialDenseVector* glstrain      (i) Green-Lagrange strains
\param  Epetra_SerialDenseVector* stress        (o) ele stress vector
\param  Epetra_SerialDenseMatrix* cmat          (o) constitutive matrix
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "neohooke.H"

extern struct _MATERIAL *mat;

MAT::NeoHooke::NeoHooke()
  : matdata_(NULL)
{
}


MAT::NeoHooke::NeoHooke(MATERIAL* matdata)
  : matdata_(matdata)
{
}

void MAT::NeoHooke::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
}


void MAT::NeoHooke::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


/*----------------------------------------------------------------------*
 |  Calculate stress and constitutive tensor (Neo Hookean law)  rm 08/07|
 *----------------------------------------------------------------------*/
void MAT::NeoHooke::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                      Epetra_SerialDenseMatrix* cmat,
                                      Epetra_SerialDenseVector* stress)
{
  // get material parameters
  double ym = matdata_->m.neohooke->youngs;    // Young's modulus
  double nu = matdata_->m.neohooke->possionratio; // Poisson's ratio

  // Identity Matrix
  Epetra_SerialDenseMatrix I(3,3);
  for (int i = 0; i < 3; ++i) I(i,i) = 1.0;
  
  // Green-Lagrange Strain Tensor
  Epetra_SerialDenseMatrix E(3,3);
  E(0,0) = (*glstrain)(0);
  E(1,1) = (*glstrain)(1);
  E(2,2) = (*glstrain)(2);
  E(0,1) = 0.5 * (*glstrain)(3);  E(1,0) = 0.5 * (*glstrain)(3);
  E(1,2) = 0.5 * (*glstrain)(4);  E(2,1) = 0.5 * (*glstrain)(4);
  E(0,2) = 0.5 * (*glstrain)(5);  E(2,0) = 0.5 * (*glstrain)(5);
  
  // Right Cauchy-Green Tensor  C = 2 * E + I
  Epetra_SerialDenseMatrix C(E);
  C.Scale(2.0);
  C += I;
   
  // Principal Invariants I1 = tr(C) and I3 = det(C)
  //double I1 = C(0,0)+C(1,1)+C(2,2); // Needed only for energy
  double I3 = C(0,0)*C(1,1)*C(2,2) + C(0,1)*C(1,2)*C(2,0)
            + C(0,2)*C(1,0)*C(2,1) - (C(0,2)*C(1,1)*C(2,0)
            + C(0,1)*C(1,0)*C(2,2) + C(0,0)*C(1,2)*C(2,1));
  
  // Calculation of C^-1 (Cinv)
  Epetra_SerialDenseMatrix Cinv = inverseTensor(C,I3);

  // Material Constants c1 and beta
  double c1 = 0.5 * ym/(2*(1+nu));
  double beta = nu/(1-2*nu);

  // Energy
  // double W = c1/beta * (pow(J,-beta) - 1) + c1 (I1-3);
  
  // PK2 Stresses
  Epetra_SerialDenseMatrix PK2(3,3);
  int i,j;  
  for (i=0; i<3; i++)
  for (j=0; j<3; j++)
  {
      PK2(i,j)= 2 * c1 * ( I(i,j) - pow(I3,-beta) * Cinv(i,j) );
  }
  
  // Transfer PK2 tensor to stress vector
  (*stress)(0) = PK2(0,0);
  (*stress)(1) = PK2(1,1);
  (*stress)(2) = PK2(2,2);
  (*stress)(3) = PK2(0,1);
  (*stress)(4) = PK2(1,2);
  (*stress)(5) = PK2(0,2);

  // Elasticity Tensor
  double delta6 = 4. * c1 * beta * pow(I3,-beta);
  double delta7 = 4. * c1 * pow(I3,-beta);
  
  int k,l;
  Epetra_SerialDenseMatrix ET(9,9);

  
  for (k=0; k<3; k++)
  for (l=0; l<3; l++)
  {
		ET[k][l]	= delta6 * (Cinv(0,0) * Cinv(k,l)) + 
					  delta7 * 0.5 *(Cinv(0,k) * Cinv(0,l) + Cinv(0,l) * Cinv(0,k));
		ET[k+3][l]	= delta6 * (Cinv(1,0) * Cinv(k,l)) + 
					  delta7 * 0.5 *(Cinv(1,k) * Cinv(0,l) + Cinv(1,l) * Cinv(0,k));
		ET[k+3][l+3]= delta6 * (Cinv(1,1) * Cinv(k,l)) +
				      delta7 * 0.5 *(Cinv(1,k) * Cinv(1,l) + Cinv(1,l) * Cinv(1,k));
		ET[k+6][l]	= delta6 * (Cinv(2,0) * Cinv(k,l)) +
					  delta7 * 0.5 *(Cinv(2,k) * Cinv(0,l) + Cinv(2,l) * Cinv(0,k));
		ET[k+6][l+3]= delta6 * (Cinv(2,1) * Cinv(k,l)) +
					  delta7 * 0.5 *(Cinv(2,k) * Cinv(1,l) + Cinv(2,l) * Cinv(1,k));
		ET[k+6][l+6]= delta6 * (Cinv(2,2) * Cinv(k,l)) +
					  delta7 * 0.5 *(Cinv(2,k) * Cinv(2,l) + Cinv(2,l) * Cinv(2,k));
  }

  (*cmat)(0,0)=ET[0][0];
  (*cmat)(0,1)=ET[1][1];
  (*cmat)(0,2)=ET[2][2];
  (*cmat)(0,3)=ET[1][0];
  (*cmat)(0,4)=ET[2][1];
  (*cmat)(0,5)=ET[2][0];
  
  (*cmat)(1,0)=ET[3][3];
  (*cmat)(1,1)=ET[4][4];
  (*cmat)(1,2)=ET[5][5];
  (*cmat)(1,3)=ET[4][3];
  (*cmat)(1,4)=ET[5][4];
  (*cmat)(1,5)=ET[5][3];
  
  (*cmat)(2,0)=ET[6][6];
  (*cmat)(2,1)=ET[7][7];
  (*cmat)(2,2)=ET[8][8];
  (*cmat)(2,3)=ET[7][6];
  (*cmat)(2,4)=ET[8][7];
  (*cmat)(2,5)=ET[8][6];
  
  (*cmat)(3,0)=ET[3][0];
  (*cmat)(3,1)=ET[4][1];
  (*cmat)(3,2)=ET[5][2];
  (*cmat)(3,3)=ET[4][0];
  (*cmat)(3,4)=ET[5][1];
  (*cmat)(3,5)=ET[5][0];
  
  (*cmat)(4,0)=ET[6][3];
  (*cmat)(4,1)=ET[7][4];
  (*cmat)(4,2)=ET[8][5];
  (*cmat)(4,3)=ET[7][3];
  (*cmat)(4,4)=ET[8][4];
  (*cmat)(4,5)=ET[8][3];
  
  (*cmat)(5,0)=ET[6][0];
  (*cmat)(5,1)=ET[7][1];
  (*cmat)(5,2)=ET[8][2];
  (*cmat)(5,3)=ET[7][0];
  (*cmat)(5,4)=ET[8][1];
  (*cmat)(5,5)=ET[8][0];
 
  return;
} // end of neohooke evaluate


/*----------------------------------------------------------------------*
 |  Calculate the inverse of a 2nd order tensor                 rm 08/07|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseMatrix MAT::NeoHooke::inverseTensor(
						const Epetra_SerialDenseMatrix M,
						const double I3)
{
  Epetra_SerialDenseMatrix Minv(3,3);
  if (I3==0) {
  	printf("Matrix nicht invertierbar!\n");
  } else {
  	Minv(0,0)= 1/I3 * (M(1,1)*M(2,2) - M(2,1)*M(1,2));
	Minv(1,0)=-1/I3 * (M(0,1)*M(2,2) - M(2,1)*M(0,2));
	Minv(2,0)= 1/I3 * (M(0,1)*M(1,2) - M(1,1)*M(0,2));
	Minv(0,1)=-1/I3 * (M(1,0)*M(2,2) - M(2,0)*M(1,2));
	Minv(1,1)= 1/I3 * (M(0,0)*M(2,2) - M(2,0)*M(0,2));
	Minv(2,1)=-1/I3 * (M(0,0)*M(1,2) - M(1,0)*M(0,2));
	Minv(0,2)= 1/I3 * (M(1,0)*M(2,1) - M(2,0)*M(1,1));
	Minv(1,2)=-1/I3 * (M(0,0)*M(2,1) - M(2,0)*M(0,1));
	Minv(2,2)= 1/I3 * (M(0,0)*M(1,1) - M(1,0)*M(0,1));
   }
return Minv;
}


/*----------------------------------------------------------------------*
 |  Returns the density of the material                         rm 08/07|
 *----------------------------------------------------------------------*/
double MAT::NeoHooke::Density()
{
  return matdata_->m.neohooke->density;  // density, returned to evaluate mass matrix
}

#endif
