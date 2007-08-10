
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "stvenantkirchhoff.H"

extern struct _MATERIAL *mat;

MAT::StVenantKirchhoff::StVenantKirchhoff()
  : matdata_(NULL)
{
}


MAT::StVenantKirchhoff::StVenantKirchhoff(MATERIAL* matdata)
  : matdata_(matdata)
{
}

void MAT::StVenantKirchhoff::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
}


void MAT::StVenantKirchhoff::Unpack(const vector<char>& data)
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

// computes isotropic eplane strain, rotational symmetry
// plane strain, rotational symmetry
void MAT::StVenantKirchhoff::SetupCmat2d(Epetra_SerialDenseMatrix* cmat)
{
  const double ym  = matdata_->m.stvenant->youngs;
  const double pv  = matdata_->m.stvenant->possionratio;

  // plane strain, rotational symmetry
  const double c1=ym/(1.0+pv);
  const double b1=c1*pv/(1.0-2.0*pv);
  const double a1=b1+c1;

  (*cmat)(0,0)=a1;
  (*cmat)(0,1)=b1;
  (*cmat)(0,2)=0.;
  (*cmat)(0,3)=b1;
  
  (*cmat)(1,0)=b1;
  (*cmat)(1,1)=a1;
  (*cmat)(1,2)=0.;
  (*cmat)(1,3)=b1;
  
  (*cmat)(2,0)=0.;
  (*cmat)(2,1)=0.;
  (*cmat)(2,2)=c1/2.;
  (*cmat)(2,3)=0.;
  
  (*cmat)(3,0)=b1;
  (*cmat)(3,1)=b1;
  (*cmat)(3,2)=0.;
  (*cmat)(3,3)=a1;
}

// computes isotropic elasticity tensor in matrix notion for 3d
void MAT::StVenantKirchhoff::SetupCmat(Epetra_SerialDenseMatrix* cmat)
{
  // get material parameters
  double Emod = matdata_->m.stvenant->youngs;    // Young's modulus (modulus of elasticity)
  double nu = matdata_->m.stvenant->possionratio;// Poisson's ratio (Querdehnzahl)

  /*--------------------------------------------------------------------*/
  /* isotropic elasticity tensor C in matrix notion */
  /*                       [ 1-nu     nu     nu |          0    0    0 ]
   *                       [        1-nu     nu |          0    0    0 ]
   *           E           [               1-nu |          0    0    0 ]
   *   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
   *       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
   *                       [                    |      (1-2*nu)/2    0 ]
   *                       [ symmetric          |           (1-2*nu)/2 ]
   */
  double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
  /* write non-zero components */
  (*cmat)(0,0) = mfac*(1.0-nu);
  (*cmat)(0,1) = mfac*nu;
  (*cmat)(0,2) = mfac*nu;
  (*cmat)(1,0) = mfac*nu;
  (*cmat)(1,1) = mfac*(1.0-nu);
  (*cmat)(1,2) = mfac*nu;
  (*cmat)(2,0) = mfac*nu;
  (*cmat)(2,1) = mfac*nu;
  (*cmat)(2,2) = mfac*(1.0-nu);
  /* ~~~ */
  (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);
}


//calculates stresses using one of the above method to evaluate the elasticity tensor
void MAT::StVenantKirchhoff::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                      Epetra_SerialDenseMatrix* cmat,
                                      Epetra_SerialDenseVector* stress)
{
  SetupCmat(cmat);
  // evaluate stresses
  (*cmat).Multiply('N',(*glstrain),(*stress));   // sigma = C . epsilon
}


double MAT::StVenantKirchhoff::Density()
{
  return matdata_->m.stvenant->density;  // density, returned to evaluate mass matrix
}


#endif
