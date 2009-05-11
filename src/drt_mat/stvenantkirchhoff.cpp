/*----------------------------------------------------------------------*/
/*!
\file stvenantkirchhoff.cpp

<pre>
Maintainer: ???
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "stvenantkirchhoff.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::StVenantKirchhoff::StVenantKirchhoff(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  thermexpans_(matdata->GetDouble("THEXPANS"))
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::StVenantKirchhoff::StVenantKirchhoff()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::StVenantKirchhoff::StVenantKirchhoff(MAT::PAR::StVenantKirchhoff* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::StVenantKirchhoff*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}

/*----------------------------------------------------------------------*
// computes isotropic eplane strain, rotational symmetry
// plane strain, rotational symmetry
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::SetupCmat2d(Epetra_SerialDenseMatrix* cmat)
{
  const double ym  = params_->youngs_;
  const double pv  = params_->poissonratio_;

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

/*----------------------------------------------------------------------*
// computes isotropic elasticity tensor in matrix notion for 3d
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::SetupCmat(LINALG::Matrix<6,6>& cmat)
{
  // get material parameters
  const double Emod = params_->youngs_;    // Young's modulus (modulus of elasticity)
  const double nu = params_->poissonratio_;// Poisson's ratio (Querdehnzahl)

/*
  if (nu == 0.5) {
    // linearly isochoric. i.e. deviatoric, isotropic elasticity tensor C in Voigt matrix notation
    //                       [  2/3   -1/3   -1/3 |   0    0    0 ]
    //                       [         2/3   -1/3 |   0    0    0 ]
    //           E           [                2/3 |   0    0    0 ]
    //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~  ~~~  ~~~ ]
    //       (1+nu)          [                    | 1/2    0    0 ]
    //                       [                    |      1/2    0 ]
    //                       [ symmetric          |           1/2 ]
    //
    const double mfac = Emod/(1.0+nu);  // 2x shear modulus
    cmat(0,0) = mfac*2.0/3.0;
    cmat(0,1) = -mfac*1.0/3.0;
    cmat(0,2) = -mfac*1.0/3.0;
    cmat(1,0) = -mfac*1.0/3.0;
    cmat(1,1) = mfac*2.0/3.0;
    cmat(1,2) = -mfac*1.0/3.0;
    cmat(2,0) = -mfac*1.0/3.0;
    cmat(2,1) = -mfac*1.0/3.0;
    cmat(2,2) = mfac*2.0/3.0;
    // ~~~
    cmat(3,3) = mfac*0.5;
    cmat(4,4) = mfac*0.5;
    cmat(5,5) = mfac*0.5;
  }
  else */ {
    // isotropic elasticity tensor C in Voigt matrix notation
    //                       [ 1-nu     nu     nu |          0    0    0 ]
    //                       [        1-nu     nu |          0    0    0 ]
    //           E           [               1-nu |          0    0    0 ]
    //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
    //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
    //                       [                    |      (1-2*nu)/2    0 ]
    //                       [ symmetric          |           (1-2*nu)/2 ]
    //
    const double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  // factor
    // write non-zero components
    cmat(0,0) = mfac*(1.0-nu);
    cmat(0,1) = mfac*nu;
    cmat(0,2) = mfac*nu;
    cmat(1,0) = mfac*nu;
    cmat(1,1) = mfac*(1.0-nu);
    cmat(1,2) = mfac*nu;
    cmat(2,0) = mfac*nu;
    cmat(2,1) = mfac*nu;
    cmat(2,2) = mfac*(1.0-nu);
    // ~~~
    cmat(3,3) = mfac*0.5*(1.0-2.0*nu);
    cmat(4,4) = mfac*0.5*(1.0-2.0*nu);
    cmat(5,5) = mfac*0.5*(1.0-2.0*nu);
  }
}


/*----------------------------------------------------------------------*
//calculates stresses using one of the above method to evaluate the elasticity tensor
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Evaluate(const Epetra_SerialDenseVector* glstrain_e,
                                      Epetra_SerialDenseMatrix* cmat_e,
                                      Epetra_SerialDenseVector* stress_e)
{
  // this is temporary as long as the material does not have a 
  // Matrix-type interface
  const LINALG::Matrix<6,1> glstrain(glstrain_e->A(),true);
        LINALG::Matrix<6,6> cmat(cmat_e->A(),true);
        LINALG::Matrix<6,1> stress(stress_e->A(),true);

  SetupCmat(cmat);
  // evaluate stresses
  stress.MultiplyNN(cmat,glstrain);  // sigma = C . epsilon
}


/*----------------------------------------------------------------------*
//calculates stresses using one of the above method to evaluate the elasticity tensor
 *----------------------------------------------------------------------*/
void MAT::StVenantKirchhoff::Evaluate(
                  const LINALG::Matrix<6,1>& glstrain,
                  LINALG::Matrix<6,6>& cmat,
                  LINALG::Matrix<6,1>& stress)
{
  SetupCmat(cmat);
  // evaluate stresses
  stress.MultiplyNN(cmat,glstrain);  // sigma = C . epsilon
}


#endif
