/*----------------------------------------------------------------------*/
/*!
\file thermostvenantkirchhoff.cpp
\brief St.Venant Kirchhoff with an additional temperature dependent term
       describing heat expansion

       example input line:
       MAT 1   MAT_Struct_ThrStVenantK   YOUNG 1.48e8 NUE 0.3 DENS 9.130e-6
       THEXPANS 1.72e-5 INITTEMP 293.15

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 02/10 |
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 | headers                                                   dano 02/10 |
 *----------------------------------------------------------------------*/
#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "thermostvenantkirchhoff.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*
 |                                                           dano 02/10 |
 *----------------------------------------------------------------------*/
MAT::PAR::ThermoStVenantKirchhoff::ThermoStVenantKirchhoff(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  thermexpans_(matdata->GetDouble("THEXPANS")),
  capa_(matdata->GetDouble("CAPA")),
  conduct_(matdata->GetDouble("CONDUCT")),
  thetainit_(matdata->GetDouble("INITTEMP"))
{
}

/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ThermoStVenantKirchhoff::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ThermoStVenantKirchhoff(this));
}


MAT::ThermoStVenantKirchhoffType MAT::ThermoStVenantKirchhoffType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::ThermoStVenantKirchhoffType::Create(
  const std::vector<char> & data
  )
{
  MAT::ThermoStVenantKirchhoff* thrstvenantk = new MAT::ThermoStVenantKirchhoff();
  thrstvenantk->Unpack(data);
  return thrstvenantk;
}

/*----------------------------------------------------------------------*
 |  constructor (public)                                     dano 02/10 |
 *----------------------------------------------------------------------*/
MAT::ThermoStVenantKirchhoff::ThermoStVenantKirchhoff()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 02/10 |
 *----------------------------------------------------------------------*/
MAT::ThermoStVenantKirchhoff::ThermoStVenantKirchhoff(
  MAT::PAR::ThermoStVenantKirchhoff* params
  )
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |  Pack (public)                                            dano 02/10 |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*
 |  Unpack (public)                                          dano 02/10 |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ThermoStVenantKirchhoff*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 | computes isotropic eplane strain, rotational symmetry     dano 02/10 |
 | plane strain, rotational symmetry                                    |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::SetupCmat2d(Epetra_SerialDenseMatrix* cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  const double ym = params_->youngs_;
  // Poisson's ratio (Querdehnzahl)
  const double pv = params_->poissonratio_;

  // plane strain, rotational symmetry
  // \f \frac{E}{1+\nu} \f
  const double c1 = ym/(1.0+pv);
  // \f  = \frac{E\,\nu}{(1-2\nu)(1+\nu)} \f
  const double b1 = c1*pv/(1.0-2.0*pv);
  // for the diagonal terms a1 = nu+lambda
  const double a1 = b1+c1;

  (*cmat)(0,0) = a1;
  (*cmat)(0,1) = b1;
  (*cmat)(0,2) = 0.;
  (*cmat)(0,3) = b1;

  (*cmat)(1,0) = b1;
  (*cmat)(1,1) = a1;
  (*cmat)(1,2) = 0.;
  (*cmat)(1,3) = b1;

  (*cmat)(2,0) = 0.;
  (*cmat)(2,1) = 0.;
  (*cmat)(2,2) = c1/2.;
  (*cmat)(2,3) = 0.;

  (*cmat)(3,0) = b1;
  (*cmat)(3,1) = b1;
  (*cmat)(3,2) = 0.;
  (*cmat)(3,3) = a1;
}


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 02/10 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::SetupCmat(LINALG::Matrix<6,6>& cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  const double Emod = params_->youngs_;
  // Poisson's ratio (Querdehnzahl)
  const double nu = params_->poissonratio_;

/*
  if (nu == 0.5) {
    // linearly isochoric. i.e. deviatoric, isotropic elasticity tensor C in
    // Voigt matrix notation
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
  else */
  {
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

    // clear the material tangent
    cmat.Clear();
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
 | calculates stresses using one of the above method to      dano 02/10 |
 | evaluate the elasticity tensor                                       |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::Evaluate(
  const Epetra_SerialDenseVector* glstrain_e,
  Epetra_SerialDenseMatrix* cmat_e,
  Epetra_SerialDenseVector* stress_e
  )
{
  // this is temporary as long as the material does not have a
  // Matrix-type interface
  const LINALG::Matrix<6,1> glstrain(glstrain_e->A(),true);
        LINALG::Matrix<6,6> cmat(cmat_e->A(),true);
        LINALG::Matrix<6,1> stress(stress_e->A(),true);

  SetupCmat(cmat);
  // evaluate stresses
  // \f \sigma = {\mathbf C}\, \varepsilon \f
  stress.MultiplyNN(cmat,glstrain);
}


/*----------------------------------------------------------------------*
 | calculates stresses using one of the above method to      dano 02/10 |
 | evaluate the elasticity tensor                                       |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress
  )
{
  SetupCmat(cmat);
  // evaluate stresses
  // \f \sigma = {\mathbf C} \,\varepsilon \f
  stress.MultiplyNN(cmat,glstrain);
}


/*----------------------------------------------------------------------*
 | computes temperature dependent isotropic eplane           dano 02/10 |
 | strain, rotational symmetry,plane strain, rotational symmetry        |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::SetupCthermo2d(
  Epetra_SerialDenseMatrix* ctemp
  )
{
  double m = STModulus();

  // add the temperature part for the stress 10.02.10
  (*ctemp)(0,0) = m;
  (*ctemp)(0,1) = 0.;
  (*ctemp)(1,0) = 0.;
  (*ctemp)(1,1) = m;
}


/*----------------------------------------------------------------------*
 | calculates stress-temperature modulus                     dano 04/10 |
 *----------------------------------------------------------------------*/
double MAT::ThermoStVenantKirchhoff::STModulus() const
{
  // initialize the parameters for the lame constants
  const double ym  = params_->youngs_;
  const double pv  = params_->poissonratio_;

  // initialize the thermal expansion coefficient
  const double thermexpans = params_->thermexpans_;

  // plane strain, rotational symmetry
  // E / (1+nu)
  const double c1 = ym/(1.0+pv);
  // (E*nu) / ((1+nu)(1-2nu))
  const double b1 = c1*pv/(1.0-2.0*pv);

  // build the lame constants
  //            E
  //   mu = --------
  //        2*(1+nu)
  //                  E*nu
  //   lambda = ----------------
  //            (1+nu)*(1-2*nu)
  //
  //  \f \mu =  \frac{E}{2(1+\nu)} \f
  const double mu = 0.5*c1;
  // lambda
  // \f \frac{E\,\nu}{(1-2\nu)(1+\nu)} \f
  const double lambda = b1;

  // stress-temperature modulus
  // \f m\, = \, -(2\,\cdot \mu \, +\, 3\cdot\lambda)\cdot\varalpha_T \f
  const double stmodulus = (-1)*(2*mu + 3*lambda)*thermexpans;

  return stmodulus;

} // STModulus()


/*----------------------------------------------------------------------*
 | computes temperature dependent isotropic                  dano 05/10 |
 | elasticity tensor in matrix notion for 3d, second(!) order tensor    |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::SetupCthermo(LINALG::Matrix<6,1>& ctemp)
{
  double m = STModulus();

  // isotropic elasticity tensor C_temp in Voigt matrix notation C_temp = m I
  //
  // Matrix-notation for 3D case
  //              [ m      0      0 ]
  //   C_temp =   [ 0      m      0 ]
  //              [ 0      0      m ]
  //
  //  in Vector notation
  //   C_temp =   [m, m, m, 0, 0, 0]^T
  //
  // write non-zero components

  // clear the material tangent
  ctemp.Clear();

  // loop over the element nodes
  for (int i=0; i<3; ++i)
    // non-zero entries only in main directions
    ctemp(i,0) = m;
  for (int i=3; i<6; ++i)
    ctemp(i,0) = 0;
}


/*----------------------------------------------------------------------*
 | calculates stresses by evaluating the temperature tangent dano 05/10 |
 *----------------------------------------------------------------------*/
void MAT::ThermoStVenantKirchhoff::Evaluate(
  const LINALG::Matrix<1,1>& Ntemp,  // shapefcts . temperatures
  LINALG::Matrix<6,1>& ctemp,
  LINALG::Matrix<6,1>& stresstemp
  )
{
  SetupCthermo(ctemp);

  LINALG::Matrix<1,1> init(true);
  const double inittemp = -1.0*(params_->thetainit_);
  // loop over the element nodes
  init(0,0) = inittemp;
  // Delta T = T - T_0
  LINALG::Matrix<1,1> deltaT(true);
  deltaT.Update(Ntemp,init);

  // temperature dependent stress
  // sigma = C_theta * Delta T = (m*I) * Delta T
  stresstemp.MultiplyNN(ctemp,deltaT);

} // Evaluate


/*----------------------------------------------------------------------*/


