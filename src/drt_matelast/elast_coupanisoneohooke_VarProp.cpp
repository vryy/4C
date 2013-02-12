/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooke_VarProp.cpp
\brief


the input line should read
  MAT 1 CoupAnisoNeoHooke_VarProp C 100 GAMMA 35.0 INIT 0 ADAPT_ANGLE 0

<pre>
Maintainer: Cristobal Bertoglio
            bertoglio@lnm.mw.tum.de
            089/289 15264
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoneohooke_VarProp.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp::CoupAnisoNeoHooke_VarProp(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  gamma_(matdata->GetDouble("GAMMA")),
  theta_(matdata->GetDouble("THETA")),
  init_(matdata->GetInt("INIT")),
  adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::CoupAnisoNeoHooke( this ) );
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::CoupAnisoNeoHooke_VarProp()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::CoupAnisoNeoHooke_VarProp(MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,a_);
  AddtoPack(data,A_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::UnpackSummand(const std::vector<char>& data,
																														std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::Setup(DRT::INPUT::LineDefinition* linedef)
{
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3,3> Id(true);
    LINALG::Matrix<3,3> locsys(true);
    for (int i=0; i<3; i++)
      Id(i,i) = 1.0;

    // To realize a full rotated fiber orientation and to keep the general structure of SetFiberVec() the
    // input of locsys has to be adapted
    // if one sets
    //               1           0                 sin(theta_)
    //  locsys := [  0       sin(theta_)               0                  ]
    //               0  sin(gamma)*cos(theta_)    cos(gamma_)*cos(theta_)
    // The call of SetFiberVec() will leed to the following fiber direction
    // a = cos(gamma_)*locsys(:,2) + sin(gamma_)*locsys(:,1)
    //         cos(gamma_)*sin(theta_)               0                          cos(gamma_)*sin(theta_)
    //   = [             0              ] + [  sin(gamma_)*sin(theta_)  ] = [   sin(gamma_)*sin(theta_) ] =: sperical coordinates
    //         cos(gamma)^2*cos(theta_)        sin(gamma_)^2*cos(theta_)               cos(theta_)
    //
    {
    	// Local initialization of spherical angles
      double theta = (params_->theta_);
      double gamma = (params_->gamma_);
      if ( gamma < 0.0 || gamma > 180.0 || abs(theta) > 180.0 )
	    {
      		dserror("Wrong choice of sherical coodinates. Correct domain is gamma in [0,180], theta in [-180, 180]");
	    }
      // conversion to radian measure
      theta = (theta*PI)/180.0;
      gamma = (gamma*PI)/180.0;
      locsys(1,1) = sin(theta);
      locsys(2,1) = sin(gamma)*cos(theta);
      locsys(0,2) = sin(theta);
      locsys(2,2) = cos(gamma)*cos(theta);
    }
    SetFiberVecs(-1.0,locsys,Id);
  }
  else if (params_->init_ == 1)
  {
    // fibers aligned in local element cosy with gamma around circumferential direction
    // -> check whether element supports local element cosy
    if (linedef->HaveNamed("RAD") and
        linedef->HaveNamed("AXI") and
        linedef->HaveNamed("CIR"))
    {
      // read local (cylindrical) cosy-directions at current element
      // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
    	std::vector<double> rad;
    	std::vector<double> axi;
    	std::vector<double> cir;
      LINALG::Matrix<3,3> locsys(true);
      linedef->ExtractDoubleVector("RAD",rad);
      linedef->ExtractDoubleVector("AXI",axi);
      linedef->ExtractDoubleVector("CIR",cir);
      double radnorm=0.; double axinorm=0.; double cirnorm=0.;

      for (int i = 0; i < 3; ++i)
      {
        radnorm += rad[i]*rad[i]; axinorm += axi[i]*axi[i]; cirnorm += cir[i]*cir[i];
      }
      radnorm = sqrt(radnorm); axinorm = sqrt(axinorm); cirnorm = sqrt(cirnorm);

      for (int i=0; i<3; ++i)
      {
        locsys(i,0) = rad[i]/radnorm;
        locsys(i,1) = axi[i]/axinorm;
        locsys(i,2) = cir[i]/cirnorm;
      }

      LINALG::Matrix<3,3> Id(true);
      for (int i=0; i<3; i++)
        Id(i,i) = 1.0;
      SetFiberVecs(0.0,locsys,Id);
    }
    // read given first fiber family
    else if ( linedef->HaveNamed("FIBER1") )
    {
    	std::vector<double> fiber1;
      LINALG::Matrix<3,3> locsys(true);
      linedef->ExtractDoubleVector("FIBER1",fiber1);
      double f1norm=0.;

      //normalization
      for (int i = 0; i < 3; ++i)
      {
      	f1norm += fiber1[i]*fiber1[i];
      }
      f1norm = sqrt(f1norm);

      // we set locsys(:,2) = fiber1, since in function SetFiberVecs
      // the fiber orientation will be calculated via
      // ca = cos(gamma)*locsys(:,2) + sin(gamma)*locsys(:,1)
      //    = locsys(:,2) for gamma=0.0
      for (int i=0; i<3; ++i)
      {
        locsys(i,2) = fiber1[i]/f1norm;
      }

      LINALG::Matrix<3,3> Id(true);
      for (int i=0; i<3; i++)
        Id(i,i) = 1.0;
      SetFiberVecs(0.0,locsys,Id);
    }
    else
    {
      dserror("Reading of element local cosy for anisotropic materials failed");
    }
  }
  else
    dserror("INIT mode not implemented");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::AddStressAnisoPrincipal(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress,
    Teuchos::ParameterList& params
)
{
   double stressFact_=params.get<double>("scalar");
   stress.Update(2*(params_->c_)*stressFact_, A_, 1.0);

   // no contribution to cmat
  // double delta = 0.0;
  // cmat.MultiplyNT(delta, A_, A_, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::SetFiberVecs(
    const double newgamma,
    const LINALG::Matrix<3,3> locsys,
    const LINALG::Matrix<3,3> defgrd
)
{
  if ((params_->gamma_<-90) || (params_->gamma_ >90)) dserror("Fiber angle not in [-90,90]");
  //convert
  double gamma = (params_->gamma_*PI)/180.;

  if (params_->adapt_angle_ && newgamma != -1.0)
  {
    if (gamma*newgamma < 0.0)
      gamma = -1.0 * newgamma;
    else
      gamma = newgamma;
  }

  LINALG::Matrix<3,1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = cos(gamma)*locsys(i,2) + sin(gamma)*locsys(i,1);
  }
  // pull back in reference configuration
  LINALG::Matrix<3,1> a_0(true);
  LINALG::Matrix<3,3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a_0.Multiply(idefgrd,ca);
  a_.Update(1./a_0.Norm2(),a_0);

  for (int i = 0; i < 3; ++i)
    A_(i) = a_(i)*a_(i);

  A_(3) = a_(0)*a_(1); A_(4) = a_(1)*a_(2); A_(5) = a_(0)*a_(2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::GetFiberVecs(
    std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

