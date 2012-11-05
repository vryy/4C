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
/* Local Changes:
11/2012 A. Nagler:
	- Added optional angle theta as input in order to use spherical coordinates as fiber initialization
	- Reorganization of Setup() and SetFiberVecs()
	- Deleted AAA initialization method via CIR, AXI, RAD
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
  init_(matdata->GetInt("INIT")),
  gamma_(matdata->GetDouble("GAMMA")),
  theta_(matdata->GetDouble("THETA"))
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
  LINALG::Matrix<3,3> Id(true);
  LINALG::Matrix<3,3> locsys(true);

  for (int i=0; i<3; i++)
    Id(i,i) = 1.0;

	// generate local fiber coordinate system from given angles
  if (params_->init_ == 0)
  {
    // setting up analytical fiber direction via the spherical coordinates gamma (azimuth angle) and theta(polar angle)
    double gamma = (params_->gamma_ * PI)/180.0;
    double theta = (params_->theta_ * PI)/180.0;

    locsys(0,0) = sin(theta)*cos(gamma);
    locsys(1,0) = sin(theta)*sin(gamma);
    locsys(2,0) = cos(theta);

    SetFiberVecs(0.0, locsys, Id);
  }
  else if (params_->init_ == 1)
  {
    if (linedef->HaveNamed("FIBER1"))
    {
    	// Reading of fiber directions
    	std::vector<double> fiber1;
      linedef->ExtractDoubleVector("FIBER1",fiber1);

      // Normalization of fiber vectors to 1
      double norm_f1=0.;
      for (int i = 0; i < 3; ++i)
      {
      	norm_f1 += fiber1[i]*fiber1[i];
      }
      norm_f1 = sqrt(norm_f1);

      // setting up of local coordinate system
      for (int i=0; i<3; ++i)
      {
        locsys(i,0) = fiber1[i]/norm_f1;
      }

      SetFiberVecs(0.0, locsys,Id);
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
   double stressFact_=params.get<double>("scalar", 1.0);
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

	// The first row of locsys is equal to the fiber
  LINALG::Matrix<3,1> a_0 (true);
  for (int i=0; i<3; i++)
  {
  	a_0(i) = locsys(i,0);
  }

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

