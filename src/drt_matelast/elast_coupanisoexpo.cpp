/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpo.H
\brief


the input line should read
  MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0  K1COMP 0.0 K2COMP 1.0 [INIT 1] [ADAPT_ANGLE No]

<pre>
Maintainer: Susanna Tinkl
            tinkl@lnm.mw.tum.de
            089/289 15265
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpo.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpo::CoupAnisoExpo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  gamma_(matdata->GetDouble("GAMMA")),
  k1comp_(matdata->GetDouble("K1COMP")),
  k2comp_(matdata->GetDouble("K2COMP")),
  init_(matdata->GetInt("INIT")),
  adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   st         03/12 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoExpo::CoupAnisoExpo(MAT::ELASTIC::PAR::CoupAnisoExpo* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,a_);
  AddtoPack(data,A_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::UnpackSummand(const std::vector<char>& data,
                                                  std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // path if fibers aren't given in .dat file
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3,3> Id(true);
    for (int i=0; i<3; i++)
      Id(i,i) = 1.0;
    SetFiberVecs(-1.0,Id,Id);
  }

  // path if fibers are given in .dat file
  else if (params_->init_ == 1)
  {

    // CIR-AXI-RAD nomenclature
    if (linedef->HaveNamed("RAD") and
        linedef->HaveNamed("AXI") and
        linedef->HaveNamed("CIR"))
    {
      // Read in of data
      LINALG::Matrix<3,3> locsys(true);
      ReadRadAxiCir(linedef, locsys);
      LINALG::Matrix<3,3> Id(true);
      for (int i=0; i<3; i++)
        Id(i,i) = 1.0;
      // final setup of fiber data
      SetFiberVecs(0.0,locsys,Id);
    }

    // FIBER1 nomenclature
    else if ( linedef->HaveNamed("FIBER1") )
    {
      // Read in of fiber data and setting fiber data
      ReadFiber(linedef, "FIBER1", a_);
      SetupStructuralTensor(a_,A_);
    }

    // error path
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
void MAT::ELASTIC::CoupAnisoExpo::AddStressAnisoPrincipal(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress,
    Teuchos::ParameterList& params
)
{
  double I4 = 0.0;
  I4 =  A_(0)*rcg(0) + A_(1)*rcg(1) + A_(2)*rcg(2)
      + A_(3)*rcg(3) + A_(4)*rcg(4) + A_(5)*rcg(5);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  double gamma = 2.*(k1*(I4-1.)* exp(k2*(I4-1.)*(I4-1.)));
  stress.Update(gamma, A_, 1.0);

  double delta = 2.*(1. + 2.*k2*(I4-1.)*(I4-1.))*2.*k1* exp(k2*(I4-1.)*(I4-1.));
  cmat.MultiplyNT(delta, A_, A_, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::GetFiberVecs(
    std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::SetFiberVecs(
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

  SetupStructuralTensor(a_,A_);
}
