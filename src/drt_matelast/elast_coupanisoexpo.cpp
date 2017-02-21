/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpo.cpp
\brief the input line should read
  MAT 1 ELAST_CoupAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0  K1COMP 0.0 K2COMP 1.0 [INIT 1] [ADAPT_ANGLE No]

\level 1

<pre>
\maintainer Fabian Braeu
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
void MAT::ELASTIC::CoupAnisoExpo::AddStrainEnergy(double& psi,
                                                  const LINALG::Matrix<3,1>& prinv,
                                                  const LINALG::Matrix<3,1>& modinv,
                                                  const LINALG::Matrix<6,1>& glstrain,
                                                  const int eleGID)
{
  // rigth Cauchy Green in strain-like Voigt notation
  LINALG::Matrix<6,1> rcg(true);

  for(int i=0;i<3;++i)
    rcg(i) = 2.0*glstrain(i)+1.0;
  rcg(3) = 2.0*glstrain(3);
  rcg(4) = 2.0*glstrain(4);
  rcg(5) = 2.0*glstrain(5);

  double I4 = 0.0;
  I4 =  A_(0)*rcg(0) + A_(1)*rcg(1) + A_(2)*rcg(2)
      + A_(3)*rcg(3) + A_(4)*rcg(4) + A_(5)*rcg(5);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (prinv(0) < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  psi += (k1/(2.0*k2))*(exp(k2*(I4-1.0)*(I4-1.0)) - 1.0);

  return;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename T >
void MAT::ELASTIC::CoupAnisoExpo::EvaluateFunc(T& psi,
                                               LINALG::TMatrix<T,3,3> const& rcg,
                                               int const eleGID) const
{
  T I4_fad = 0.0;
  static LINALG::TMatrix<T,6,1> Av_T(true);
  for(int i=0;i<6;++i) Av_T(i) = A_(i);
  I4_fad =  Av_T(0)*rcg(0,0) + Av_T(1)*rcg(1,1) + Av_T(2)*rcg(2,2)
      + Av_T(3)*(rcg(0,1)+rcg(1,0)) + Av_T(4)*(rcg(1,2)+rcg(2,1)) + Av_T(5)*(rcg(0,2)+rcg(2,0));

  T k1 = params_->k1_;
  T k2 = params_->k2_;

  if (I4_fad < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  psi += (k1/(2.0*k2))*(exp(k2*(I4_fad-1.0)*(I4_fad-1.0)) - 1.0);

  return;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template< typename T >
void MAT::ELASTIC::CoupAnisoExpo::GetDerivativesAniso(LINALG::TMatrix<T,2,1>& dPI_aniso,
                                                      LINALG::TMatrix<T,3,1>& ddPII_aniso,
                                                      LINALG::TMatrix<T,4,1>& dddPIII_aniso,
                                                      LINALG::TMatrix<T,3,3> const& rcg,
                                                      const int eleGID) const
{
  T I4 = 0.0;
  LINALG::TMatrix<T,3,3> AM(true);
  for(int i=0;i<3;++i) AM(i,i) = A_(i);
  AM(0,1) = AM(1,0) = A_(3);
  AM(2,1) = AM(1,2) = A_(4);
  AM(0,2) = AM(2,0) = A_(5);

  I4 =  AM.Dot(rcg);

  T k1 = params_->k1_;
  T k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }


  dPI_aniso(0) = k1*(I4-1.0)*exp(k2*(I4-1.0)*(I4-1.0));

  ddPII_aniso(0) = (1.0+2.0*k2*(I4-1.0)*(I4-1.0))*k1*exp(k2*(I4-1.0)*(I4-1.0));

  dddPIII_aniso(0) = (3.0 + 2.0*k2*(I4-1.0)*(I4-1.0))*2.0*k1*k2*(I4-1.0)*exp(k2*(I4-1.0)*(I4-1.0));

  return;
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpo::AddStressAnisoPrincipal(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress,
    Teuchos::ParameterList& params,
    const int eleGID
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


// explicit instantiation of template functions
template void MAT::ELASTIC::CoupAnisoExpo::GetDerivativesAniso<double>(LINALG::TMatrix<double,2,1>&,
                                                                       LINALG::TMatrix<double,3,1>&,
                                                                       LINALG::TMatrix<double,4,1>&,
                                                                       LINALG::TMatrix<double,3,3> const&,
                                                                       const int) const;
template void MAT::ELASTIC::CoupAnisoExpo::GetDerivativesAniso<FAD>(LINALG::TMatrix<FAD,2,1>&,
                                                                       LINALG::TMatrix<FAD,3,1>&,
                                                                       LINALG::TMatrix<FAD,4,1>&,
                                                                       LINALG::TMatrix<FAD,3,3> const&,
                                                                       const int) const;
template void MAT::ELASTIC::CoupAnisoExpo::EvaluateFunc<double>(double&,
                                                                LINALG::TMatrix<double,3,3> const&,
                                                                const int) const;
template void MAT::ELASTIC::CoupAnisoExpo::EvaluateFunc<FAD>(FAD&,
                                                             LINALG::TMatrix<FAD,3,3> const&,
                                                             const int) const;
