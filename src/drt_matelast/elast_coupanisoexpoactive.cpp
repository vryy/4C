/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpoactive.cpp
\brief the input line should read
     MAT 1 ELAST_CoupAnisoExpoActive K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0 S 54000 LAMBDAMAX 1.4 LAMBDA0 0.8 DENS 1050

\level 3

<pre>
\maintainer Fabian Braeu
            braeu@lnm.mw.tum.de
            089/289 15236
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpoactive.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpoActive::CoupAnisoExpoActive(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  gamma_(matdata->GetDouble("GAMMA")),
  k1comp_(matdata->GetDouble("K1COMP")),
  k2comp_(matdata->GetDouble("K2COMP")),
  init_(matdata->GetInt("INIT")),
  adapt_angle_(matdata->GetInt("ADAPT_ANGLE")),
  s_(matdata->GetDouble("S")),
  lambdamax_(matdata->GetDouble("LAMBDAMAX")),
  lambda0_(matdata->GetDouble("LAMBDA0")),
  dens_(matdata->GetDouble("DENS"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   fb         07/16 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoExpoActive::CoupAnisoExpoActive(MAT::ELASTIC::PAR::CoupAnisoExpoActive* params)
  : params_(params)
{
  dPIact_ = 0.0;
  lambdaact_ = 1.0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,a_);
  AddtoPack(data,A_);
  AddtoPack(data,dPIact_);
  AddtoPack(data,lambdaact_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::UnpackSummand(const std::vector<char>& data,
                                                      std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
  ExtractfromPack(position,data,dPIact_);
  ExtractfromPack(position,data,lambdaact_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::Setup(DRT::INPUT::LineDefinition* linedef)
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

  // setup first derivative of active fiber potential w.r.t. active fiber stretch (const during the whole simulation)
  lambdaact_ = 1.0;

  dPIact_ = params_->s_/params_->dens_*(1.0 - std::pow(params_->lambdamax_ - lambdaact_,2.0)/std::pow(params_->lambdamax_ - params_->lambda0_,2.0));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::AddStrainEnergy(double& psi,
                                                        const LINALG::Matrix<3,1>& prinv,
                                                        const LINALG::Matrix<3,1>& modinv,
                                                        const LINALG::Matrix<6,1> glstrain,
                                                        const int eleGID)
{
  // rigth Cauchy Green
  LINALG::Matrix<6,1> rcg(true);

  for(int i=0;i<3;++i)
    rcg(i) = 2.0*glstrain(i,i)+1.0;
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

  // passive contribution
  psi += (k1/(2.0*k2))*(exp(k2*(I4-1.0)*(I4-1.0)) - 1.0);

  // active contribution
  psi += params_->s_/params_->dens_*(lambdaact_ + (1./3.)*(std::pow(params_->lambdamax_ - lambdaact_,3.0)/std::pow(params_->lambdamax_ - params_->lambda0_,2.0)));

  return;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::EvaluateActiveStressCmatAniso(const LINALG::Matrix<3,3> defgrd,
                                                                      LINALG::Matrix<6,6>& cmat,
                                                                      LINALG::Matrix<6,1>& stress,
                                                                      const int eleGID)
{
  double lambda = 0.0;
  lambda =  A_(0)*defgrd(0,0) + A_(1)*defgrd(1,1) + A_(2)*defgrd(2,2)
      + A_(3)*defgrd(0,1) + A_(4)*defgrd(1,2) + A_(5)*defgrd(0,2)
      + A_(3)*defgrd(1,0) + A_(4)*defgrd(2,1) + A_(5)*defgrd(2,0);

  stress.Update(dPIact_*1./(lambda*lambda),A_,0.0);
  cmat.MultiplyNT(-2.0*dPIact_*1./(lambda*lambda*lambda*lambda),A_,A_,0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::GetDerivativesAniso(LINALG::Matrix<2,1>& dPI_aniso,
                                                            LINALG::Matrix<3,1>& ddPII_aniso,
                                                            LINALG::Matrix<4,1>& dddPIII_aniso,
                                                            const LINALG::Matrix<3,3> rcg,
                                                            const int eleGID)
{
  double I4 = 0.0;
  I4 =  A_(0)*rcg(0,0) + A_(1)*rcg(1,1) + A_(2)*rcg(2,2)
      + A_(3)*rcg(0,1) + A_(4)*rcg(1,2) + A_(5)*rcg(0,2)
      + A_(3)*rcg(1,0) + A_(4)*rcg(2,1) + A_(5)*rcg(2,0);

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  // passive contribution
  dPI_aniso(0) = k1*(I4-1.0)*exp(k2*(I4-1.0)*(I4-1.0));

  ddPII_aniso(0) = (1.0+2.0*k2*(I4-1.0)*(I4-1.0))*k1*exp(k2*(I4-1.0)*(I4-1.0));

  dddPIII_aniso(0) = (3.0 + 2.0*k2*(I4-1.0)*(I4-1.0))*2.0*k1*k2*(I4-1.0)*exp(k2*(I4-1.0)*(I4-1.0));

  return;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::AddStressAnisoPrincipal(const LINALG::Matrix<6,1> rcg,
                                                                LINALG::Matrix<6,6>& cmat,
                                                                LINALG::Matrix<6,1>& stress,
                                                                Teuchos::ParameterList& params,
                                                                const int eleGID)
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
void MAT::ELASTIC::CoupAnisoExpoActive::GetFiberVecs(std::vector<LINALG::Matrix<3,1> >& fibervecs)
{
  fibervecs.push_back(a_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoActive::SetFiberVecs(const double newgamma,
                                                     const LINALG::Matrix<3,3> locsys,
                                                     const LINALG::Matrix<3,3> defgrd)
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

