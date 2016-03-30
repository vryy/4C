/*----------------------------------------------------------------------*/
/*!
\file elast_isoanisoexpodispersion.cpp
\brief


the input line should read
  MAT 1 ELAST_IsoAnisoExpoDispersion K1 10.0 K2 1.0 GAMMA 35.0  K1COMP 0.0 K2COMP 1.0 [INIT 1] [ADAPT_ANGLE No] KAPPA 0.0


// This worked on the small example testcase_aspect10f_haskettD.dat, but then somehow didn't work on a bigger example.
// I think it still lies with how the fibre-contribution is excluded.
// If J4 < 1
// -- switching k1, k2 to 0 is wrong - see Matlab figures (WhereAreYouGaussPoint.m)
// -- and currently not sure about making J4 == 1, what should happen in R, and what happens with that other weird projection
// term that also includes J4?


<pre>
\maintainer Hamman de Vaal
            devaal@lnm.mw.tum.de
            089/289 15251
</pre>
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoanisoexpodispersion.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/material.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoAnisoExpoDispersion::IsoAnisoExpoDispersion(
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
  kappa_(matdata->GetDouble("KAPPA"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   st         03/12 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoAnisoExpoDispersion::IsoAnisoExpoDispersion(MAT::ELASTIC::PAR::IsoAnisoExpoDispersion* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpoDispersion::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,a_);
  AddtoPack(data,A_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpoDispersion::UnpackSummand(
  const std::vector<char>& data,
  std::vector<char>::size_type& position
  )
{
  ExtractfromPack(position,data,a_);
  ExtractfromPack(position,data,A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpoDispersion::Setup(DRT::INPUT::LineDefinition* linedef)
{
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3,3> Id(true);
    for (int i=0; i<3; i++)
      Id(i,i) = 1.0;
    SetFiberVecs(-1.0,Id,Id);
  }
  else if (params_->init_ == 1)
  {
    // fibers aligned in local element cosy with gamma around circumferential direction
    // -> check whether element supports local element cosy
    if (linedef->HaveNamed("RAD") and
        linedef->HaveNamed("AXI") and
        linedef->HaveNamed("CIR"))
    {
      std::vector<double> rad;
      std::vector<double> axi;
      std::vector<double> cir;
      // read local (cylindrical) cosy-directions at current element
      // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
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
      SetFiberVecs(-1.0,locsys,Id);
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
void MAT::ELASTIC::IsoAnisoExpoDispersion::AddStressAnisoModified(
    const LINALG::Matrix<6,1> rcg,
    LINALG::Matrix<6,1> icg,
    LINALG::Matrix<6,6>& cmat,
    LINALG::Matrix<6,1>& stress,
    double I3,
    const int eleGID
)
{
  double incJ = std::pow(I3,-1./3.);  // J^{-2/3}

  // modified 1st and 4th invariants
  double J1 = incJ*(rcg(0) + rcg(1) + rcg(2)); //J1 = J^{-2/3} I1
  // just temporary for debugging purposes
  double I4 = 0.0;
  I4 =  A_(0)*rcg(0) + A_(1)*rcg(1) + A_(2)*rcg(2)
      + A_(3)*rcg(3) + A_(4)*rcg(4) + A_(5)*rcg(5);
  //double J4 = incJ * ( A_(0)*rcg(0) + A_(1)*rcg(1) + A_(2)*rcg(2) + A_(3)*rcg(3) + A_(4)*rcg(4) + A_(5)*rcg(5)); //J4 = J^{-2/3} I4
  double J4 = incJ * I4;

  double k1 = params_->k1_;
  double k2 = params_->k2_;
  double kappa = params_->kappa_;

  if (kappa == 0.33) kappa = 1./3.;
  else if (kappa > 0.33) dserror("kappa must be [0:1/3]. By setting KAPPA 0.33, you will get a KAPPA of (1./3.)");
  else if (kappa < 0.0) dserror("kappa must be [0:1/3]");

  // The inbetween term when differentiating the exponent
  double R;

  if (J4 < 1.0)
  {
    //TODO: confirm how the fibre-contribution is neglected when stretch in fibre direction < 1
    // At this stage, the compression fibre-stiffness parameters are not considered here as in elast_isoanisoexpo,
    // but rather the J4 contribution is zeroed out
    // It would be wrong to allow zero k1 and k2 when kappa > 0, because it would then cause a
    // discontinuity in the SEF when J1 != 3 at the border between J4 < 1 and J4 > 1.
    // How exactly to do this is not clearly defined.

    //k1 = params_->k1comp_;
    //k2 = params_->k2comp_;

    R = kappa*J1 + (1. - 3.*kappa)*1. - 1.;
    //R = kappa*J1 - 1.; // this does not work
    //J4 = 1.; // not sure if this should be the case
  }
  else
  {
      // The inbetween term when differentiating the exponent
      R = kappa*J1 + (1. - 3.*kappa)*J4 - 1.;
  }

  // The exponent term: Q = k2*(((kappa*I1 + (1. - 3.*kappa)*I4) - 1.)*((kappa*I1 + (1. - 3.*kappa)*I4) - 1.));
  double Q = exp(k2*R*R);
  // Term reoccurring in Cbar coefficients
  double V = 2.*k2*R*R + 1.;


  // The gamma and delta calculations for the anisotropic parts are for the time-being
  // located in the material itself, and not within elasthyper.cpp
  LINALG::Matrix<2,1> gammabar(true);
  LINALG::Matrix<4,1> deltabar(true);

  // build Cartesian identity 2-tensor I_{AB} -- could perhaps read in from somewhere else
  LINALG::Matrix<6,1> id2_temp(true);
  for (int i=0; i<3; i++) id2_temp(i) = 1.0;

  // 1st invariant, trace -- could perhaps read in from somewhere else
  const double I1 = rcg(0) + rcg(1) + rcg(2);

  /*std::cout << "\n" << std::endl;
  std::cout << "this is rcg(0) and rcg(1) and rcg(2): " << rcg(0) << " " << rcg(1) << " " << rcg(2) << std::endl;
  std::cout << "this is k1, k2 and kappa: " << k1 << " " << k2 << " " << kappa << std::endl;
  std::cout << "this is J1, J4: " << J1 << " " << J4 << std::endl;
  std::cout << "this is I1, I4: " << I1 << " " << I4 << std::endl;
  std::cout << "this is R and Q: " << R << " " << Q << std::endl;*/

  if (std::isinf(Q))
  {
    std::cout << "\n" << std::endl;
    std::cout << "this is rcg(0) and rcg(1) and rcg(2): " << rcg(0) << " " << rcg(1) << " " << rcg(2) << std::endl;
    std::cout << "this is k1, k2 and kappa: " << k1 << " " << k2 << " " << kappa << std::endl;
    std::cout << "this is J1, J4: " << J1 << " " << J4 << std::endl;
    std::cout << "this is I1, I4: " << I1 << " " << I4 << std::endl;
    std::cout << "this is R and Q: " << R << " " << Q << std::endl;
    dserror("The exponential Q became infinite");
  }
  if (std::isnan(Q))
  {
    std::cout << "\n" << std::endl;
    std::cout << "this is rcg(0) and rcg(1) and rcg(2): " << rcg(0) << " " << rcg(1) << " " << rcg(2) << std::endl;
    std::cout << "this is k1, k2 and kappa: " << k1 << " " << k2 << " " << kappa << std::endl;
    std::cout << "this is I1, I4: " << I1 << " " << I4 << std::endl;
    std::cout << "this is R and Q: " << R << " " << Q << std::endl;
    dserror("The exponential Q is NaN");
  }

  // Sfbar = 2 dW/dJ1 1_ + 2 dW/dJ4 A_
  // first compute scalar coefficients of Sfbar
  gammabar(0) = 2.*k1*kappa        *R*Q; // gammabar(0) = 2 dW/dJ1
  gammabar(1) = 2.*k1*(1.-3.*kappa)*R*Q; // gammabar(1) = 2 dW/dJ4

  // Now build Sfbar
  LINALG::Matrix<6,1> Saniso(A_);
  Saniso.Update(gammabar(0),id2_temp,gammabar(1)); //Sfbar = 2 dW/dJ1 id2 + 2 dW/dJ4 A_

  double traceCSfbar =  Saniso(0)*rcg(0) + Saniso(1)*rcg(1) + Saniso(2)*rcg(2)
                      + 1.*(Saniso(3)*rcg(3) + Saniso(4)*rcg(4) + Saniso(5)*rcg(5));
  Saniso.Update(-incJ/3.*traceCSfbar,icg,incJ); //S_isoch (eq. 6.90)

  LINALG::Matrix<6,6>  Psl(true);        // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl,icg,1.0);
  Psl.MultiplyNT(-1./3.,icg,icg,1.0);

  // Cbar/J^{-4/3} = 4 d(dW)/(dJ1*dJ1) id2 \otimes id2 + 4 d(dW)/(dJ1*dJ4) A_ \otimes id2
  //        + 4 d(dW)/(dJ4*dJ1) id2 \otimes A_ + 4 d(dW)/(dJ4*dJ4) A_ \otimes A_
  // first compute scalar coefficients of Cbar
  deltabar(0) = 4.*k1*kappa        *kappa      *Q*V; // 4 incJ*incJ d(dW)/(dJ1*dJ1)
  deltabar(1) = 4.*k1*(1.-3.*kappa)*kappa      *Q*V; // 4 d(dW)/(dJ1*dJ4)
  deltabar(2) = deltabar(1);                                  // 4 d(dW)/(dJ4*dJ1)
  deltabar(3) = 4.*k1*(1.-3.*kappa)*(1-3*kappa)*Q*V; // 4 d(dW)/(dJ4*dJ4)

  /*std::cout << "gammabar(0), gammabar(1): " << gammabar(0) << " " << gammabar(1) << std::endl;
  std::cout << "deltabar(0), deltabar(1), deltabar(2), deltabar(3): " << deltabar(0) << " " << deltabar(1) << " " << deltabar(2)<< " " << deltabar(3) << std::endl;*/

  // A_hat (see derivation)
  LINALG::Matrix<6,1> Aiso(A_);
  Aiso.Update(-J4/3.0,icg,incJ);
  // id2_hat (see derivation)
  LINALG::Matrix<6,1> id2_hat(id2_temp);
  id2_hat.Update(-J1/3.,icg,incJ);

  // Now build Cbar
  LINALG::Matrix<6,6> cmataniso(true); // isochoric elastic cmat

  // term: |P : Cbar : |P^T
  // contribution: deltabar(0) (id2_hat \otimes id2_hat)
  cmataniso.MultiplyNT(deltabar(0),id2_hat,id2_hat);
  // contribution: deltabar(1) (A_hat \otimes id2_hat)
  cmataniso.MultiplyNT(deltabar(1),Aiso,id2_hat,1.0);
  // contribution: deltabar(2) (id2_temp \otimes A_hat)
  cmataniso.MultiplyNT(deltabar(2),id2_hat,Aiso,1.0);
  // contribution: deltabar(3) (A_hat \otimes A_hat)
  cmataniso.MultiplyNT(deltabar(3),Aiso,Aiso,1.0);

  // term: (2/3) Tr(incJ Sfbar) Psl -- (same as in elast_isoanisoexpo.cpp)
  cmataniso.Update(2./3.*incJ*traceCSfbar,Psl,1.0);

  // term: -(2/3) (icg \otimes S_isoch + S_isoch \otimes icg) -- (same as in elast_isoanisoexpo.cpp)
  cmataniso.MultiplyNT(-2./3.,icg,Saniso,1.0);
  cmataniso.MultiplyNT(-2./3.,Saniso,icg,1.0);

  stress.Update(1.0,Saniso,1.0);
  cmat.Update(1.0,cmataniso,1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpoDispersion::SetFiberVecs(
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
void MAT::ELASTIC::IsoAnisoExpoDispersion::GetFiberVecs(
    std::vector<LINALG::Matrix<3,1> >& fibervecs ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}
