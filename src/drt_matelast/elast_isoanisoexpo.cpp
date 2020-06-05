/*----------------------------------------------------------------------*/
/*! \file
\brief
The input line should read
  MAT 1 ELAST_IsoAnisoExpo K1 10.0 K2 1.0 GAMMA 35.0 K1COMP 0.0 K2COMP 1.0 INIT 0 ADAPT_ANGLE 0

\level 1

\maintainer Amadeus Gebauer
*/


/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoanisoexpo.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/material.H"
#include "../drt_mat/material_service.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoAnisoExpo::IsoAnisoExpo(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ParameterAniso(matdata),
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
MAT::ELASTIC::IsoAnisoExpo::IsoAnisoExpo(MAT::ELASTIC::PAR::IsoAnisoExpo* params) : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpo::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, a_);
  AddtoPack(data, A_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpo::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a_);
  ExtractfromPack(position, data, A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpo::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    LINALG::Matrix<3, 3> Id(true);
    for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
    SetFiberVecs(-1.0, Id, Id);
  }

  // path if fibers are given in .dat file
  else if (params_->init_ == 1)
  {
    // CIR-AXI-RAD nomenclature
    if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
    {
      // Read in of data
      LINALG::Matrix<3, 3> locsys(true);
      ReadRadAxiCir(linedef, locsys);

      LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      // final setup of fiber data
      SetFiberVecs(0.0, locsys, Id);
    }

    // FIBER1 nomenclature
    else if (linedef->HaveNamed("FIBER1"))
    {
      // Read in of data
      ReadFiber(linedef, "FIBER1", a_);
      params_->StructuralTensorStrategy()->SetupStructuralTensor(a_, A_);
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
void MAT::ELASTIC::IsoAnisoExpo::AddStressAnisoModified(const LINALG::Matrix<6, 1>& rcg,
    const LINALG::Matrix<6, 1>& icg, LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress,
    double I3, const int gp, const int eleGID)
{
  double incJ = std::pow(I3, -1.0 / 3.0);  // J^{-2/3}

  double J4 = incJ * (A_(0) * rcg(0) + A_(1) * rcg(1) + A_(2) * rcg(2) + A_(3) * rcg(3) +
                         A_(4) * rcg(4) + A_(5) * rcg(5));  // J4 = J^{-2/3} I4

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (J4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  LINALG::Matrix<6, 1> Saniso(A_);  // first compute Sfbar = 2 dW/dJ4 A_
  double gammabar = 2. * (k1 * (J4 - 1.) * exp(k2 * (J4 - 1.) * (J4 - 1.)));  // 2 dW/dJ4
  Saniso.Scale(gammabar);                                                     // Sfbar

  double traceCSfbar = Saniso(0) * rcg(0) + Saniso(1) * rcg(1) + Saniso(2) * rcg(2) +
                       1. * (Saniso(3) * rcg(3) + Saniso(4) * rcg(4) + Saniso(5) * rcg(5));
  Saniso.Update(-incJ / 3. * traceCSfbar, icg, incJ);

  LINALG::Matrix<6, 6> Psl(true);  // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl, icg, 1.0);
  Psl.MultiplyNT(-1. / 3., icg, icg, 1.0);

  LINALG::Matrix<6, 1> Aiso(A_);
  Aiso.Update(-J4 / 3.0, icg, incJ);
  LINALG::Matrix<6, 6> cmataniso(true);  // isochoric elastic cmat
  double deltabar = 2. * (1. + 2. * k2 * (J4 - 1.) * (J4 - 1.)) * 2. * k1 *
                    exp(k2 * (J4 - 1.) * (J4 - 1.));  // 4 d^2Wf/dJ4dJ4
  cmataniso.MultiplyNT(deltabar, Aiso, Aiso);
  cmataniso.Update(2. / 3. * incJ * traceCSfbar, Psl, 1.0);
  cmataniso.MultiplyNT(-2. / 3., icg, Saniso, 1.0);
  cmataniso.MultiplyNT(-2. / 3., Saniso, icg, 1.0);

  stress.Update(1.0, Saniso, 1.0);
  cmat.Update(1.0, cmataniso, 1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpo::GetDerivativesAniso(LINALG::Matrix<2, 1>& dPI_aniso,
    LINALG::Matrix<3, 1>& ddPII_aniso, LINALG::Matrix<4, 1>& dddPIII_aniso, const double I4,
    const int gp, const int eleGID)
{
  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (I4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }


  dPI_aniso(0) = k1 * (I4 - 1.0) * exp(k2 * (I4 - 1.0) * (I4 - 1.0));

  ddPII_aniso(0) =
      (1.0 + 2.0 * k2 * (I4 - 1.0) * (I4 - 1.0)) * k1 * exp(k2 * (I4 - 1.0) * (I4 - 1.0));

  dddPIII_aniso(0) = (3.0 + 2.0 * k2 * (I4 - 1.0) * (I4 - 1.0)) * 2.0 * k1 * k2 * (I4 - 1.0) *
                     exp(k2 * (I4 - 1.0) * (I4 - 1.0));

  return;
};


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpo::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoAnisoExpo::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  if ((params_->gamma_ < -90) || (params_->gamma_ > 90)) dserror("Fiber angle not in [-90,90]");
  // convert
  double gamma = (params_->gamma_ * PI) / 180.;

  if (params_->adapt_angle_ && newgamma != -1.0)
  {
    if (gamma * newgamma < 0.0)
      gamma = -1.0 * newgamma;
    else
      gamma = newgamma;
  }

  LINALG::Matrix<3, 1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = cos(gamma) * locsys(i, 2) + sin(gamma) * locsys(i, 1);
  }
  // pull back in reference configuration
  LINALG::Matrix<3, 1> a_0(true);
  LINALG::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a_0.Multiply(idefgrd, ca);
  a_.Update(1. / a_0.Norm2(), a_0);

  params_->StructuralTensorStrategy()->SetupStructuralTensor(a_, A_);
}
