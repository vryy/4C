/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of the isochoric contribution of a anisotropic exponential fiber material

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_isoanisoexpo.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::IsoAnisoExpo::IsoAnisoExpo(const Teuchos::RCP<Core::Mat::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      k1_(matdata->Get<double>("K1")),
      k2_(matdata->Get<double>("K2")),
      gamma_(matdata->Get<double>("GAMMA")),
      k1comp_(matdata->Get<double>("K1COMP")),
      k2comp_(matdata->Get<double>("K2COMP")),
      init_(matdata->Get<int>("INIT")),
      adapt_angle_(matdata->Get<bool>("ADAPT_ANGLE"))
{
}

Mat::Elastic::IsoAnisoExpo::IsoAnisoExpo(Mat::Elastic::PAR::IsoAnisoExpo* params) : params_(params)
{
}

void Mat::Elastic::IsoAnisoExpo::PackSummand(Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, a_);
  add_to_pack(data, structural_tensor_);
}

void Mat::Elastic::IsoAnisoExpo::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  extract_from_pack(position, data, a_);
  extract_from_pack(position, data, structural_tensor_);
}

void Mat::Elastic::IsoAnisoExpo::Setup(int numgp, Input::LineDefinition* linedef)
{
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    Core::LinAlg::Matrix<3, 3> Id(true);
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
      Core::LinAlg::Matrix<3, 3> locsys(true);
      ReadRadAxiCir(linedef, locsys);

      Core::LinAlg::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      // final setup of fiber data
      SetFiberVecs(0.0, locsys, Id);
    }

    // FIBER1 nomenclature
    else if (linedef->HaveNamed("FIBER1"))
    {
      // Read in of data
      ReadFiber(linedef, "FIBER1", a_);
      params_->structural_tensor_strategy()->setup_structural_tensor(a_, structural_tensor_);
    }

    // error path
    else
    {
      FOUR_C_THROW("Reading of element local cosy for anisotropic materials failed");
    }
  }
  else
    FOUR_C_THROW("INIT mode not implemented");
}

void Mat::Elastic::IsoAnisoExpo::add_stress_aniso_modified(const Core::LinAlg::Matrix<6, 1>& rcg,
    const Core::LinAlg::Matrix<6, 1>& icg, Core::LinAlg::Matrix<6, 6>& cmat,
    Core::LinAlg::Matrix<6, 1>& stress, double I3, const int gp, const int eleGID,
    Teuchos::ParameterList& params)
{
  double incJ = std::pow(I3, -1.0 / 3.0);  // J^{-2/3}

  double J4 = incJ * (structural_tensor_(0) * rcg(0) + structural_tensor_(1) * rcg(1) +
                         structural_tensor_(2) * rcg(2) + structural_tensor_(3) * rcg(3) +
                         structural_tensor_(4) * rcg(4) +
                         structural_tensor_(5) * rcg(5));  // J4 = J^{-2/3} I4

  double k1 = params_->k1_;
  double k2 = params_->k2_;

  if (J4 < 1.0)
  {
    k1 = params_->k1comp_;
    k2 = params_->k2comp_;
  }

  Core::LinAlg::Matrix<6, 1> Saniso(structural_tensor_);  // first compute Sfbar = 2 dW/dJ4 A_
  double gammabar = 2. * (k1 * (J4 - 1.) * exp(k2 * (J4 - 1.) * (J4 - 1.)));  // 2 dW/dJ4
  Saniso.Scale(gammabar);                                                     // Sfbar

  double traceCSfbar = Saniso(0) * rcg(0) + Saniso(1) * rcg(1) + Saniso(2) * rcg(2) +
                       1. * (Saniso(3) * rcg(3) + Saniso(4) * rcg(4) + Saniso(5) * rcg(5));
  Saniso.Update(-incJ / 3. * traceCSfbar, icg, incJ);

  Core::LinAlg::Matrix<6, 6> Psl(true);  // Psl = Cinv o Cinv - 1/3 Cinv x Cinv
  AddtoCmatHolzapfelProduct(Psl, icg, 1.0);
  Psl.MultiplyNT(-1. / 3., icg, icg, 1.0);

  Core::LinAlg::Matrix<6, 1> Aiso(structural_tensor_);
  Aiso.Update(-J4 / 3.0, icg, incJ);
  Core::LinAlg::Matrix<6, 6> cmataniso(true);  // isochoric elastic cmat
  double deltabar = 2. * (1. + 2. * k2 * (J4 - 1.) * (J4 - 1.)) * 2. * k1 *
                    exp(k2 * (J4 - 1.) * (J4 - 1.));  // 4 d^2Wf/dJ4dJ4
  cmataniso.MultiplyNT(deltabar, Aiso, Aiso);
  cmataniso.Update(2. / 3. * incJ * traceCSfbar, Psl, 1.0);
  cmataniso.MultiplyNT(-2. / 3., icg, Saniso, 1.0);
  cmataniso.MultiplyNT(-2. / 3., Saniso, icg, 1.0);

  stress.Update(1.0, Saniso, 1.0);
  cmat.Update(1.0, cmataniso, 1.0);
}

void Mat::Elastic::IsoAnisoExpo::GetDerivativesAniso(Core::LinAlg::Matrix<2, 1>& dPI_aniso,
    Core::LinAlg::Matrix<3, 1>& ddPII_aniso, Core::LinAlg::Matrix<4, 1>& dddPIII_aniso,
    const double I4, const int gp, const int eleGID)
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
};

void Mat::Elastic::IsoAnisoExpo::GetFiberVecs(
    std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

void Mat::Elastic::IsoAnisoExpo::SetFiberVecs(const double newgamma,
    const Core::LinAlg::Matrix<3, 3>& locsys, const Core::LinAlg::Matrix<3, 3>& defgrd)
{
  if ((params_->gamma_ < -90) || (params_->gamma_ > 90))
    FOUR_C_THROW("Fiber angle not in [-90,90]");
  // convert
  double gamma = (params_->gamma_ * M_PI) / 180.;

  if (params_->adapt_angle_ && newgamma != -1.0)
  {
    if (gamma * newgamma < 0.0)
      gamma = -1.0 * newgamma;
    else
      gamma = newgamma;
  }

  Core::LinAlg::Matrix<3, 1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = cos(gamma) * locsys(i, 2) + sin(gamma) * locsys(i, 1);
  }
  // pull back in reference configuration
  Core::LinAlg::Matrix<3, 1> a_0(true);
  Core::LinAlg::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a_0.Multiply(idefgrd, ca);
  a_.Update(1. / a_0.Norm2(), a_0);

  params_->structural_tensor_strategy()->setup_structural_tensor(a_, structural_tensor_);
}
FOUR_C_NAMESPACE_CLOSE
