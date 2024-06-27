/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a coupled anisotropic neo-Hooke material with one fiber direction

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoneohooke.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::CoupAnisoNeoHooke::CoupAnisoNeoHooke(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ParameterAniso(matdata),
      c_(matdata.parameters.get<double>("C")),
      gamma_(matdata.parameters.get<double>("GAMMA")),
      init_(matdata.parameters.get<int>("INIT")),
      adapt_angle_(matdata.parameters.get<bool>("ADAPT_ANGLE"))
{
}

Mat::Elastic::CoupAnisoNeoHooke::CoupAnisoNeoHooke(Mat::Elastic::PAR::CoupAnisoNeoHooke* params)
    : params_(params)
{
}

void Mat::Elastic::CoupAnisoNeoHooke::PackSummand(Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, a_);
  add_to_pack(data, structural_tensor_);
}

void Mat::Elastic::CoupAnisoNeoHooke::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  extract_from_pack(position, data, a_);
  extract_from_pack(position, data, structural_tensor_);
}

void Mat::Elastic::CoupAnisoNeoHooke::setup(int numgp, Input::LineDefinition* linedef)
{
  // warning message
  std::cout << "Material does not respect a stress free reference state" << std::endl;

  // path if fibers aren't given in .dat file
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
    if (linedef->has_named("RAD") and linedef->has_named("AXI") and linedef->has_named("CIR"))
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
    else if (linedef->has_named("FIBER1"))
    {
      // Read in of fiber data and setting fiber data
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

void Mat::Elastic::CoupAnisoNeoHooke::add_stress_aniso_principal(
    const Core::LinAlg::Matrix<6, 1>& rcg, Core::LinAlg::Matrix<6, 6>& cmat,
    Core::LinAlg::Matrix<6, 1>& stress, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  double c = params_->c_;

  double gamma = 2. * c;
  stress.update(gamma, structural_tensor_, 1.0);

  // no contribution to cmat
  // double delta = 0.0;
  // cmat.multiply_nt(delta, A_, A_, 1.0);
}

void Mat::Elastic::CoupAnisoNeoHooke::GetFiberVecs(
    std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

void Mat::Elastic::CoupAnisoNeoHooke::SetFiberVecs(const double newgamma,
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
  idefgrd.invert(defgrd);

  a_0.multiply(idefgrd, ca);
  a_.update(1. / a_0.norm2(), a_0);

  params_->structural_tensor_strategy()->setup_structural_tensor(a_, structural_tensor_);
}
FOUR_C_NAMESPACE_CLOSE
