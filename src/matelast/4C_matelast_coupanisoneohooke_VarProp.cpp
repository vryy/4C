/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a coupled anisotropic neo-Hooke material with one fiber direction with
space-time varying coefficients

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoneohooke_VarProp.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::Elastic::PAR::CoupAnisoNeoHookeVarProp::CoupAnisoNeoHookeVarProp(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ParameterAniso(matdata),
      c_(matdata.parameters.get<double>("C")),
      sourceactiv_(matdata.parameters.get<int>("SOURCE_ACTIVATION")),
      gamma_(matdata.parameters.get<double>("GAMMA")),
      theta_(matdata.parameters.get<double>("THETA")),
      init_(matdata.parameters.get<int>("INIT")),
      adapt_angle_(matdata.parameters.get<bool>("ADAPT_ANGLE"))
{
}

Mat::Elastic::CoupAnisoNeoHookeVarProp::CoupAnisoNeoHookeVarProp(
    Mat::Elastic::PAR::CoupAnisoNeoHookeVarProp* params)
    : params_(params)
{
}

void Mat::Elastic::CoupAnisoNeoHookeVarProp::PackSummand(
    Core::Communication::PackBuffer& data) const
{
  add_to_pack(data, a_);
  add_to_pack(data, structural_tensor_);
}

void Mat::Elastic::CoupAnisoNeoHookeVarProp::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  extract_from_pack(position, data, a_);
  extract_from_pack(position, data, structural_tensor_);
}

void Mat::Elastic::CoupAnisoNeoHookeVarProp::setup(int numgp, Input::LineDefinition* linedef)
{
  // path if fibers aren't given in .dat file
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    Core::LinAlg::Matrix<3, 3> Id(true);
    Core::LinAlg::Matrix<3, 3> locsys(true);
    for (int i = 0; i < 3; i++) Id(i, i) = 1.0;

    // To realize a full rotated fiber orientation and to keep the general structure of
    // SetFiberVec() the input of locsys has to be adapted if one sets
    //               1           0                 sin(theta_)
    //  locsys := [  0       sin(theta_)               0                  ]
    //               0  sin(gamma)*cos(theta_)    cos(gamma_)*cos(theta_)
    // The call of SetFiberVec() will leed to the following fiber direction
    // a = cos(gamma_)*locsys(:,2) + sin(gamma_)*locsys(:,1)
    //         cos(gamma_)*sin(theta_)               0 cos(gamma_)*sin(theta_)
    //   = [             0              ] + [  sin(gamma_)*sin(theta_)  ] = [
    //   sin(gamma_)*sin(theta_) ] =: sperical coordinates
    //         cos(gamma)^2*cos(theta_)        sin(gamma_)^2*cos(theta_)               cos(theta_)
    //
    {
      // Local initialization of spherical angles
      double theta = (params_->theta_);
      double gamma = (params_->gamma_);
      if (gamma < 0.0 || gamma > 180.0 || abs(theta) > 180.0)
      {
        FOUR_C_THROW(
            "Wrong choice of sherical coodinates. Correct domain is gamma in [0,180], theta in "
            "[-180, 180]");
      }
      // conversion to radian measure
      theta = (theta * M_PI) / 180.0;
      gamma = (gamma * M_PI) / 180.0;
      locsys(1, 1) = sin(theta);
      locsys(2, 1) = sin(gamma) * cos(theta);
      locsys(0, 2) = sin(theta);
      locsys(2, 2) = cos(gamma) * cos(theta);
    }
    SetFiberVecs(-1.0, locsys, Id);
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

void Mat::Elastic::CoupAnisoNeoHookeVarProp::add_stress_aniso_principal(
    const Core::LinAlg::Matrix<6, 1>& rcg, Core::LinAlg::Matrix<6, 6>& cmat,
    Core::LinAlg::Matrix<6, 1>& stress, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  double time_ = params.get<double>("total time", 0.0);
  const auto& element_center_coordinates_ref =
      params.get<Core::LinAlg::Matrix<3, 1>>("elecenter_coords_ref");
  double stressFact_ =
      Global::Problem::Instance()
          ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(params_->sourceactiv_ - 1)
          .evaluate(element_center_coordinates_ref.data(), time_, 0);


  // double stressFact_=params.get<double>("scalar");
  stress.update(2 * (params_->c_) * stressFact_, structural_tensor_, 1.0);

  // no contribution to cmat
  // double delta = 0.0;
  // cmat.multiply_nt(delta, A_, A_, 1.0);
}

void Mat::Elastic::CoupAnisoNeoHookeVarProp::GetFiberVecs(
    std::vector<Core::LinAlg::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

void Mat::Elastic::CoupAnisoNeoHookeVarProp::SetFiberVecs(const double newgamma,
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
