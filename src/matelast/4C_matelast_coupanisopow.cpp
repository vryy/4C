/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a pow-like anisotropic material

\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_matelast_coupanisopow.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::ELASTIC::PAR::CoupAnisoPow::CoupAnisoPow(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      k_(matdata->Get<double>("K")),
      d1_(matdata->Get<double>("D1")),
      d2_(matdata->Get<double>("D2")),
      fibernumber_(matdata->Get<int>("FIBER")),
      activethres_(matdata->Get<double>("ACTIVETHRES")),
      gamma_(matdata->Get<double>("GAMMA")),
      init_(matdata->Get<int>("INIT")),
      adapt_angle_(matdata->Get<bool>("ADAPT_ANGLE"))
{
}

MAT::ELASTIC::CoupAnisoPow::CoupAnisoPow(MAT::ELASTIC::PAR::CoupAnisoPow* params) : params_(params)
{
}

void MAT::ELASTIC::CoupAnisoPow::PackSummand(CORE::COMM::PackBuffer& data) const
{
  AddtoPack(data, a_);
  AddtoPack(data, structural_tensor_);
}

void MAT::ELASTIC::CoupAnisoPow::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a_);
  ExtractfromPack(position, data, structural_tensor_);
}

void MAT::ELASTIC::CoupAnisoPow::Setup(int numgp, INPUT::LineDefinition* linedef)
{
  // path if fibers aren't given in .dat file
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    CORE::LINALG::Matrix<3, 3> Id(true);
    for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
    SetFiberVecs(-1.0, Id, Id);
  }

  // path if fibers are given in .dat file
  else if (params_->init_ == 1)
  {
    std::ostringstream ss;
    ss << params_->fibernumber_;
    std::string fibername = "FIBER" + ss.str();  // FIBER Name
    // CIR-AXI-RAD nomenclature
    if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
    {
      // Read in of data
      CORE::LINALG::Matrix<3, 3> locsys(true);
      ReadRadAxiCir(linedef, locsys);
      CORE::LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      // final setup of fiber data
      SetFiberVecs(0.0, locsys, Id);
    }
    // FIBERi nomenclature
    else if (linedef->HaveNamed(fibername))
    {
      // Read in of data
      ReadFiber(linedef, fibername, a_);
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

void MAT::ELASTIC::CoupAnisoPow::add_stress_aniso_principal(const CORE::LINALG::Matrix<6, 1>& rcg,
    CORE::LINALG::Matrix<6, 6>& cmat, CORE::LINALG::Matrix<6, 1>& stress,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // load params
  double k = params_->k_;
  double d1 = params_->d1_;
  double d2 = params_->d2_;
  double activethres = params_->activethres_;

  if (d2 <= 1.0)
  {
    FOUR_C_THROW(
        "exponential factor D2 should be greater than 1.0, since otherwise one can't achieve a "
        "stress free reference state");
  }

  // calc invariant I4
  double I4 = 0.0;
  I4 = structural_tensor_(0) * rcg(0) + structural_tensor_(1) * rcg(1) +
       structural_tensor_(2) * rcg(2) + structural_tensor_(3) * rcg(3) +
       structural_tensor_(4) * rcg(4) + structural_tensor_(5) * rcg(5);

  double lambda4 = pow(I4, 0.5);
  double pow_I4_d1 = pow(I4, d1);
  double pow_I4_d1m1 = pow(I4, d1 - 1.0);
  double pow_I4_d1m2 = pow(I4, d1 - 2.0);
  // Compute stress and material update
  // Beware that the fiber will be turned off in case of compression under activethres.
  // Hence, some compression (i.e. activethres<1.0) could be allow since the fibers are embedded in
  // the matrix and at usually at the microscale not fibers are allowed in the same given direction
  // by FIBER1
  double gamma = 0.0;
  double delta = 0.0;
  if (lambda4 > activethres)
  {
    // Coefficient for residual
    if (pow_I4_d1 > 1.0)
    {
      gamma = 2.0 * k * d2 * d1 * pow_I4_d1m1 * pow(pow_I4_d1 - 1.0, d2 - 1.0);
      // Coefficient for matrix
      delta = 4.0 * k * d2 * (d2 - 1) * d1 * pow_I4_d1m1 * d1 * pow_I4_d1m1 *
                  pow(pow_I4_d1 - 1.0, d2 - 2.0) +
              4.0 * k * d2 * d1 * (d1 - 1.0) * pow_I4_d1m2 * pow(pow_I4_d1 - 1.0, d2 - 1.0);
    }
    else
    {
      gamma = -2.0 * k * d2 * d1 * pow_I4_d1m1 *
              pow(1.0 - pow_I4_d1, d2 - 1.0);  // Note minus sign at the beginning
      // Coefficient for matrix
      delta = 4.0 * k * d2 * (d2 - 1) * d1 * pow_I4_d1m1 * d1 * pow_I4_d1m1 *
                  pow(1.0 - pow_I4_d1, d2 - 2.0) -  // Note minus sign
              4.0 * k * d2 * d1 * (d1 - 1.0) * pow_I4_d1m2 * pow(1.0 - pow_I4_d1, d2 - 1.0);
    }
  }
  stress.Update(gamma, structural_tensor_, 1.0);
  cmat.MultiplyNT(delta, structural_tensor_, structural_tensor_, 1.0);
}

void MAT::ELASTIC::CoupAnisoPow::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

void MAT::ELASTIC::CoupAnisoPow::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
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

  CORE::LINALG::Matrix<3, 1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = cos(gamma) * locsys(i, 2) + sin(gamma) * locsys(i, 1);
  }
  // pull back in reference configuration
  CORE::LINALG::Matrix<3, 1> a_0(true);
  CORE::LINALG::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a_0.Multiply(idefgrd, ca);
  a_.Update(1. / a_0.Norm2(), a_0);

  params_->structural_tensor_strategy()->setup_structural_tensor(a_, structural_tensor_);
}
FOUR_C_NAMESPACE_CLOSE
