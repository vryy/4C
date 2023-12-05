/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a coupled anisotropic neo-Hooke material with one fiber direction

\level 3
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_coupanisoneohooke.H"

#include "baci_io_linedefinition.H"
#include "baci_mat_par_material.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"


MAT::ELASTIC::PAR::CoupAnisoNeoHooke::CoupAnisoNeoHooke(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      c_(matdata->GetDouble("C")),
      gamma_(matdata->GetDouble("GAMMA")),
      init_(matdata->GetInt("INIT")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}

MAT::ELASTIC::CoupAnisoNeoHooke::CoupAnisoNeoHooke(MAT::ELASTIC::PAR::CoupAnisoNeoHooke* params)
    : params_(params)
{
}

void MAT::ELASTIC::CoupAnisoNeoHooke::PackSummand(CORE::COMM::PackBuffer& data) const
{
  AddtoPack(data, a_);
  AddtoPack(data, A_);
}

void MAT::ELASTIC::CoupAnisoNeoHooke::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a_);
  ExtractfromPack(position, data, A_);
}

void MAT::ELASTIC::CoupAnisoNeoHooke::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // warning message
  std::cout << "Material does not respect a stress free reference state" << std::endl;

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

    // FIBER1 nomenclature
    else if (linedef->HaveNamed("FIBER1"))
    {
      // Read in of fiber data and setting fiber data
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

void MAT::ELASTIC::CoupAnisoNeoHooke::AddStressAnisoPrincipal(const CORE::LINALG::Matrix<6, 1>& rcg,
    CORE::LINALG::Matrix<6, 6>& cmat, CORE::LINALG::Matrix<6, 1>& stress,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  double c = params_->c_;

  double gamma = 2. * c;
  stress.Update(gamma, A_, 1.0);

  // no contribution to cmat
  // double delta = 0.0;
  // cmat.MultiplyNT(delta, A_, A_, 1.0);
}

void MAT::ELASTIC::CoupAnisoNeoHooke::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

void MAT::ELASTIC::CoupAnisoNeoHooke::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  if ((params_->gamma_ < -90) || (params_->gamma_ > 90)) dserror("Fiber angle not in [-90,90]");
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

  params_->StructuralTensorStrategy()->SetupStructuralTensor(a_, A_);
}