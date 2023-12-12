/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a coupled anisotropic neo-Hooke material with one fiber direction with
space-time varying coefficients

\level 3
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_coupanisoneohooke_VarProp.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_globalproblem.H"
#include "baci_mat_par_material.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"
#include "baci_utils_function.H"

BACI_NAMESPACE_OPEN


MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp::CoupAnisoNeoHooke_VarProp(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      c_(matdata->GetDouble("C")),
      sourceactiv_(matdata->GetInt("SOURCE_ACTIVATION")),
      gamma_(matdata->GetDouble("GAMMA")),
      theta_(matdata->GetDouble("THETA")),
      init_(matdata->GetInt("INIT")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}

MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::CoupAnisoNeoHooke_VarProp(
    MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp* params)
    : params_(params)
{
}

void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::PackSummand(CORE::COMM::PackBuffer& data) const
{
  AddtoPack(data, a_);
  AddtoPack(data, A_);
}

void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a_);
  ExtractfromPack(position, data, A_);
}

void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // path if fibers aren't given in .dat file
  if (params_->init_ == 0)
  {
    // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
    CORE::LINALG::Matrix<3, 3> Id(true);
    CORE::LINALG::Matrix<3, 3> locsys(true);
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
        dserror(
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

void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::AddStressAnisoPrincipal(
    const CORE::LINALG::Matrix<6, 1>& rcg, CORE::LINALG::Matrix<6, 6>& cmat,
    CORE::LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  double time_ = params.get<double>("total time", 0.0);
  Teuchos::RCP<std::vector<double>> pos_ =
      params.get<Teuchos::RCP<std::vector<double>>>("position");
  const double* coordgpref_ = pos_->data();
  double stressFact_ =
      DRT::Problem::Instance()
          ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(params_->sourceactiv_ - 1)
          .Evaluate(coordgpref_, time_, 0);


  // double stressFact_=params.get<double>("scalar");
  stress.Update(2 * (params_->c_) * stressFact_, A_, 1.0);

  // no contribution to cmat
  // double delta = 0.0;
  // cmat.MultiplyNT(delta, A_, A_, 1.0);
}

void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

void MAT::ELASTIC::CoupAnisoNeoHooke_VarProp::SetFiberVecs(const double newgamma,
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

BACI_NAMESPACE_CLOSE
