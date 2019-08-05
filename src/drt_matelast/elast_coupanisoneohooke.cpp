/*----------------------------------------------------------------------*/
/*!
\brief the input line should read
  MAT 1 CoupAnisoNeoHooke C 100 GAMMA 35.0 INIT 0 ADAPT_ANGLE 0

\level 3

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoneohooke.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoNeoHooke::CoupAnisoNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ParameterAniso(matdata),
      c_(matdata->GetDouble("C")),
      gamma_(matdata->GetDouble("GAMMA")),
      init_(matdata->GetInt("INIT")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoNeoHooke::CoupAnisoNeoHooke(MAT::ELASTIC::PAR::CoupAnisoNeoHooke* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, a_);
  AddtoPack(data, A_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a_);
  ExtractfromPack(position, data, A_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::Setup(DRT::INPUT::LineDefinition* linedef)
{
  // warning message
  std::cout << "Material does not respect a stress free reference state" << std::endl;

  // path if fibers aren't given in .dat file
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
    const int eleGID)
{
  double c = params_->c_;

  double gamma = 2. * c;
  stress.Update(gamma, A_, 1.0);

  // no contribution to cmat
  // double delta = 0.0;
  // cmat.MultiplyNT(delta, A_, A_, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoNeoHooke::SetFiberVecs(
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
