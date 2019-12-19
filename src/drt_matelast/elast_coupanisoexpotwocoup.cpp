/*----------------------------------------------------------------------*/
/*! \file
\brief the input line should read MAT 1 ELAST_CoupAnisoExpoTwoCoup A4 18.472 B4 16.026 A6 2.481
B6 11.120 A8 0.216 B8 11.436 GAMMA 0.0 [INIT 1] [FIB_COMP Yes] [ADAPT_ANGLE No]

\maintainer Amadeus Gebauer

\level 2

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpotwocoup.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/material_service.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      A4_(matdata->GetDouble("A4")),
      B4_(matdata->GetDouble("B4")),
      A6_(matdata->GetDouble("A6")),
      B6_(matdata->GetDouble("B6")),
      A8_(matdata->GetDouble("A8")),
      B8_(matdata->GetDouble("B8")),
      gamma_(matdata->GetDouble("GAMMA")),
      init_(matdata->GetInt("INIT")),
      fib_comp_(matdata->GetInt("FIB_COMP")),
      adapt_angle_(matdata->GetInt("ADAPT_ANGLE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                                         |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params)
    : Anisotropy(2, params->init_, params->StructuralTensorStrategy()), params_(params), a1a2_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, a1a2_);
  AddtoPack(data, A1A2_);

  PackAnisotropy(data);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a1a2_);
  ExtractfromPack(position, data, A1A2_);

  UnpackAnisotropy(data, position);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  Summand::Setup(numgp, linedef);

  Anisotropy::SetNumberOfGaussPoints(numgp);
  Anisotropy::ReadAnisotropyFromElement(linedef);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::PostSetup(Teuchos::ParameterList& params)
{
  Summand::PostSetup(params);
  Anisotropy::ReadAnisotropyFromParameterList(params);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
    const int eleGID)
{
  int gp = GetGPId(params);

  LINALG::Matrix<6, 1> A1 = GetStructuralTensor_stress(gp, 0);
  LINALG::Matrix<6, 1> A2 = GetStructuralTensor_stress(gp, 1);
  LINALG::Matrix<6, 1> A1A2 = A1A2_[gp];
  double a1a2 = a1a2_[gp];

  double I4 = A1.Dot(rcg);
  double I6 = A2.Dot(rcg);
  double I8;
  I8 = A1A2(0) * rcg(0) + A1A2(1) * rcg(1) + A1A2(2) * rcg(2) + A1A2(3) * rcg(3) +
       A1A2(4) * rcg(4) + A1A2(5) * rcg(5);

  double A4 = params_->A4_;
  double B4 = params_->B4_;
  double A6 = params_->A6_;
  double B6 = params_->B6_;
  double A8 = params_->A8_;
  double B8 = params_->B8_;

  // check if fibers should support compression or not - if not, set the multipliers in front of
  // their strain-energy contribution to zero when the square of their stretches (fiber invariants
  // I4, I6) is smaller than one, respectively - mhv 03/14
  if ((params_->fib_comp_) == 0)
  {
    if (I4 < 1.0) A4 = 0.;
    if (I6 < 1.0) A6 = 0.;
  }

  double gamma = 2.0 * A4 * (I4 - 1.0) * exp(B4 * (I4 - 1.0) * (I4 - 1.0));
  stress.Update(gamma, A1, 1.0);
  gamma = 2.0 * A6 * (I6 - 1.0) * exp(B6 * (I6 - 1.0) * (I6 - 1.0));
  stress.Update(gamma, A2, 1.0);
  gamma = 2.0 * A8 * (I8 - a1a2) * exp(B8 * (I8 - a1a2) * (I8 - a1a2));
  stress.Update(gamma, A1A2, 1.0);

  double delta = 2.0 * (1.0 + 2.0 * B4 * (I4 - 1.0) * (I4 - 1.0)) * 2.0 * A4 *
                 exp(B4 * (I4 - 1.0) * (I4 - 1.0));
  cmat.MultiplyNT(delta, A1, A1, 1.0);
  delta = 2.0 * (1.0 + 2.0 * B6 * (I6 - 1.0) * (I6 - 1.0)) * 2.0 * A6 *
          exp(B6 * (I6 - 1.0) * (I6 - 1.0));
  cmat.MultiplyNT(delta, A2, A2, 1.0);
  delta =
      4.0 * A8 * exp(B8 * (I8 - a1a2) * (I8 - a1a2)) * (1 + 2.0 * B8 * (I8 - a1a2) * (I8 - a1a2));
  cmat.MultiplyNT(delta, A1A2, A1A2, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  if (params_->init_ == INIT_MODE_NODAL_FIBERS)
  {
    // This method expects constant fibers within this element but the init mode is such that
    // fibers are defined on the Gauss points
    // We therefore cannot return sth here.

    // ToDo: This may needs improvements later on if needed!
    return;
  }

  fibervecs.push_back(GetFiber(GPDEFAULT, 0));
  fibervecs.push_back(GetFiber(GPDEFAULT, 1));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  if (params_->init_ != INIT_MODE_EXTERNAL and params_->init_ != INIT_MODE_ELEMENT_FIBERS)
  {
    dserror("Setting the fiber vectors is only possible for external mode and element fibers");
  }
  LINALG::Matrix<3, 1> ca1(true);
  LINALG::Matrix<3, 1> ca2(true);

  // Fiber direction derived from local cosy
  if (params_->init_ == 0 || params_->init_ == 1)
  {
    // alignment angles gamma_i are read from first entry of then unnecessary vectors a1 and a2
    if ((params_->gamma_ < 0) || (params_->gamma_ > 90)) dserror("Fiber angle not in [0,90]");
    // convert
    double gamma = (params_->gamma_ * PI) / 180.;

    if (params_->adapt_angle_ != 0 && newgamma != -1.0)
    {
      if (gamma * newgamma < 0.0)
      {
        gamma = -1.0 * newgamma;
      }
      else
      {
        gamma = newgamma;
      }
    }

    // TODO a1 and a2 aren't orthogonal in the current definition; leads to non-zero value for I8
    // and therefore
    // a non-zero stress value at the reference configuration
    for (int i = 0; i < 3; ++i)
    {
      // a1 = cos gamma e3 + sin gamma e2
      ca1(i) = cos(gamma) * locsys(i, 2) + sin(gamma) * locsys(i, 1);
      // a2 = cos gamma e3 - sin gamma e2
      ca2(i) = cos(gamma) * locsys(i, 2) - sin(gamma) * locsys(i, 1);
    }
  }
  // INIT = 2 : Fiber direction aligned to local cosy
  else if (params_->init_ == 2)
  {
    for (int i = 0; i < 3; ++i)
    {
      ca1(i) = locsys(i, 0);
      ca2(i) = locsys(i, 1);
    }
  }
  else
  {
    dserror("Problem with fiber initialization");
  }

  // pull back in reference configuration
  LINALG::Matrix<3, 1> a1_0(true);
  LINALG::Matrix<3, 1> a2_0(true);
  LINALG::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);

  a1_0.Multiply(idefgrd, ca1);
  a1_0.Scale(1.0 / a1_0.Norm2());

  a2_0.Multiply(idefgrd, ca2);
  a2_0.Scale(1.0 / a2_0.Norm2());

  std::vector<LINALG::Matrix<3, 1>> fibers(0);
  fibers.emplace_back(a1_0);
  fibers.emplace_back(a2_0);

  SetFibers(GPDEFAULT, fibers);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::SetupFiberByCosy(LINALG::Matrix<3, 3>& locsys)
{
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);

  // final setup of fiber data
  SetFiberVecs(0.0, locsys, Id);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::DoFiberInitialization()
{
  // Setup coupling tensor for anisotropy
  if (params_->init_ != 0)
  {
    dserror(
        "The fibers should only be set up manually if INIT is 0, otherwise we (should) get them "
        "from the input file.");
  }

  // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
  LINALG::Matrix<3, 3> Id(false);
  MAT::IdentityMatrix(Id);
  SetFiberVecs(-1.0, Id, Id);
}

void MAT::ELASTIC::CoupAnisoExpoTwoCoup::OnFibersInitialized()
{
  Anisotropy::OnFibersInitialized();

  // Setup structural tensor of the coupling part
  const int fibersperele = GetFibersPerElement();

  A1A2_.resize(fibersperele);
  a1a2_.resize(fibersperele);

  for (int gp = 0; gp < fibersperele; ++gp)
  {
    LINALG::Matrix<3, 1> a1 = GetFiber(gp, 0);
    LINALG::Matrix<3, 1> a2 = GetFiber(gp, 1);
    A1A2_[gp](0) = a1(0) * a2(0);
    A1A2_[gp](1) = a1(1) * a2(1);
    A1A2_[gp](2) = a1(2) * a2(2);
    A1A2_[gp](3) = 0.5 * (a1(0) * a2(1) + a1(1) * a2(0));
    A1A2_[gp](4) = 0.5 * (a1(1) * a2(2) + a1(2) * a2(1));
    A1A2_[gp](5) = 0.5 * (a1(0) * a2(2) + a1(2) * a2(0));

    a1a2_[gp] = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      a1a2_[gp] += a1(i) * a2(i);
    }
  }
}
