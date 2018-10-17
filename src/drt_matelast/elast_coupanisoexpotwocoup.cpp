/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoexpotwocoup.cpp

\brief the input line should read MAT 1 ELAST_CoupAnisoExpoTwoCoup A4 18.472 B4 16.026 A6 2.481
B6 11.120 A8 0.216 B8 11.436 GAMMA 0.0 [INIT 1] [FIB_COMP Yes] [ADAPT_ANGLE No]

\maintainer Marc Hirschvogel, originally by A. Nagler

\level 2

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpotwocoup.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup::CoupAnisoExpoTwoCoup(
    Teuchos::RCP<MAT::PAR::Material> matdata)
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
    : params_(params), a1a2_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, a1_);
  AddtoPack(data, a2_);
  AddtoPack(data, a1a2_);
  AddtoPack(data, A1_);
  AddtoPack(data, A2_);
  AddtoPack(data, A1A2_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, a1_);
  ExtractfromPack(position, data, a2_);
  ExtractfromPack(position, data, a1a2_);
  ExtractfromPack(position, data, A1_);
  ExtractfromPack(position, data, A2_);
  ExtractfromPack(position, data, A1A2_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::Setup(DRT::INPUT::LineDefinition* linedef)
{
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
      ReadFiber(linedef, "FIBER1", a1_);
      ReadFiber(linedef, "FIBER2", a2_);
      params_->StructuralTensorStrategy()->SetupStructuralTensor(a1_, A1_);
      params_->StructuralTensorStrategy()->SetupStructuralTensor(a2_, A2_);
      A1A2_.MultiplyNT(a1_, a2_);
    }

    // error path
    else
    {
      dserror("Reading of element local cosy for anisotropic materials failed");
    }
  }
  // fibers defined on nodes
  else if (params_->init_ == 3)
  {
    // nothing to do here. gp fibers are passed from element at evluate
  }
  else
    dserror("INIT mode not implemented");

  a1a2_ = a1_(0) * a2_(0) + a1_(1) * a2_(1) + a1_(2) * a2_(2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::AddStressAnisoPrincipal(const LINALG::Matrix<6, 1>& rcg,
    LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress, Teuchos::ParameterList& params,
    const int eleGID)
{
  if (params_->init_ == 3)
  {
    if (params.isParameter("gpfiber1") && params.isParameter("gpfiber2"))
    {
      a1_ = params.get<LINALG::Matrix<3, 1>>("gpfiber1");
      a2_ = params.get<LINALG::Matrix<3, 1>>("gpfiber2");
      params_->StructuralTensorStrategy()->SetupStructuralTensor(a1_, A1_);
      params_->StructuralTensorStrategy()->SetupStructuralTensor(a2_, A2_);
      A1A2_.MultiplyNT(a1_, a2_);
      a1a2_ = a1_(0) * a2_(0) + a1_(1) * a2_(1) + a1_(2) * a2_(2);
    }
    else
      dserror("No fiber at gauss point available.");
  }

  double I4 = 0.0;
  I4 = A1_(0) * rcg(0) + A1_(1) * rcg(1) + A1_(2) * rcg(2) + A1_(3) * rcg(3) + A1_(4) * rcg(4) +
       A1_(5) * rcg(5);
  double I6 = 0.0;
  I6 = A2_(0) * rcg(0) + A2_(1) * rcg(1) + A2_(2) * rcg(2) + A2_(3) * rcg(3) + A2_(4) * rcg(4) +
       A2_(5) * rcg(5);
  double I8 = 0.0;
  I8 = A1A2_(0, 0) * rcg(0) + A1A2_(1, 1) * rcg(1) + A1A2_(2, 2) * rcg(2) +
       0.5 * (A1A2_(0, 1) * rcg(3) + A1A2_(1, 2) * rcg(4) + A1A2_(0, 2) * rcg(5)) +
       0.5 * (A1A2_(1, 0) * rcg(3) + A1A2_(2, 1) * rcg(4) + A1A2_(2, 0) * rcg(5));

  // build Voigt (stress-like) version of a1 \otimes a2 + a2 \otimes a1
  LINALG::Matrix<6, 1> A1A2sym;
  A1A2sym(0) = 2 * A1A2_(0, 0);
  A1A2sym(1) = 2 * A1A2_(1, 1);
  A1A2sym(2) = 2 * A1A2_(2, 2);
  A1A2sym(3) = A1A2_(0, 1) + A1A2_(1, 0);
  A1A2sym(4) = A1A2_(1, 2) + A1A2_(2, 1);
  A1A2sym(5) = A1A2_(0, 2) + A1A2_(2, 0);

  double A4 = params_->A4_;
  double B4 = params_->B4_;
  double A6 = params_->A6_;
  double B6 = params_->B6_;
  double A8 = params_->A8_;
  double B8 = params_->B8_;

  // check if fibers should support compression or not - if not, set the multipliers infront of
  // their strain-energy contribution to zero when the square of their stretches (fiber invariants
  // I4, I6) is smaller than one, respectively - mhv 03/14
  if (!(params_->fib_comp_))
  {
    if (I4 < 1.0) A4 = 0.;
    if (I6 < 1.0) A6 = 0.;
  }

  double gamma = 2.0 * A4 * (I4 - 1.0) * exp(B4 * (I4 - 1.0) * (I4 - 1.0));
  stress.Update(gamma, A1_, 1.0);
  gamma = 2.0 * A6 * (I6 - 1.0) * exp(B6 * (I6 - 1.0) * (I6 - 1.0));
  stress.Update(gamma, A2_, 1.0);
  gamma = A8 * (I8 - a1a2_) * exp(B8 * (I8 - a1a2_) * (I8 - a1a2_));
  stress.Update(gamma, A1A2sym, 1.0);

  double delta = 2.0 * (1.0 + 2.0 * B4 * (I4 - 1.0) * (I4 - 1.0)) * 2.0 * A4 *
                 exp(B4 * (I4 - 1.0) * (I4 - 1.0));
  cmat.MultiplyNT(delta, A1_, A1_, 1.0);
  delta = 2.0 * (1.0 + 2.0 * B6 * (I6 - 1.0) * (I6 - 1.0)) * 2.0 * A6 *
          exp(B6 * (I6 - 1.0) * (I6 - 1.0));
  cmat.MultiplyNT(delta, A2_, A2_, 1.0);
  delta = A8 * exp(B8 * (I8 - a1a2_) * (I8 - a1a2_)) * (1 + 2.0 * B8 * (I8 - a1a2_) * (I8 - a1a2_));
  cmat.MultiplyNT(delta, A1A2sym, A1A2sym, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs  ///< vector of all fiber vectors
)
{
  fibervecs.push_back(a1_);
  fibervecs.push_back(a2_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupAnisoExpoTwoCoup::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  LINALG::Matrix<3, 1> ca1(true);
  LINALG::Matrix<3, 1> ca2(true);

  // Fiber direction derived from local cosy
  if (params_->init_ == 0 || params_->init_ == 1)
  {
    // alignment angles gamma_i are read from first entry of then unnecessary vectors a1 and a2
    if ((params_->gamma_ < 0) || (params_->gamma_ > 90)) dserror("Fiber angle not in [0,90]");
    // convert
    double gamma = (params_->gamma_ * PI) / 180.;

    if (params_->adapt_angle_ && newgamma != -1.0)
    {
      if (gamma * newgamma < 0.0)
        gamma = -1.0 * newgamma;
      else
        gamma = newgamma;
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
  a1_.Update(1. / a1_0.Norm2(), a1_0);
  a2_0.Multiply(idefgrd, ca2);
  a2_.Update(1. / a2_0.Norm2(), a2_0);

  params_->StructuralTensorStrategy()->SetupStructuralTensor(a1_, A1_);
  params_->StructuralTensorStrategy()->SetupStructuralTensor(a2_, A2_);
  A1A2_.MultiplyNT(a1_, a2_);
}
