/*----------------------------------------------------------------------*/
/*!
\file elast_couptransverselyisotropic.cpp

\brief Summand for the transversely hyperelastic isotropic material model
\level 3

\maintainer Michael Hiermeier

\date Sep 20, 2017

*/
/*----------------------------------------------------------------------*/

#include "elast_couptransverselyisotropic.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

#include "../drt_mat/material_service.H"
#include "../drt_io/io_pstream.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupTransverselyIsotropic::CoupTransverselyIsotropic(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : ParameterAniso(matdata),
      alpha_(matdata->GetDouble("ALPHA")),
      beta_(matdata->GetDouble("BETA")),
      gamma_(matdata->GetDouble("GAMMA")),
      angle_(matdata->GetDouble("ANGLE")),
      fiber_gid_(matdata->GetInt("FIBER")),
      init_(matdata->GetInt("INIT"))
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::PAR::CoupTransverselyIsotropic::Print() const
{
  IO::cout << "--- Material parameters of CoupTransverselyIsotropic\n";
  IO::cout << "ALPHA           = " << alpha_ << IO::endl;
  IO::cout << "BETA            = " << beta_ << IO::endl;
  IO::cout << "GAMMA           = " << gamma_ << IO::endl;
  IO::cout << "ANGLE           = " << angle_ << IO::endl;
  IO::cout << "GLOBAL FIBER ID = " << fiber_gid_ << IO::endl;
  IO::cout << "INIT            = " << init_ << IO::endl;
  IO::cout << "----------------------------------------------------\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
MAT::ELASTIC::CoupTransverselyIsotropic::CoupTransverselyIsotropic(my_params* params)
    : params_(params)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::Setup(DRT::INPUT::LineDefinition* linedef)
{
  switch (params_->init_)
  {
    // path if fibers aren't given in .dat file
    case 0:
    {
      // fibers aligned in YZ-plane with gamma around Z in global cartesian cosy
      LINALG::Matrix<3, 3> Id(true);
      for (int i = 0; i < 3; i++) Id(i, i) = 1.0;
      SetFiberVecs(-1.0, Id, Id);

      break;
    }
    // path if fibers are given in .dat file
    case 1:
    {
      std::ostringstream ss;
      ss << params_->fiber_gid_;
      std::string fibername = "FIBER" + ss.str();  // FIBER Name
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
      // FIBERi nomenclature
      else if (linedef->HaveNamed(fibername))
      {
        // Read in of data
        ReadFiber(linedef, fibername, A_);
        params_->StructuralTensorStrategy()->SetupStructuralTensor(A_, AA_);
      }
      // error path
      else
      {
        dserror("Reading of element local cosy for anisotropic materials failed");
      }

      break;
    }
    default:
    {
      dserror("INIT mode not implemented");
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, A_);
  AddtoPack(data, AA_);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, A_);
  ExtractfromPack(position, data, AA_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::GetFiberVecs(
    std::vector<LINALG::Matrix<3, 1>>& fibervecs)
{
  fibervecs.push_back(A_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  if ((params_->angle_ < -90) || (params_->angle_ > 90))
    dserror("Fiber angle not in [-90,90]! Given angle = %d", params_->angle_);
  // convert
  const double angle = (params_->gamma_ * PI) / 180.;

  LINALG::Matrix<3, 1> ca(true);
  for (int i = 0; i < 3; ++i)
  {
    // a = cos gamma e3 + sin gamma e2
    ca(i) = std::cos(angle) * locsys(i, 2) + std::sin(angle) * locsys(i, 1);
  }
  // pull back in reference configuration
  LINALG::Matrix<3, 1> A_0(true);
  LINALG::Matrix<3, 3> idefgrd(true);
  idefgrd.Invert(defgrd);

  A_0.Multiply(idefgrd, ca);
  A_.Update(1. / A_0.Norm2(), A_0);

  params_->StructuralTensorStrategy()->SetupStructuralTensor(A_, AA_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::AddStrainEnergy(double& psi,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& modinv,
    const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<6, 1> identity(true);
  std::fill(identity.A(), identity.A() + 3, 1.0);

  // convert Green-Lagrange strain to right Cauchy-Green Tensor
  // C_{AB} = 2 * E_{AB} + I_{AB} [ REMARK: strain-like 6-Voigt vector ]
  LINALG::Matrix<6, 1> rcg(true);
  rcg.Update(2.0, glstrain, 1.0);
  rcg.Update(1.0, identity, 1.0);

  ResetInvariants(rcg);

  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;

  psi += (alpha + 0.5 * beta * std::log(prinv(2)) + gamma * (I4_ - 1.0)) * (I4_ - 1.0) -
         0.5 * alpha * (I5_ - 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::AddStressAnisoPrincipal(
    const LINALG::Matrix<6, 1>& rcg, LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& stress,
    Teuchos::ParameterList& params, const int eleGID)
{
  ResetInvariants(rcg);

  // switch to stress notation
  LINALG::Matrix<6, 1> rcg_s(rcg);
  MAT::ConvertStrainToStressNotation(rcg_s);

  LINALG::Matrix<6, 1> rcg_inv_s(false);
  UpdateSecondPiolaKirchhoffStress(stress, rcg_s, rcg_inv_s);

  UpdateElasticityTensor(cmat, rcg_inv_s);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::UpdateElasticityTensor(
    LINALG::Matrix<6, 6>& cmat, const LINALG::Matrix<6, 1>& rcg_inv_s) const
{
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;

  const double identity[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};

  // (0) contribution
  {
    const double delta = 8.0 * gamma;
    cmat.MultiplyNT(delta, AA_, AA_, 1.0);
  }

  // (1) contribution
  {
    const double delta = 2.0 * beta;
    cmat.MultiplyNT(delta, rcg_inv_s, AA_, 1.0);
    cmat.MultiplyNT(delta, AA_, rcg_inv_s, 1.0);
  }

  // (2) contribution
  {
    const double delta = -alpha;
    for (unsigned a = 0; a < 6; ++a)
    {
      for (unsigned b = 0; b < 6; ++b)
      {
        const unsigned i = MAT::IMap::fourth_[a][b][0];
        const unsigned j = MAT::IMap::fourth_[a][b][1];
        const unsigned k = MAT::IMap::fourth_[a][b][2];
        const unsigned l = MAT::IMap::fourth_[a][b][3];

        cmat(a, b) += delta * (A_(i) * A_(l) * identity[MAT::IMap::second_[j][k]] +
                                  A_(i) * A_(k) * identity[MAT::IMap::second_[j][l]] +
                                  A_(k) * A_(j) * identity[MAT::IMap::second_[i][l]] +
                                  A_(l) * A_(j) * identity[MAT::IMap::second_[i][k]]);
      }
    }
  }

  // (3) contribution
  {
    const double delta = -beta * (I4_ - 1.0);
    for (unsigned a = 0; a < 6; ++a)
    {
      for (unsigned b = 0; b < 6; ++b)
      {
        const unsigned i = MAT::IMap::fourth_[a][b][0];
        const unsigned j = MAT::IMap::fourth_[a][b][1];
        const unsigned k = MAT::IMap::fourth_[a][b][2];
        const unsigned l = MAT::IMap::fourth_[a][b][3];

        cmat(a, b) +=
            delta * (rcg_inv_s(MAT::IMap::second_[i][k]) * rcg_inv_s(MAT::IMap::second_[j][l]) +
                        rcg_inv_s(MAT::IMap::second_[i][l]) * rcg_inv_s(MAT::IMap::second_[j][k]));
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::ResetInvariants(const LINALG::Matrix<6, 1>& rcg)
{
  // calculate the square root of the third invariant alias the determinant
  // of the deformation gradient
  const double I3 = rcg(0) * rcg(1) * rcg(2) + 0.25 * rcg(3) * rcg(4) * rcg(5) -
                    0.25 * rcg(1) * rcg(5) * rcg(5) - 0.25 * rcg(2) * rcg(3) * rcg(3) -
                    0.25 * rcg(0) * rcg(4) * rcg(4);
  if (I3 < 0.0) dserror("I3 is negative!");

  // jacobian determinant
  j_ = std::sqrt(I3);

  // calculate pseudo invariant I4 ( strain measure in fiber direction )
  I4_ = AA_(0) * rcg(0) + AA_(1) * rcg(1) + AA_(2) * rcg(2) + AA_(3) * rcg(3) + AA_(4) * rcg(4) +
        AA_(5) * rcg(5);

  // calculate pseudo invariant I5 ( quad. strain measure in fiber direction )
  LINALG::Matrix<6, 1> rcg_quad(false);
  MAT::VStrainUtils::PowerOfSymmetricTensor(2, rcg, rcg_quad);
  I5_ = AA_(0) * (rcg_quad(0)) + AA_(1) * (rcg_quad(1)) + AA_(2) * (rcg_quad(2)) +
        AA_(3) * (rcg_quad(3)) + AA_(4) * (rcg_quad(4)) + AA_(5) * (rcg_quad(5));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void MAT::ELASTIC::CoupTransverselyIsotropic::UpdateSecondPiolaKirchhoffStress(
    LINALG::Matrix<6, 1>& stress, const LINALG::Matrix<6, 1>& rcg_s,
    LINALG::Matrix<6, 1>& rcg_inv_s) const
{
  const double alpha = params_->alpha_;
  const double beta = params_->beta_;
  const double gamma = params_->gamma_;

  // compute inverse right Cauchy Green tensor
  MAT::VStressUtils::InverseTensor(rcg_s, rcg_inv_s);

  // (0) contribution
  {
    const double fac = beta * (I4_ - 1.0);
    stress.Update(fac, rcg_inv_s, 1.0);
  }

  // (1) contribution
  {
    const double fac = 2.0 * (alpha + beta * std::log(j_) + 2.0 * gamma * (I4_ - 1.0));
    stress.Update(fac, AA_, 1.0);
  }

  // (2) contribution
  {
    LINALG::Matrix<3, 1> ca(true);
    MAT::VStressUtils::MultiplyTensorVector(rcg_s, A_, ca);

    LINALG::Matrix<6, 1> caa_aac(true);
    MAT::VStressUtils::SymmetricOuterProduct(ca, A_, caa_aac);

    const double fac = -alpha;
    stress.Update(fac, caa_aac, 1.0);
  }
}
