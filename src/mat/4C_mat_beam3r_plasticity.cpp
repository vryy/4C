/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_mat_beam3r_plasticity.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_beam_elasthyper_parameter.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_fad.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamReissnerElastPlasticMaterialParams::BeamReissnerElastPlasticMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamReissnerElastHyperMaterialParams(matdata),
      yield_stress_n_(*matdata->Get<double>("YIELDN")),
      yield_stress_m_(*matdata->Get<double>("YIELDM")),
      isohard_modulus_n_(*matdata->Get<double>("ISOHARDN")),
      isohard_modulus_m_(*matdata->Get<double>("ISOHARDM")),
      torsion_plasticity_(*matdata->Get<bool>("TORSIONPLAST"))
{
  if (yield_stress_n_ == -1.0 && yield_stress_m_ == -1.0 && isohard_modulus_n_ == -1.0 &&
      isohard_modulus_m_ == -1.0)
    FOUR_C_THROW("no plasticity material parameter is given; use elastic material instead");

  if (isohard_modulus_n_ <= 0.0) yield_stress_n_ = -1.0;

  if (isohard_modulus_m_ <= 0.0) yield_stress_m_ = -1.0;

  if (torsion_plasticity_ && std::abs(GetYoungsModulus() - 2.0 * GetShearModulus()) > 1e-9)
  {
    FOUR_C_THROW(
        "Young's modulus must be equal to two times the shear modulus if plasticity for torsional "
        "moments is turned on");
  }

  if (yield_stress_m_ >= 0 && std::abs(GetMomentInertia2() - GetMomentInertia3()) > 1e-9)
    FOUR_C_THROW("area moment of inertia 2 and 3 need to be equal");

  if (torsion_plasticity_ && std::abs(GetMomentInertiaPolar() - 2.0 * GetMomentInertia2()) > 1e-9)
  {
    FOUR_C_THROW(
        "polar area moment of inertia needs to be assigned twice the value of area moment of "
        "inertia 2");
  }
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::BeamReissnerElastPlasticMaterialParams::CreateMaterial()
{
  /* all the different parameter sets (Reissner/Kirchhoff/..., 'classic'/'by modes') are used to
   * parameterize the same constitutive relations based on a hyperelastic stored energy function
   * formulated for cross-section resultants which are implemented in BeamElastHyperMaterial */
  Teuchos::RCP<MAT::Material> matobject;

  if (Uses_FAD())
  {
    FOUR_C_THROW(
        "The elastoplastic beam material is not yet implemented to be used with automatic "
        "differentiation!");
  }
  else
    matobject = Teuchos::rcp(new MAT::BeamPlasticMaterial<double>(this));
  return matobject;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
CORE::COMM::ParObject* MAT::BeamElastPlasticMaterialType<T>::Create(const std::vector<char>& data)
{
  // create material from packed data
  MAT::BeamPlasticMaterial<T>* matobject = new MAT::BeamPlasticMaterial<T>();
  matobject->Unpack(data);
  return matobject;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
MAT::BeamElastPlasticMaterialType<T> MAT::BeamElastPlasticMaterialType<T>::instance_;

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
MAT::BeamPlasticMaterial<T>::BeamPlasticMaterial(
    MAT::PAR::BeamReissnerElastPlasticMaterialParams* params)
    : MAT::BeamElastHyperMaterial<T>(params)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

template <typename T>
void MAT::BeamPlasticMaterial<T>::Setup(int numgp_force, int numgp_moment)
{
  c_n_eff_.resize(numgp_force);
  c_m_eff_.resize(numgp_moment);
  gammaplastconv_.resize(numgp_force);
  gammaplastnew_.resize(numgp_force);
  gammaplastaccum_.resize(numgp_force);
  kappaplastconv_.resize(numgp_moment);
  kappaplastnew_.resize(numgp_moment);
  kappaplastaccum_.resize(numgp_moment);
  effyieldstress_n_.resize(numgp_force);
  effyieldstress_m_.resize(numgp_moment);
  delta_kappaplast_.resize(numgp_moment);
  normstress_m_.resize(numgp_moment);
  deltastress_m_.resize(numgp_moment);
  kappaelast_.resize(numgp_moment);
  kappaelastflow_.resize(numgp_moment);
  elastic_curvature_.resize(numgp_moment);
  delta_gammaplast_.resize(numgp_force);
  deltastress_n_.resize(numgp_force);
  stress_n_.resize(numgp_force);

  for (int gp = 0; gp < numgp_force; gp++)
  {
    c_n_eff_[gp] = CORE::LINALG::Matrix<3, 3, T>(true);
    gammaplastconv_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    gammaplastnew_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    gammaplastaccum_[gp] = 0;
    effyieldstress_n_[gp] = 0;
    delta_kappaplast_[gp] = 0;
    delta_gammaplast_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    deltastress_n_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    stress_n_[gp] = 0;
  }

  for (int gp = 0; gp < numgp_moment; gp++)
  {
    c_m_eff_[gp] = CORE::LINALG::Matrix<3, 3, T>(true);
    kappaplastconv_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    kappaplastnew_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    kappaplastaccum_[gp] = 0;
    effyieldstress_m_[gp] = 0;
    normstress_m_[gp] = 0;
    deltastress_m_[gp] = 0;
    kappaelast_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    kappaelastflow_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
    elastic_curvature_[gp] = CORE::LINALG::Matrix<3, 1, T>(true);
  }

  numgp_force_ = numgp_force;
  numgp_moment_ = numgp_moment;
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
// Pack data
template <typename T>
void MAT::BeamPlasticMaterial<T>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  this->AddtoPack(data, type);

  // Pack material id
  int matid = -1;
  if (this->Parameter() != nullptr)
    matid = this->Params().Id();  // in case we are in post-process mode

  this->AddtoPack(data, matid);

  // Pack all internal variables
  this->AddtoPack(data, numgp_force_);
  this->AddtoPack(data, numgp_moment_);
  this->AddtoPack(data, gammaplastaccum_);
  this->AddtoPack(data, gammaplastconv_);
  this->AddtoPack(data, kappaplastaccum_);
  this->AddtoPack(data, kappaplastconv_);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
// Unpack data
template <typename T>
void MAT::BeamPlasticMaterial<T>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  int matid;
  this->ExtractfromPack(position, data, matid);

  this->ExtractfromPack(position, data, numgp_force_);
  this->ExtractfromPack(position, data, numgp_moment_);
  this->Setup(numgp_force_, numgp_moment_);
  this->ExtractfromPack(position, data, gammaplastaccum_);
  this->ExtractfromPack(position, data, gammaplastconv_);
  this->ExtractfromPack(position, data, kappaplastaccum_);
  this->ExtractfromPack(position, data, kappaplastconv_);

  this->SetParameter(nullptr);

  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);


      this->SetParameter(static_cast<MAT::PAR::BeamReissnerElastPlasticMaterialParams*>(mat));
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::EvaluateForceContributionsToStress(
    CORE::LINALG::Matrix<3, 1, T>& stressN, const CORE::LINALG::Matrix<3, 3, T>& CN,
    const CORE::LINALG::Matrix<3, 1, T>& Gamma, const unsigned int gp)
{
  //*************Begin: Plasticity of strains in axial direction

  // If no yielding parameter was given, the material is modeled in a purely elastic manner
  if (this->Params().GetYieldStressN() < 0.0)
  {
    // compute material stresses by multiplying strains with constitutive matrix
    MAT::BeamElastHyperMaterial<T>::EvaluateForceContributionsToStress(stressN, CN, Gamma, gp);
  }
  else
  {
    // material elastic strain
    CORE::LINALG::Matrix<3, 1, T> Gammaelast(true);

    // compute elastic strain
    for (int i = 0; i < 3; i++)
    {
      Gammaelast(i) = Gamma(i) - gammaplastconv_[gp](i);
    }
    // compute resulting stress
    stressN.Multiply(CN, Gammaelast);

    // check if yield stress is surpassed
    if (std::abs(stressN(0)) > effyieldstress_n_[gp])
    {
      // compute plastic strain increment
      deltastress_n_[gp](0) = std::abs(stressN(0)) - effyieldstress_n_[gp];

      delta_gammaplast_[gp](0) =
          ((CN(0, 0) - c_n_eff_[gp](0, 0)) / CN(0, 0) * deltastress_n_[gp](0) / CN(0, 0)) *
          CORE::FADUTILS::Signum(stressN(0));

      gammaplastnew_[gp](0) = gammaplastconv_[gp](0) + delta_gammaplast_[gp](0);

      // update elastic strain and stress
      for (int i = 0; i < 3; i++)
      {
        Gammaelast(i) = Gamma(i) - gammaplastnew_[gp](i);
      }

      stressN.Multiply(CN, Gammaelast);
    }
    stress_n_[gp] = std::abs(stressN(0));
  }
  //*************End: Plasticity of strains in axial direction
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::EvaluateMomentContributionsToStress(
    CORE::LINALG::Matrix<3, 1, T>& stressM, const CORE::LINALG::Matrix<3, 3, T>& CM,
    const CORE::LINALG::Matrix<3, 1, T>& Cur, const unsigned int gp)
{
  //*************Begin: Plasticity of curvatures

  // If no yielding parameter was given, the material is modeled in a purely elastic manner
  if (this->Params().GetYieldStressM() < 0.0)
  {
    // compute material stresses by multiplying curvature with constitutive matrix
    MAT::BeamElastHyperMaterial<T>::EvaluateMomentContributionsToStress(stressM, CM, Cur, gp);
  }
  else
  {
    //! copy of material curvature K (but first entry is 0 if torsional plasticity is turned off)
    CORE::LINALG::Matrix<3, 1, T> kappa{true};

    // If torsional plasticity is turned on, use full curvature vector for plasticity,
    // else, continue with reduced curvature vector (first entry is zero)
    if (this->Params().GetTorsionPlasticity())
    {
      kappa(0) = Cur(0);
    }
    kappa(1) = Cur(1);
    kappa(2) = Cur(2);


    // return-mapping algorithm

    // compute elastic curvature
    for (int i = 0; i < 3; i++)
    {
      kappaelast_[gp](i) = kappa(i) - kappaplastconv_[gp](i);
    }

    // compute resulting moments
    stressM.Multiply(CM, kappaelast_[gp]);

    // compute norm of moment vector
    normstress_m_[gp] = stressM.Norm2();

    // compute fraction that exceeds the current yield moment
    deltastress_m_[gp] = normstress_m_[gp] - effyieldstress_m_[gp];

    // check if yield moment is surpassed
    if (deltastress_m_[gp] > 0.0)
    {
      // compute plastic curvature increment
      delta_kappaplast_[gp] =
          (CM(1, 1) - c_m_eff_[gp](1, 1)) / CM(1, 1) * deltastress_m_[gp] / CM(1, 1);

      // update plastic curvature
      for (int i = 0; i < 3; i++)
      {
        kappaplastnew_[gp](i) = kappaplastconv_[gp](i) + delta_kappaplast_[gp] *
                                                             kappaelast_[gp](i) /
                                                             kappaelast_[gp].Norm2();
      }

      // update elastic curvature
      for (int i = 0; i < 3; i++)
      {
        kappaelast_[gp](i) = kappa(i) - kappaplastnew_[gp](i);
      }

      // update moment vector and its norm
      stressM.Multiply(CM, kappaelast_[gp]);
      normstress_m_[gp] = stressM.Norm2();
    }

    // if torsional plasticity is turned off, the moment needs to be recomputed using the full
    // elastic curvature (kappaelast(0) is zero in this case)
    if (!this->Params().GetTorsionPlasticity())
    {
      kappaelast_[gp](0) = Cur(0);
      stressM.Multiply(CM, kappaelast_[gp]);
    }
  }

  //*************End: Plasticity of curvatures
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::ComputeConstitutiveParameter(
    CORE::LINALG::Matrix<3, 3, T>& C_N, CORE::LINALG::Matrix<3, 3, T>& C_M)
{
  MAT::BeamElastHyperMaterial<T>::ComputeConstitutiveParameter(C_N, C_M);

  for (unsigned int gp = 0; gp < numgp_force_; gp++)
  {
    // If plasticity for axial strains is enabled, get hardening constitutive parameters
    if (this->Params().GetYieldStressN() >= 0)
    {
      GetHardeningConstitutiveMatrixOfForcesMaterialFrame(c_n_eff_[gp]);
      GetEffectiveYieldStressN(effyieldstress_n_[gp], this->Params().GetYieldStressN(), C_N(0, 0),
          c_n_eff_[gp](0, 0), gp);
    }
  }
  for (unsigned int gp = 0; gp < numgp_moment_; gp++)
  {
    // If plasticity for curvatures is enabled, get hardening constitutive parameters
    if (this->Params().GetYieldStressM() >= 0)
    {
      GetHardeningConstitutiveMatrixOfMomentsMaterialFrame(c_m_eff_[gp]);
      GetEffectiveYieldStressM(effyieldstress_m_[gp], this->Params().GetYieldStressM(), C_M(1, 1),
          c_m_eff_[gp](1, 1), gp);
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

template <typename T>
void MAT::BeamPlasticMaterial<T>::Update()
{
  for (unsigned int gp = 0; gp < numgp_force_; gp++)
  {
    gammaplastaccum_[gp] += std::abs(gammaplastconv_[gp](0) - gammaplastnew_[gp](0));
    gammaplastconv_[gp] = gammaplastnew_[gp];
    c_n_eff_[gp].putScalar(0.0);
    effyieldstress_n_[gp] = 0.0;
    delta_gammaplast_[gp].putScalar(0.0);
    deltastress_n_[gp].putScalar(0.0);
  }
  for (unsigned int gp = 0; gp < numgp_moment_; gp++)
  {
    kappaplastaccum_[gp] += std::sqrt((kappaplastconv_[gp](0) - kappaplastnew_[gp](0)) *
                                          (kappaplastconv_[gp](0) - kappaplastnew_[gp](0)) +
                                      (kappaplastconv_[gp](1) - kappaplastnew_[gp](1)) *
                                          (kappaplastconv_[gp](1) - kappaplastnew_[gp](1)) +
                                      (kappaplastconv_[gp](2) - kappaplastnew_[gp](2)) *
                                          (kappaplastconv_[gp](2) - kappaplastnew_[gp](2)));
    kappaplastconv_[gp] = kappaplastnew_[gp];
    c_m_eff_[gp].putScalar(0.0);
    effyieldstress_m_[gp] = 0.0;
    kappaelast_[gp].putScalar(0.0);
    kappaelastflow_[gp].putScalar(0.0);
    elastic_curvature_[gp].putScalar(0.0);
    delta_kappaplast_[gp] = 0.0;
    normstress_m_[gp] = 0.0;
    deltastress_m_[gp] = 0.0;
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::Reset()
{
  gammaplastnew_ = gammaplastconv_;
  kappaplastnew_ = kappaplastconv_;
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetConstitutiveMatrixOfForcesMaterialFrame(
    CORE::LINALG::Matrix<3, 3, T>& C_N) const
{
  // defining material constitutive matrix CN between Gamma and N
  // according to Jelenic 1999, section 2.4
  C_N.Clear();

  C_N(0, 0) = this->Params().GetAxialRigidity();
  C_N(1, 1) = this->Params().GetShearRigidity2();
  C_N(2, 2) = this->Params().GetShearRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetConstitutiveMatrixOfMomentsMaterialFrame(
    CORE::LINALG::Matrix<3, 3, T>& C_M) const
{
  // defining material constitutive matrix CM between curvature and moment
  // according to Jelenic 1999, section 2.4
  C_M.Clear();

  C_M(0, 0) = this->Params().GetTorsionalRigidity();
  C_M(1, 1) = this->Params().GetBendingRigidity2();
  C_M(2, 2) = this->Params().GetBendingRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
double MAT::BeamPlasticMaterial<T>::GetTranslationalMassInertiaFactor() const
{
  return this->Params().GetTranslationalMassInertia();
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
double MAT::BeamPlasticMaterial<T>::GetInteractionRadius() const
{
  return this->Params().GetInteractionRadius();
}


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetHardeningConstitutiveMatrixOfForcesMaterialFrame(
    CORE::LINALG::Matrix<3, 3, T>& CN_eff) const
{
  CN_eff.Clear();

  CN_eff(0, 0) = this->Params().GetHardeningAxialRigidity();
  CN_eff(1, 1) = this->Params().GetHardeningShearRigidity2();
  CN_eff(2, 2) = this->Params().GetHardeningShearRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetHardeningConstitutiveMatrixOfMomentsMaterialFrame(
    CORE::LINALG::Matrix<3, 3, T>& CM_eff) const
{
  CM_eff.Clear();
  if (this->Params().GetTorsionPlasticity())
    CM_eff(0, 0) = this->Params().GetHardeningMomentalRigidity();
  else
    CM_eff(0, 0) = this->Params().GetTorsionalRigidity();

  CM_eff(1, 1) = this->Params().GetHardeningMomentalRigidity();
  CM_eff(2, 2) = CM_eff(1, 1);
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetEffectiveYieldStressN(
    T& eff_yieldN, T init_yieldN, T CN_0, T CN_eff_0, const unsigned int gp) const
{
  eff_yieldN = init_yieldN + (CN_0 * CN_eff_0) / (CN_0 - CN_eff_0) * gammaplastaccum_[gp];
}


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetEffectiveYieldStressM(
    T& eff_yieldM, T init_yieldM, T CM_1, T CM_eff_1, const unsigned int gp) const
{
  eff_yieldM = init_yieldM + (CM_1 * CM_eff_1) / (CM_1 - CM_eff_1) * kappaplastaccum_[gp];
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetStiffnessMatrixOfMoments(
    CORE::LINALG::Matrix<3, 3, T>& stiffM, const CORE::LINALG::Matrix<3, 3, T>& C_M, const int gp)
{
  /* compute spatial stresses and constitutive matrix from material ones according to Jelenic
   * 1999, page 148, paragraph between (2.22) and (2.23) and Romero 2004, (3.10)*/

  if (this->Params().GetYieldStressM() < 0 || normstress_m_[gp] + 10e-10 < effyieldstress_m_[gp])
    MAT::BeamElastHyperMaterial<T>::GetStiffnessMatrixOfMoments(stiffM, C_M, gp);
  else
  {
    // Compute stiffness matrix for plastic regime:

    // norm of kappaelastflow
    T normKappaelastflow;

    // starting index for vector assignments
    int i_start = 0;

    // Set starting index s.t. the first entry of kappaelastflow and e will be 0 if
    // torsion plasticity is turned off
    if (!(this->Params().GetTorsionPlasticity()))
    {
      i_start = 1;
    }

    // initialize kappaelastflow, which points in the direction of the plastic increment
    for (int i = i_start; i < 3; i++)
    {
      kappaelastflow_[gp](i) = kappaelast_[gp](i);
    }
    normKappaelastflow = kappaelastflow_[gp].Norm2();

    // compute e, which is the unit vector in the direction of kappaelastflow
    for (int i = i_start; i < 3; i++)
    {
      elastic_curvature_[gp](i) = kappaelastflow_[gp](i) / normKappaelastflow;
    }

    // compute the stiffness matrix of moments
    for (int i = i_start; i < 3; i++)
    {
      for (int j = i_start; j < 3; j++)
      {
        stiffM(i, j) =
            C_M(1, 1) *
            ((delta_kappaplast_[gp] / normKappaelastflow - (1.0 - c_m_eff_[gp](1, 1) / C_M(1, 1))) *
                elastic_curvature_[gp](i) * elastic_curvature_[gp](j));

        if (i == j)
          stiffM(i, j) =
              stiffM(i, j) + C_M(1, 1) * (1.0 - delta_kappaplast_[gp] / normKappaelastflow);
      }
    }
    // if torsional plasticity is turned off, the first entry of the stiffness matrix is that of
    // cM
    if (i_start == 1)
    {
      stiffM(0, 0) = C_M(1, 1);
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetStiffnessMatrixOfForces(
    CORE::LINALG::Matrix<3, 3, T>& stiffN, const CORE::LINALG::Matrix<3, 3, T>& C_N, const int gp)
{
  if (this->Params().GetYieldStressN() < 0.0 || stress_n_[gp] < effyieldstress_n_[gp])
  {
    stiffN = C_N;
  }
  else
  {
    stiffN = c_n_eff_[gp];
  }
}


// explicit template instantiations
template class MAT::BeamPlasticMaterial<double>;

template class MAT::BeamElastPlasticMaterialType<double>;

FOUR_C_NAMESPACE_CLOSE
