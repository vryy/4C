/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam3r_plasticity.H"
#include "beam_elasthyper_parameter.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../headers/FAD_utils.H"
#include <cmath>

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamReissnerElastPlasticMaterialParams::BeamReissnerElastPlasticMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamReissnerElastHyperMaterialParams(matdata),
      yield_stress_n_(matdata->GetDouble("YIELDN")),
      yield_stress_m_(matdata->GetDouble("YIELDM")),
      isohard_modulus_n_(matdata->GetDouble("ISOHARDN")),
      isohard_modulus_m_(matdata->GetDouble("ISOHARDM")),
      torsion_plasticity_((bool)matdata->GetDouble("TORSIONPLAST"))
{
  if (yield_stress_n_ == -1.0 && yield_stress_m_ == -1.0 && isohard_modulus_n_ == -1.0 &&
      isohard_modulus_m_ == -1.0)
    dserror("no plasticity material parameter is given; use elastic material instead");

  if (isohard_modulus_n_ <= 0.0) yield_stress_n_ = -1.0;

  if (isohard_modulus_m_ <= 0.0) yield_stress_m_ = -1.0;

  if (torsion_plasticity_ && std::abs(GetYoungsModulus() - 2.0 * GetShearModulus()) > 1e-9)
  {
    dserror(
        "Young's modulus must be equal to two times the shear modulus if plasticity for torsional "
        "moments is turned on");
  }

  if (yield_stress_m_ >= 0 && std::abs(GetMomentInertia2() - GetMomentInertia3()) > 1e-9)
    dserror("area moment of inertia 2 and 3 need to be equal");

  if (torsion_plasticity_ && std::abs(GetMomentInertiaPolar() - 2.0 * GetMomentInertia2()) > 1e-9)
  {
    dserror(
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
    dserror(
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
DRT::ParObject* MAT::BeamElastPlasticMaterialType<T>::Create(const std::vector<char>& data)
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
// Pack data
template <typename T>
void MAT::BeamPlasticMaterial<T>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
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
  int type = 0;
  this->ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  int matid;
  this->ExtractfromPack(position, data, matid);

  this->ExtractfromPack(position, data, gammaplastaccum_);
  this->ExtractfromPack(position, data, gammaplastconv_);
  this->ExtractfromPack(position, data, kappaplastaccum_);
  this->ExtractfromPack(position, data, kappaplastconv_);

  this->SetParameter(nullptr);

  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);


      this->SetParameter(static_cast<MAT::PAR::BeamReissnerElastPlasticMaterialParams*>(mat));
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::EvaluateForceContributionsToStress(
    LINALG::Matrix<3, 1, T>& stressN, const LINALG::Matrix<3, 3, T>& CN,
    const LINALG::Matrix<3, 1, T>& Gamma)
{
  //*************Begin: Plasticity of strains in axial direction

  // If no yielding parameter was given, the material is modeled in a purely elastic manner
  if (this->Params().GetYieldStressN() < 0.0)
  {
    // compute material stresses by multiplying strains with constitutive matrix
    MAT::BeamElastHyperMaterial<T>::EvaluateForceContributionsToStress(stressN, CN, Gamma);
  }
  else
  {
    // material elastic strain
    LINALG::Matrix<3, 1, T> Gammaelast(true);

    // compute elastic strain
    for (int i = 0; i < 3; i++)
    {
      Gammaelast(i) = Gamma(i) - gammaplastconv_(i);
    }
    // compute resulting stress
    stressN.Multiply(CN, Gammaelast);

    // check if yield stress is surpassed
    if (std::abs(stressN(0)) > effyieldstressN_)
    {
      // compute plastic strain increment
      deltastressN_(0) = std::abs(stressN(0)) - effyieldstressN_;

      deltaGammaplast_(0) = ((CN(0, 0) - cN_eff_(0, 0)) / CN(0, 0) * deltastressN_(0) / CN(0, 0)) *
                            FADUTILS::Signum(stressN(0));

      gammaplastnew_(0) = gammaplastconv_(0) + deltaGammaplast_(0);

      // update elastic strain and stress
      for (int i = 0; i < 3; i++)
      {
        Gammaelast(i) = Gamma(i) - gammaplastnew_(i);
      }

      stressN.Multiply(CN, Gammaelast);
    }
    stressN_ = std::abs(stressN(0));
  }
  //*************End: Plasticity of strains in axial direction
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::EvaluateMomentContributionsToStress(
    LINALG::Matrix<3, 1, T>& stressM, const LINALG::Matrix<3, 3, T>& CM,
    const LINALG::Matrix<3, 1, T>& Cur)
{
  //*************Begin: Plasticity of curvatures

  // If no yielding parameter was given, the material is modeled in a purely elastic manner
  if (this->Params().GetYieldStressM() < 0.0)
  {
    // compute material stresses by multiplying curvature with constitutive matrix
    MAT::BeamElastHyperMaterial<T>::EvaluateMomentContributionsToStress(stressM, CM, Cur);
  }
  else
  {
    // If torsional plasticity is turned on, use full curvature vector for plasticity,
    // else, continue with reduced curvature vector (first entry is zero)
    if (this->Params().GetTorsionPlasticity())
    {
      kappa_(0) = Cur(0);
    }
    kappa_(1) = Cur(1);
    kappa_(2) = Cur(2);


    // return-mapping algorithm

    // compute elastic curvature
    for (int i = 0; i < 3; i++)
    {
      kappaelast_(i) = kappa_(i) - kappaplastconv_(i);
    }

    // compute resulting moments
    stressM.Multiply(CM, kappaelast_);

    // compute norm of moment vector
    normstressM_ = stressM.Norm2();

    // compute fraction that exceeds the current yield moment
    deltastressM_ = normstressM_ - effyieldstressM_;

    // check if yield moment is surpassed
    if (deltastressM_ > 0.0)
    {
      // compute plastic curvature increment
      deltaKappaplast_ = (CM(1, 1) - cM_eff_(1, 1)) / CM(1, 1) * deltastressM_ / CM(1, 1);

      // update plastic curvature
      for (int i = 0; i < 3; i++)
      {
        kappaplastnew_(i) =
            kappaplastconv_(i) + deltaKappaplast_ * kappaelast_(i) / kappaelast_.Norm2();
      }

      // update elastic curvature
      for (int i = 0; i < 3; i++)
      {
        kappaelast_(i) = kappa_(i) - kappaplastnew_(i);
      }

      // update moment vector and its norm
      stressM.Multiply(CM, kappaelast_);
      normstressM_ = stressM.Norm2();
    }

    // if torsional plasticity is turned off, the moment needs to be recomputed using the full
    // elastic curvature (kappaelast(0) is zero in this case)
    if (!this->Params().GetTorsionPlasticity())
    {
      kappaelast_(0) = Cur(0);
      stressM.Multiply(CM, kappaelast_);
    }
  }

  //*************End: Plasticity of curvatures
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::ComputeConstitutiveParameter(
    LINALG::Matrix<3, 3, T>& C_N, LINALG::Matrix<3, 3, T>& C_M)
{
  MAT::BeamElastHyperMaterial<T>::ComputeConstitutiveParameter(C_N, C_M);

  // If plasticity for axial strains is enabled, get hardening constitutive parameters
  if (this->Params().GetYieldStressN() >= 0)
  {
    GetHardeningConstitutiveMatrixOfForcesMaterialFrame(cN_eff_);
    gammaplastnew_ = gammaplastconv_;
    GetEffectiveYieldStressN(
        effyieldstressN_, this->Params().GetYieldStressN(), C_N(0, 0), cN_eff_(0, 0));
  }

  // If plasticity for curvatures is enabled, get hardening constitutive parameters
  if (this->Params().GetYieldStressM() >= 0)
  {
    GetHardeningConstitutiveMatrixOfMomentsMaterialFrame(cM_eff_);
    kappaplastnew_ = kappaplastconv_;
    GetEffectiveYieldStressM(
        effyieldstressM_, this->Params().GetYieldStressM(), C_M(1, 1), cM_eff_(1, 1));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

template <typename T>
void MAT::BeamPlasticMaterial<T>::Update()
{
  gammaplastaccum_ += std::abs(gammaplastconv_(0) - gammaplastnew_(0));
  gammaplastconv_ = gammaplastnew_;
  kappaplastaccum_ += std::sqrt(
      (kappaplastconv_(0) - kappaplastnew_(0)) * (kappaplastconv_(0) - kappaplastnew_(0)) +
      (kappaplastconv_(1) - kappaplastnew_(1)) * (kappaplastconv_(1) - kappaplastnew_(1)) +
      (kappaplastconv_(2) - kappaplastnew_(2)) * (kappaplastconv_(2) - kappaplastnew_(2)));
  kappaplastconv_ = kappaplastnew_;
  cN_eff_.Scale(0.0);
  cM_eff_.Scale(0.0);
  kappa_.Scale(0.0);
  kappaelast_.Scale(0.0);
  kappaelastflow_.Scale(0.0);
  elastic_curvature_.Scale(0.0);
  deltaGammaplast_.Scale(0.0);
  deltastressN_.Scale(0.0);
  deltaKappaplast_ = 0.0;
  normstressM_ = 0.0;
  deltastressM_ = 0.0;
  effyieldstressN_ = 0.0;
  effyieldstressM_ = 0.0;
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
    LINALG::Matrix<3, 3, T>& C_N) const
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
    LINALG::Matrix<3, 3, T>& C_M) const
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
    LINALG::Matrix<3, 3, T>& CN_eff) const
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
    LINALG::Matrix<3, 3, T>& CM_eff) const
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
    T& eff_yieldN, T init_yieldN, T CN_0, T CN_eff_0) const
{
  eff_yieldN = init_yieldN + (CN_0 * CN_eff_0) / (CN_0 - CN_eff_0) * gammaplastaccum_;
}


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetEffectiveYieldStressM(
    T& eff_yieldM, T init_yieldM, T CM_1, T CM_eff_1) const
{
  eff_yieldM = init_yieldM + (CM_1 * CM_eff_1) / (CM_1 - CM_eff_1) * kappaplastaccum_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamPlasticMaterial<T>::GetStiffnessMatrixOfMoments(
    LINALG::Matrix<3, 3, T>& stiffM, const LINALG::Matrix<3, 3, T>& C_M)
{
  /* compute spatial stresses and constitutive matrix from material ones according to Jelenic
   * 1999, page 148, paragraph between (2.22) and (2.23) and Romero 2004, (3.10)*/
  if (this->Params().GetYieldStressM() < 0 || normstressM_ + 10e-10 < effyieldstressM_)
    MAT::BeamElastHyperMaterial<T>::GetStiffnessMatrixOfMoments(stiffM, C_M);
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
      kappaelastflow_(i) = kappaelast_(i);
    }
    normKappaelastflow = kappaelastflow_.Norm2();

    // compute e, which is the unit vector in the direction of kappaelastflow
    for (int i = i_start; i < 3; i++)
    {
      elastic_curvature_(i) = kappaelastflow_(i) / normKappaelastflow;
    }

    // compute the stiffness matrix of moments
    for (int i = i_start; i < 3; i++)
    {
      for (int j = i_start; j < 3; j++)
      {
        stiffM(i, j) =
            C_M(1, 1) *
            ((deltaKappaplast_ / normKappaelastflow - (1.0 - cM_eff_(1, 1) / C_M(1, 1))) *
                elastic_curvature_(i) * elastic_curvature_(j));

        if (i == j)
          stiffM(i, j) = stiffM(i, j) + C_M(1, 1) * (1.0 - deltaKappaplast_ / normKappaelastflow);
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
    LINALG::Matrix<3, 3, T>& stiffN, const LINALG::Matrix<3, 3, T>& C_N)
{
  if (this->Params().GetYieldStressN() < 0.0 || stressN_ < effyieldstressN_)
  {
    stiffN = C_N;
  }
  else
  {
    stiffN = cN_eff_;
  }
}


// explicit template instantiations
template class MAT::BeamPlasticMaterial<double>;

template class MAT::BeamElastPlasticMaterialType<double>;
