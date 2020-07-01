/*----------------------------------------------------------------------*/
/*! \file
\brief Implentation for an exponential strain energy function for fibers

\level 3

*/

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexposhear.H"
#include "elast_aniso_structuraltensor_strategy.H"

#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/voigt_notation.H"


MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::CoupAnisoExpoShearAnisotropyExtension(
    const int init_mode, const std::array<int, 2> fiber_ids)
    : init_mode_(init_mode), fiber_ids_(fiber_ids)
{
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::PackAnisotropy(
    DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, scalarProducts_);
  DRT::ParObject::AddtoPack(data, structuralTensors_stress_);
  DRT::ParObject::AddtoPack(data, structuralTensors_);
  DRT::ParObject::AddtoPack(data, isInitialized_);
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DRT::ParObject::ExtractfromPack(position, data, scalarProducts_);
  DRT::ParObject::ExtractfromPack(position, data, structuralTensors_stress_);
  DRT::ParObject::ExtractfromPack(position, data, structuralTensors_);
  isInitialized_ = static_cast<bool>(DRT::ParObject::ExtractInt(position, data));
}

double MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::GetScalarProduct(int gp) const
{
  if (!isInitialized_)
  {
    dserror("Fibers have not been initialized yet.");
  }

  if (scalarProducts_.size() == 1)
  {
    return scalarProducts_[0];
  }

  return scalarProducts_[gp];
}

const LINALG::Matrix<3, 3>&
MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::GetStructuralTensor(int gp) const
{
  if (!isInitialized_)
  {
    dserror("Fibers have not been initialized yet.");
  }

  if (structuralTensors_.size() == 1)
  {
    return structuralTensors_[0];
  }

  return structuralTensors_[gp];
}

const LINALG::Matrix<6, 1>&
MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::GetStructuralTensor_stress(int gp) const
{
  if (!isInitialized_)
  {
    dserror("Fibers have not been initialized yet.");
  }

  if (structuralTensors_stress_.size() == 1)
  {
    return structuralTensors_stress_[0];
  }

  return structuralTensors_stress_[gp];
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::OnGlobalDataInitialized()
{
  // do nothing
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::OnGlobalElementDataInitialized()
{
  if (init_mode_ == DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_EXTERNAL ||
      init_mode_ == DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_FIBERS)
  {
    // this is the initalization method for element fibers, so do nothing here
    return;
  }

  if (init_mode_ == DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_EXTERNAL)
  {
    dserror(
        "This material only supports the fiber prescription with the FIBER1 FIBER2 notation and "
        "INIT modes %d and %d.",
        DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_FIBERS,
        DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_FIBERS);
  }

  if (GetAnisotropy()->GetElementFibers().empty())
  {
    dserror("No element fibers are given with the FIBER1 FIBER2 notation");
  }

  scalarProducts_.resize(1);
  structuralTensors_.resize(1);
  structuralTensors_stress_.resize(1);
  scalarProducts_[0] = GetAnisotropy()
                           ->GetElementFiber(fiber_ids_[0])
                           .Dot(GetAnisotropy()->GetElementFiber(fiber_ids_[1]));

  LINALG::Matrix<3, 3> fiber1fiber2T(false);
  fiber1fiber2T.MultiplyNT(GetAnisotropy()->GetElementFiber(fiber_ids_[0]),
      GetAnisotropy()->GetElementFiber(fiber_ids_[1]));

  structuralTensors_[0].Update(0.5, fiber1fiber2T);
  structuralTensors_[0].UpdateT(0.5, fiber1fiber2T, 1.0);
  UTILS::VOIGT::Stresses::MatrixToVector(structuralTensors_[0], structuralTensors_stress_[0]);

  isInitialized_ = true;
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::OnGlobalGPDataInitialized()
{
  if (init_mode_ == DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_EXTERNAL ||
      init_mode_ == DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_FIBERS)
  {
    // this is the initalization method for nodal fibers, so do nothing here
    return;
  }

  if (init_mode_ == DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_EXTERNAL)
  {
    dserror(
        "This material only supports the fiber prescription with the FIBER1 FIBER2 notation and "
        "INIT modes %d and %d.",
        DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_FIBERS,
        DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_FIBERS);
  }

  if (GetAnisotropy()->GetNumberOfGPFibers() == 0)
  {
    dserror("No element fibers are given with the FIBER1 FIBER2 notation");
  }

  scalarProducts_.resize(GetAnisotropy()->GetNumberOfGaussPoints());
  structuralTensors_.resize(GetAnisotropy()->GetNumberOfGaussPoints());
  structuralTensors_stress_.resize(GetAnisotropy()->GetNumberOfGaussPoints());

  for (auto gp = 0; gp < GetAnisotropy()->GetNumberOfGaussPoints(); ++gp)
  {
    scalarProducts_[gp] = GetAnisotropy()
                              ->GetGPFiber(gp, fiber_ids_[0])
                              .Dot(GetAnisotropy()->GetGPFiber(gp, fiber_ids_[1]));

    LINALG::Matrix<3, 3> fiber1fiber2T(false);
    fiber1fiber2T.MultiplyNT(GetAnisotropy()->GetGPFiber(gp, fiber_ids_[0]),
        GetAnisotropy()->GetGPFiber(gp, fiber_ids_[1]));

    structuralTensors_[gp].Update(0.5, fiber1fiber2T);
    structuralTensors_[gp].UpdateT(0.5, fiber1fiber2T, 1.0);
    UTILS::VOIGT::Stresses::MatrixToVector(structuralTensors_[gp], structuralTensors_stress_[gp]);
  }

  isInitialized_ = true;
}

MAT::ELASTIC::PAR::CoupAnisoExpoShear::CoupAnisoExpoShear(Teuchos::RCP<MAT::PAR::Material> matdata)
    : MAT::PAR::Parameter(matdata), MAT::ELASTIC::PAR::CoupAnisoExpoBase(matdata)
{
  std::copy_n(matdata->Get<std::vector<int>>("FIBER_IDS")->begin(), 2, fiber_id_.begin());

  for (unsigned i = 0; i < fiber_id_.size(); ++i) fiber_id_[i] -= 1;
}

MAT::ELASTIC::CoupAnisoExpoShear::CoupAnisoExpoShear(MAT::ELASTIC::PAR::CoupAnisoExpoShear* params)
    : MAT::ELASTIC::CoupAnisoExpoBase(params),
      params_(params),
      anisotropyExtension_(params_->init_, params->fiber_id_)
{
}

void MAT::ELASTIC::CoupAnisoExpoShear::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropyExtension_);
}

void MAT::ELASTIC::CoupAnisoExpoShear::PackSummand(DRT::PackBuffer& data) const
{
  anisotropyExtension_.PackAnisotropy(data);
}

void MAT::ELASTIC::CoupAnisoExpoShear::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  anisotropyExtension_.UnpackAnisotropy(data, position);
}

void MAT::ELASTIC::CoupAnisoExpoShear::GetFiberVecs(std::vector<LINALG::Matrix<3, 1>>& fibervecs)
{
  // no fibers to export here
}

void MAT::ELASTIC::CoupAnisoExpoShear::SetFiberVecs(
    const double newgamma, const LINALG::Matrix<3, 3>& locsys, const LINALG::Matrix<3, 3>& defgrd)
{
  dserror("This function is not implemented for this summand!");
}
