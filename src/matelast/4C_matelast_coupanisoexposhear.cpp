/*----------------------------------------------------------------------*/
/*! \file
\brief Implentation for an exponential strain energy function for fibers

\level 3
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coupanisoexposhear.hpp"

#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_par_material.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"

FOUR_C_NAMESPACE_OPEN


MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::CoupAnisoExpoShearAnisotropyExtension(
    const int init_mode, const std::array<int, 2> fiber_ids)
    : init_mode_(init_mode), fiber_ids_(fiber_ids)
{
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::PackAnisotropy(
    CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::ParObject::AddtoPack(data, scalar_products_);
  CORE::COMM::ParObject::AddtoPack(data, structural_tensors_stress_);
  CORE::COMM::ParObject::AddtoPack(data, structural_tensors_);
  CORE::COMM::ParObject::AddtoPack(data, is_initialized_);
}

void MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  CORE::COMM::ParObject::ExtractfromPack(position, data, scalar_products_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, structural_tensors_stress_);
  CORE::COMM::ParObject::ExtractfromPack(position, data, structural_tensors_);
  is_initialized_ = static_cast<bool>(CORE::COMM::ParObject::ExtractInt(position, data));
}

double MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::GetScalarProduct(int gp) const
{
  if (!is_initialized_)
  {
    FOUR_C_THROW("Fibers have not been initialized yet.");
  }

  if (scalar_products_.size() == 1)
  {
    return scalar_products_[0];
  }

  return scalar_products_[gp];
}

const CORE::LINALG::Matrix<3, 3>&
MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::GetStructuralTensor(int gp) const
{
  if (!is_initialized_)
  {
    FOUR_C_THROW("Fibers have not been initialized yet.");
  }

  if (structural_tensors_.size() == 1)
  {
    return structural_tensors_[0];
  }

  return structural_tensors_[gp];
}

const CORE::LINALG::Matrix<6, 1>&
MAT::ELASTIC::CoupAnisoExpoShearAnisotropyExtension::GetStructuralTensor_stress(int gp) const
{
  if (!is_initialized_)
  {
    FOUR_C_THROW("Fibers have not been initialized yet.");
  }

  if (structural_tensors_stress_.size() == 1)
  {
    return structural_tensors_stress_[0];
  }

  return structural_tensors_stress_[gp];
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
    FOUR_C_THROW(
        "This material only supports the fiber prescription with the FIBER1 FIBER2 notation and "
        "INIT modes %d and %d.",
        DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_FIBERS,
        DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_FIBERS);
  }

  if (GetAnisotropy()->GetElementFibers().empty())
  {
    FOUR_C_THROW("No element fibers are given with the FIBER1 FIBER2 notation");
  }

  scalar_products_.resize(1);
  structural_tensors_.resize(1);
  structural_tensors_stress_.resize(1);
  scalar_products_[0] = GetAnisotropy()
                            ->GetElementFiber(fiber_ids_[0])
                            .Dot(GetAnisotropy()->GetElementFiber(fiber_ids_[1]));

  CORE::LINALG::Matrix<3, 3> fiber1fiber2T(false);
  fiber1fiber2T.MultiplyNT(GetAnisotropy()->GetElementFiber(fiber_ids_[0]),
      GetAnisotropy()->GetElementFiber(fiber_ids_[1]));

  structural_tensors_[0].Update(0.5, fiber1fiber2T);
  structural_tensors_[0].UpdateT(0.5, fiber1fiber2T, 1.0);
  CORE::LINALG::VOIGT::Stresses::MatrixToVector(
      structural_tensors_[0], structural_tensors_stress_[0]);

  is_initialized_ = true;
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
    FOUR_C_THROW(
        "This material only supports the fiber prescription with the FIBER1 FIBER2 notation and "
        "INIT modes %d and %d.",
        DefaultAnisotropyExtension<2>::INIT_MODE_ELEMENT_FIBERS,
        DefaultAnisotropyExtension<2>::INIT_MODE_NODAL_FIBERS);
  }

  if (GetAnisotropy()->GetNumberOfGPFibers() == 0)
  {
    FOUR_C_THROW("No element fibers are given with the FIBER1 FIBER2 notation");
  }

  scalar_products_.resize(GetAnisotropy()->GetNumberOfGaussPoints());
  structural_tensors_.resize(GetAnisotropy()->GetNumberOfGaussPoints());
  structural_tensors_stress_.resize(GetAnisotropy()->GetNumberOfGaussPoints());

  for (auto gp = 0; gp < GetAnisotropy()->GetNumberOfGaussPoints(); ++gp)
  {
    scalar_products_[gp] = GetAnisotropy()
                               ->GetGPFiber(gp, fiber_ids_[0])
                               .Dot(GetAnisotropy()->GetGPFiber(gp, fiber_ids_[1]));

    CORE::LINALG::Matrix<3, 3> fiber1fiber2T(false);
    fiber1fiber2T.MultiplyNT(GetAnisotropy()->GetGPFiber(gp, fiber_ids_[0]),
        GetAnisotropy()->GetGPFiber(gp, fiber_ids_[1]));

    structural_tensors_[gp].Update(0.5, fiber1fiber2T);
    structural_tensors_[gp].UpdateT(0.5, fiber1fiber2T, 1.0);
    CORE::LINALG::VOIGT::Stresses::MatrixToVector(
        structural_tensors_[gp], structural_tensors_stress_[gp]);
  }

  is_initialized_ = true;
}

MAT::ELASTIC::PAR::CoupAnisoExpoShear::CoupAnisoExpoShear(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MAT::PAR::Parameter(matdata), MAT::ELASTIC::PAR::CoupAnisoExpoBase(matdata)
{
  std::copy_n(matdata->Get<std::vector<int>>("FIBER_IDS").begin(), 2, fiber_id_.begin());

  for (int& i : fiber_id_) i -= 1;
}

MAT::ELASTIC::CoupAnisoExpoShear::CoupAnisoExpoShear(MAT::ELASTIC::PAR::CoupAnisoExpoShear* params)
    : MAT::ELASTIC::CoupAnisoExpoBase(params),
      params_(params),
      anisotropy_extension_(params_->init_, params->fiber_id_)
{
}

void MAT::ELASTIC::CoupAnisoExpoShear::RegisterAnisotropyExtensions(MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropy_extension_);
}

void MAT::ELASTIC::CoupAnisoExpoShear::PackSummand(CORE::COMM::PackBuffer& data) const
{
  anisotropy_extension_.PackAnisotropy(data);
}

void MAT::ELASTIC::CoupAnisoExpoShear::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  anisotropy_extension_.UnpackAnisotropy(data, position);
}

void MAT::ELASTIC::CoupAnisoExpoShear::GetFiberVecs(
    std::vector<CORE::LINALG::Matrix<3, 1>>& fibervecs)
{
  // no fibers to export here
}

void MAT::ELASTIC::CoupAnisoExpoShear::SetFiberVecs(const double newgamma,
    const CORE::LINALG::Matrix<3, 3>& locsys, const CORE::LINALG::Matrix<3, 3>& defgrd)
{
  FOUR_C_THROW("This function is not implemented for this summand!");
}

FOUR_C_NAMESPACE_CLOSE
