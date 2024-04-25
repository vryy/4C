/*----------------------------------------------------------------------*/
/*! \file

\brief

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mat_anisotropy_extension.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_mat_anisotropy_utils.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

template <unsigned int numfib>
MAT::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension(
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy)
    : fibers_(0),
      fiber_structural_tensors_stress_(0),
      fiber_structural_tensors_(0),
      structural_tensor_strategy_(stucturalTensorStrategy)
{
}

template <unsigned int numfib>
MAT::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension()
    : fibers_(0),
      fiber_structural_tensors_stress_(0),
      fiber_structural_tensors_(0),
      structural_tensor_strategy_(Teuchos::null)
{
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::ComputeStructuralTensors_stress()
{
  MAT::ComputeStructuralTensors<CORE::LINALG::Matrix<6, 1>, numfib>(
      fibers_, fiber_structural_tensors_stress_, structural_tensor_strategy_);
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::ComputeStructuralTensors()
{
  MAT::ComputeStructuralTensors<CORE::LINALG::Matrix<3, 3>, numfib>(
      fibers_, fiber_structural_tensors_, structural_tensor_strategy_);
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::SetFibers(
    int gp, const std::array<CORE::LINALG::Matrix<3, 1>, numfib>& fibers)
{
  if (gp >= GetAnisotropy()->GetNumberOfGaussPoints())
  {
    FOUR_C_THROW(
        "The current Gauss point %i is out of range of the expected number of Gauss points %i.", gp,
        GetAnisotropy()->GetNumberOfGaussPoints());
  }

  if (fibers_.empty())
  {
    fibers_.resize(GetFibersPerElement());
  }

  // Store fibers
  for (std::vector<CORE::LINALG::Matrix<3, 1>>::size_type i = 0; i < fibers.size(); ++i)
  {
    fibers_[gp].at(i).Update(fibers.at(i));
  }

  // Compute needed structural tensors
  ComputeNeededStructuralTensors();

  // Call the notifier method
  OnFibersInitialized();
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::SetFibers(
    const std::vector<std::array<CORE::LINALG::Matrix<3, 1>, numfib>>& fibers)
{
  fibers_ = fibers;

  // Compute needed structural tensors
  ComputeNeededStructuralTensors();

  // Call the notifier method
  OnFibersInitialized();
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::ComputeNeededStructuralTensors()
{
  if (static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR))
  {
    ComputeStructuralTensors();
  }
  if (static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR_STRESS))
  {
    ComputeStructuralTensors_stress();
  }
}

template <unsigned int numfib>
const CORE::LINALG::Matrix<3, 1>& MAT::FiberAnisotropyExtension<numfib>::GetFiber(
    int gp, int i) const
{
  switch (fiber_location_)
  {
    case FiberLocation::ElementFibers:
      return fibers_[GPDEFAULT][i];
    case FiberLocation::GPFibers:
      return fibers_[gp][i];
    default:
      FOUR_C_THROW(
          "You have not specified, whether you want fibers on GP level or on element level.");
  }

  // just for compilation reasons. We will never land here because of the FOUR_C_THROW() above
  std::abort();
}

template <unsigned int numfib>
const CORE::LINALG::Matrix<6, 1>& MAT::FiberAnisotropyExtension<numfib>::GetStructuralTensor_stress(
    int gp, int i) const
{
  if (not static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR_STRESS))
  {
    FOUR_C_THROW("You have not specified that you need the fiber vector.");
  }
  switch (fiber_location_)
  {
    case FiberLocation::ElementFibers:
      return fiber_structural_tensors_stress_[GPDEFAULT][i];
    case FiberLocation::GPFibers:
      return fiber_structural_tensors_stress_[gp][i];
    default:
      FOUR_C_THROW(
          "You have not specified, whether you want fibers on GP level or on element level.");
  }

  // just for compilation reasons. We will never land here because of the FOUR_C_THROW() above
  std::abort();
}

template <unsigned int numfib>
const CORE::LINALG::Matrix<3, 3>& MAT::FiberAnisotropyExtension<numfib>::GetStructuralTensor(
    int gp, int i) const
{
  if (not static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR))
  {
    FOUR_C_THROW("You have not specified that you need the structural tensor.");
  }
  switch (fiber_location_)
  {
    case FiberLocation::ElementFibers:
      return fiber_structural_tensors_[GPDEFAULT][i];
    case FiberLocation::GPFibers:
      return fiber_structural_tensors_[gp][i];
    default:
      FOUR_C_THROW(
          "You have not specified, whether you want fibers on GP level or on element level.");
  }

  // just for compilation reasons. We will never land here because of the FOUR_C_THROW() above
  std::abort();
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::PackAnisotropy(CORE::COMM::PackBuffer& data) const
{
  PackFiberArray<CORE::LINALG::Matrix<3, 1>, numfib>(data, fibers_);
  PackFiberArray<CORE::LINALG::Matrix<6, 1>, numfib>(data, fiber_structural_tensors_stress_);
  PackFiberArray<CORE::LINALG::Matrix<3, 3>, numfib>(data, fiber_structural_tensors_);
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(tensor_flags_));
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(fiber_location_));
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  UnpackFiberArray<CORE::LINALG::Matrix<3, 1>, numfib>(position, data, fibers_);
  UnpackFiberArray<CORE::LINALG::Matrix<6, 1>, numfib>(
      position, data, fiber_structural_tensors_stress_);
  UnpackFiberArray<CORE::LINALG::Matrix<3, 3>, numfib>(position, data, fiber_structural_tensors_);
  tensor_flags_ = static_cast<uint_fast8_t>(CORE::COMM::ParObject::ExtractInt(position, data));
  fiber_location_ = static_cast<FiberLocation>(CORE::COMM::ParObject::ExtractInt(position, data));
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::SetFiberLocation(FiberLocation location)
{
  fiber_location_ = location;
}

template <unsigned int numfib>
int MAT::FiberAnisotropyExtension<numfib>::GetVirtualGaussPoint(const int gp) const
{
  if (fiber_location_ == FiberLocation::ElementFibers)
  {
    return GPDEFAULT;
  }

  return gp;
}

template <unsigned int numfib>
int MAT::FiberAnisotropyExtension<numfib>::GetFibersPerElement() const
{
  if (fiber_location_ == FiberLocation::ElementFibers)
  {
    return 1;
  }

  return GetAnisotropy()->GetNumberOfGaussPoints();
}

// explicit instatiations of template classes
template class MAT::FiberAnisotropyExtension<1u>;
template class MAT::FiberAnisotropyExtension<2u>;
FOUR_C_NAMESPACE_CLOSE
