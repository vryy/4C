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
Mat::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension(
    const Teuchos::RCP<Elastic::StructuralTensorStrategyBase>& stucturalTensorStrategy)
    : fibers_(0),
      fiber_structural_tensors_stress_(0),
      fiber_structural_tensors_(0),
      structural_tensor_strategy_(stucturalTensorStrategy)
{
}

template <unsigned int numfib>
Mat::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension()
    : fibers_(0),
      fiber_structural_tensors_stress_(0),
      fiber_structural_tensors_(0),
      structural_tensor_strategy_(Teuchos::null)
{
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::compute_structural_tensors_stress()
{
  Mat::compute_structural_tensors<Core::LinAlg::Matrix<6, 1>, numfib>(
      fibers_, fiber_structural_tensors_stress_, structural_tensor_strategy_);
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::compute_structural_tensors()
{
  Mat::compute_structural_tensors<Core::LinAlg::Matrix<3, 3>, numfib>(
      fibers_, fiber_structural_tensors_, structural_tensor_strategy_);
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::set_fibers(
    int gp, const std::array<Core::LinAlg::Matrix<3, 1>, numfib>& fibers)
{
  if (gp >= get_anisotropy()->get_number_of_gauss_points())
  {
    FOUR_C_THROW(
        "The current Gauss point %i is out of range of the expected number of Gauss points %i.", gp,
        get_anisotropy()->get_number_of_gauss_points());
  }

  if (fibers_.empty())
  {
    fibers_.resize(GetFibersPerElement());
  }

  // Store fibers
  for (std::vector<Core::LinAlg::Matrix<3, 1>>::size_type i = 0; i < fibers.size(); ++i)
  {
    fibers_[gp].at(i).Update(fibers.at(i));
  }

  // Compute needed structural tensors
  compute_needed_structural_tensors();

  // Call the notifier method
  on_fibers_initialized();
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::set_fibers(
    const std::vector<std::array<Core::LinAlg::Matrix<3, 1>, numfib>>& fibers)
{
  fibers_ = fibers;

  // Compute needed structural tensors
  compute_needed_structural_tensors();

  // Call the notifier method
  on_fibers_initialized();
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::compute_needed_structural_tensors()
{
  if (static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR))
  {
    compute_structural_tensors();
  }
  if (static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR_STRESS))
  {
    compute_structural_tensors_stress();
  }
}

template <unsigned int numfib>
const Core::LinAlg::Matrix<3, 1>& Mat::FiberAnisotropyExtension<numfib>::GetFiber(
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
const Core::LinAlg::Matrix<6, 1>&
Mat::FiberAnisotropyExtension<numfib>::get_structural_tensor_stress(int gp, int i) const
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
const Core::LinAlg::Matrix<3, 3>& Mat::FiberAnisotropyExtension<numfib>::GetStructuralTensor(
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
void Mat::FiberAnisotropyExtension<numfib>::PackAnisotropy(
    Core::Communication::PackBuffer& data) const
{
  PackFiberArray<Core::LinAlg::Matrix<3, 1>, numfib>(data, fibers_);
  PackFiberArray<Core::LinAlg::Matrix<6, 1>, numfib>(data, fiber_structural_tensors_stress_);
  PackFiberArray<Core::LinAlg::Matrix<3, 3>, numfib>(data, fiber_structural_tensors_);
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(tensor_flags_));
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(fiber_location_));
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  UnpackFiberArray<Core::LinAlg::Matrix<3, 1>, numfib>(position, data, fibers_);
  UnpackFiberArray<Core::LinAlg::Matrix<6, 1>, numfib>(
      position, data, fiber_structural_tensors_stress_);
  UnpackFiberArray<Core::LinAlg::Matrix<3, 3>, numfib>(position, data, fiber_structural_tensors_);
  tensor_flags_ =
      static_cast<uint_fast8_t>(Core::Communication::ParObject::ExtractInt(position, data));
  fiber_location_ =
      static_cast<FiberLocation>(Core::Communication::ParObject::ExtractInt(position, data));
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::set_fiber_location(FiberLocation location)
{
  fiber_location_ = location;
}

template <unsigned int numfib>
int Mat::FiberAnisotropyExtension<numfib>::get_virtual_gauss_point(const int gp) const
{
  if (fiber_location_ == FiberLocation::ElementFibers)
  {
    return GPDEFAULT;
  }

  return gp;
}

template <unsigned int numfib>
int Mat::FiberAnisotropyExtension<numfib>::GetFibersPerElement() const
{
  if (fiber_location_ == FiberLocation::ElementFibers)
  {
    return 1;
  }

  return get_anisotropy()->get_number_of_gauss_points();
}

// explicit instatiations of template classes
template class Mat::FiberAnisotropyExtension<1u>;
template class Mat::FiberAnisotropyExtension<2u>;
FOUR_C_NAMESPACE_CLOSE
