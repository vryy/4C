// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_anisotropy_extension.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_mat_anisotropy_utils.hpp"

#include <Teuchos_SerialDenseSolver.hpp>

FOUR_C_NAMESPACE_OPEN

template <unsigned int numfib>
Mat::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension(
    const std::shared_ptr<Elastic::StructuralTensorStrategyBase>& stucturalTensorStrategy)
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
      structural_tensor_strategy_(nullptr)
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
    fibers_.resize(get_fibers_per_element());
  }

  // Store fibers
  for (std::vector<Core::LinAlg::Matrix<3, 1>>::size_type i = 0; i < fibers.size(); ++i)
  {
    fibers_[gp].at(i).update(fibers.at(i));
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
const Core::LinAlg::Matrix<3, 1>& Mat::FiberAnisotropyExtension<numfib>::get_fiber(
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
const Core::LinAlg::Matrix<3, 3>& Mat::FiberAnisotropyExtension<numfib>::get_structural_tensor(
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
void Mat::FiberAnisotropyExtension<numfib>::pack_anisotropy(
    Core::Communication::PackBuffer& data) const
{
  pack_fiber_array<Core::LinAlg::Matrix<3, 1>, numfib>(data, fibers_);
  pack_fiber_array<Core::LinAlg::Matrix<6, 1>, numfib>(data, fiber_structural_tensors_stress_);
  pack_fiber_array<Core::LinAlg::Matrix<3, 3>, numfib>(data, fiber_structural_tensors_);
  add_to_pack(data, tensor_flags_);
  add_to_pack(data, fiber_location_);
}

template <unsigned int numfib>
void Mat::FiberAnisotropyExtension<numfib>::unpack_anisotropy(
    Core::Communication::UnpackBuffer& buffer)
{
  unpack_fiber_array<Core::LinAlg::Matrix<3, 1>, numfib>(buffer, fibers_);
  unpack_fiber_array<Core::LinAlg::Matrix<6, 1>, numfib>(buffer, fiber_structural_tensors_stress_);
  unpack_fiber_array<Core::LinAlg::Matrix<3, 3>, numfib>(buffer, fiber_structural_tensors_);
  extract_from_pack(buffer, tensor_flags_);
  extract_from_pack(buffer, fiber_location_);
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
int Mat::FiberAnisotropyExtension<numfib>::get_fibers_per_element() const
{
  if (fiber_location_ == FiberLocation::ElementFibers)
  {
    return 1;
  }

  return get_anisotropy()->get_number_of_gauss_points();
}

// explicit instantiations of template classes
template class Mat::FiberAnisotropyExtension<1u>;
template class Mat::FiberAnisotropyExtension<2u>;
FOUR_C_NAMESPACE_CLOSE
