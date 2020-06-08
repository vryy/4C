/*----------------------------------------------------------------------*/
/*! \file

\brief

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include <Epetra_SerialDenseSolver.h>
#include "anisotropy_extension.H"
#include "anisotropy_utils.H"
#include "../drt_lib/drt_parobject.H"

template <unsigned int numfib>
MAT::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension(
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy)
    : fibers_(0),
      fiberStructuralTensors_stress_(0),
      fiberStructuralTensors_(0),
      structuralTensorStrategy_(stucturalTensorStrategy)
{
}

template <unsigned int numfib>
MAT::FiberAnisotropyExtension<numfib>::FiberAnisotropyExtension()
    : fibers_(0),
      fiberStructuralTensors_stress_(0),
      fiberStructuralTensors_(0),
      structuralTensorStrategy_(Teuchos::null)
{
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::ComputeStructuralTensors_stress()
{
  MAT::ComputeStructuralTensors<LINALG::Matrix<6, 1>, numfib>(
      fibers_, fiberStructuralTensors_stress_, structuralTensorStrategy_);
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::ComputeStructuralTensors()
{
  MAT::ComputeStructuralTensors<LINALG::Matrix<3, 3>, numfib>(
      fibers_, fiberStructuralTensors_, structuralTensorStrategy_);
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::SetFibers(
    int gp, const std::array<LINALG::Matrix<3, 1>, numfib>& fibers)
{
  if (gp >= GetAnisotropy()->GetNumberOfGaussPoints())
  {
    dserror("The current Gauss point %i is out of range of the expected number of Gauss points %i.",
        gp, GetAnisotropy()->GetNumberOfGaussPoints());
  }

  if (fibers_.empty())
  {
    fibers_.resize(GetFibersPerElement());
  }

  // Store fibers
  for (std::vector<LINALG::Matrix<3, 1>>::size_type i = 0; i < fibers.size(); ++i)
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
    const std::vector<std::array<LINALG::Matrix<3, 1>, numfib>>& fibers)
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
const LINALG::Matrix<3, 1>& MAT::FiberAnisotropyExtension<numfib>::GetFiber(int gp, int i) const
{
  switch (fiberLocation_)
  {
    case FiberLocation::ElementFibers:
      return fibers_[GPDEFAULT][i];
    case FiberLocation::GPFibers:
      return fibers_[gp][i];
    default:
      dserror("You have not specified, whether you want fibers on GP level or on element level.");
  }

  // just for compilation reasons. We will never land here because of the dserror() above
  std::abort();
}

template <unsigned int numfib>
const LINALG::Matrix<6, 1>& MAT::FiberAnisotropyExtension<numfib>::GetStructuralTensor_stress(
    int gp, int i) const
{
  if (not static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR_STRESS))
  {
    dserror("You have not specified that you need the fiber vector.");
  }
  switch (fiberLocation_)
  {
    case FiberLocation::ElementFibers:
      return fiberStructuralTensors_stress_[GPDEFAULT][i];
    case FiberLocation::GPFibers:
      return fiberStructuralTensors_stress_[gp][i];
    default:
      dserror("You have not specified, whether you want fibers on GP level or on element level.");
  }

  // just for compilation reasons. We will never land here because of the dserror() above
  std::abort();
}

template <unsigned int numfib>
const LINALG::Matrix<3, 3>& MAT::FiberAnisotropyExtension<numfib>::GetStructuralTensor(
    int gp, int i) const
{
  if (not static_cast<bool>(tensor_flags_ & STRUCTURAL_TENSOR))
  {
    dserror("You have not specified that you need the structural tensor.");
  }
  switch (fiberLocation_)
  {
    case FiberLocation::ElementFibers:
      return fiberStructuralTensors_[GPDEFAULT][i];
    case FiberLocation::GPFibers:
      return fiberStructuralTensors_[gp][i];
    default:
      dserror("You have not specified, whether you want fibers on GP level or on element level.");
  }

  // just for compilation reasons. We will never land here because of the dserror() above
  std::abort();
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::PackAnisotropy(DRT::PackBuffer& data) const
{
  PackFiberArray<LINALG::Matrix<3, 1>, numfib>(data, fibers_);
  PackFiberArray<LINALG::Matrix<6, 1>, numfib>(data, fiberStructuralTensors_stress_);
  PackFiberArray<LINALG::Matrix<3, 3>, numfib>(data, fiberStructuralTensors_);
  DRT::ParObject::AddtoPack(data, static_cast<int>(tensor_flags_));
  DRT::ParObject::AddtoPack(data, static_cast<int>(fiberLocation_));
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  UnpackFiberArray<LINALG::Matrix<3, 1>, numfib>(position, data, fibers_);
  UnpackFiberArray<LINALG::Matrix<6, 1>, numfib>(position, data, fiberStructuralTensors_stress_);
  UnpackFiberArray<LINALG::Matrix<3, 3>, numfib>(position, data, fiberStructuralTensors_);
  tensor_flags_ = static_cast<uint_fast8_t>(DRT::ParObject::ExtractInt(position, data));
  fiberLocation_ = static_cast<FiberLocation>(DRT::ParObject::ExtractInt(position, data));
}

template <unsigned int numfib>
void MAT::FiberAnisotropyExtension<numfib>::SetFiberLocation(FiberLocation location)
{
  fiberLocation_ = location;
}

template <unsigned int numfib>
int MAT::FiberAnisotropyExtension<numfib>::GetVirtualGaussPoint(const int gp) const
{
  if (fiberLocation_ == FiberLocation::ElementFibers)
  {
    return GPDEFAULT;
  }

  return gp;
}

template <unsigned int numfib>
int MAT::FiberAnisotropyExtension<numfib>::GetFibersPerElement() const
{
  if (fiberLocation_ == FiberLocation::ElementFibers)
  {
    return 1;
  }

  return GetAnisotropy()->GetNumberOfGaussPoints();
}

// explicit instatiations of template classes
template class MAT::FiberAnisotropyExtension<1u>;
template class MAT::FiberAnisotropyExtension<2u>;