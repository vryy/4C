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

MAT::FiberAnisotropyExtension::FiberAnisotropyExtension(
    const Teuchos::RCP<ELASTIC::StructuralTensorStrategyBase>& stucturalTensorStrategy)
    : tensor_flags_(0),
      fibers_(0),
      fiberStructuralTensors_stress_(0),
      fiberStructuralTensors_(0),
      structuralTensorStrategy_(stucturalTensorStrategy)
{
}

MAT::FiberAnisotropyExtension::FiberAnisotropyExtension()
    : tensor_flags_(0),
      fibers_(0),
      fiberStructuralTensors_stress_(0),
      fiberStructuralTensors_(0),
      structuralTensorStrategy_(Teuchos::null)
{
}

void MAT::FiberAnisotropyExtension::ComputeStructuralTensors_stress()
{
  MAT::ComputeStructuralTensors<>(
      fibers_, fiberStructuralTensors_stress_, structuralTensorStrategy_);
}

void MAT::FiberAnisotropyExtension::ComputeStructuralTensors()
{
  MAT::ComputeStructuralTensors<>(fibers_, fiberStructuralTensors_, structuralTensorStrategy_);
}

void MAT::FiberAnisotropyExtension::SetFibers(
    int gp, const std::vector<LINALG::Matrix<3, 1>>& fibers)
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
    if (fibers_[gp].empty())
    {
      fibers_[gp].resize(fibers.size());
    }

    fibers_[gp][i].Update(fibers[i]);
  }

  // Compute needed structural tensors
  ComputeNeededStructuralTensors();

  // Call the notifier method
  OnFibersInitialized();
}

void MAT::FiberAnisotropyExtension::SetFibers(
    const std::vector<std::vector<LINALG::Matrix<3, 1>>>& fibers)
{
  fibers_ = fibers;

  // Compute needed structural tensors
  ComputeNeededStructuralTensors();

  // Call the notifier method
  OnFibersInitialized();
}

void MAT::FiberAnisotropyExtension::ComputeNeededStructuralTensors()
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

const LINALG::Matrix<3, 1>& MAT::FiberAnisotropyExtension::GetFiber(int gp, int i) const
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

  // This return statement should never be reached. However, the compiler enforces me to have it
  // here
  return fibers_[0][0];
}

const LINALG::Matrix<6, 1>& MAT::FiberAnisotropyExtension::GetStructuralTensor_stress(
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

  return fiberStructuralTensors_stress_[0][0];
}

const LINALG::Matrix<3, 3>& MAT::FiberAnisotropyExtension::GetStructuralTensor(int gp, int i) const
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

  return fiberStructuralTensors_[0][0];
}

void MAT::FiberAnisotropyExtension::PackAnisotropy(DRT::PackBuffer& data) const
{
  PackFiberVector(data, fibers_);
  PackFiberVector(data, fiberStructuralTensors_stress_);
  PackFiberVector(data, fiberStructuralTensors_);
  DRT::ParObject::AddtoPack(data, static_cast<int>(tensor_flags_));
  DRT::ParObject::AddtoPack(data, static_cast<int>(fiberLocation_));
}

void MAT::FiberAnisotropyExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  UnpackFiberVector<LINALG::Matrix<3, 1>>(position, data, fibers_);
  UnpackFiberVector<LINALG::Matrix<6, 1>>(position, data, fiberStructuralTensors_stress_);
  UnpackFiberVector<LINALG::Matrix<3, 3>>(position, data, fiberStructuralTensors_);
  tensor_flags_ = static_cast<uint_fast8_t>(DRT::ParObject::ExtractInt(position, data));
  fiberLocation_ = static_cast<FiberLocation>(DRT::ParObject::ExtractInt(position, data));
}

void MAT::FiberAnisotropyExtension::SetFiberLocation(FiberLocation location)
{
  fiberLocation_ = location;
}

int MAT::FiberAnisotropyExtension::GetVirtualGaussPoint(const Teuchos::ParameterList& params) const
{
  if (fiberLocation_ == FiberLocation::ElementFibers)
  {
    return GPDEFAULT;
  }

  return params.get<int>("gp");
}

int MAT::FiberAnisotropyExtension::GetVirtualGaussPoint(const int gp) const
{
  if (fiberLocation_ == FiberLocation::ElementFibers)
  {
    return GPDEFAULT;
  }

  return gp;
}

int MAT::FiberAnisotropyExtension::GetFibersPerElement() const
{
  if (fiberLocation_ == FiberLocation::ElementFibers)
  {
    return 1;
  }

  return GetAnisotropy()->GetNumberOfGaussPoints();
}