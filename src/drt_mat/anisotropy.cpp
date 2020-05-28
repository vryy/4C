/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of common functionality for anisotropic materials

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/
#include "anisotropy.H"
#include "material_service.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_fiber/nodal_fiber_holder.H"
#include "anisotropy_utils.H"
#include "anisotropy_extension.H"

MAT::Anisotropy::Anisotropy()
    : element_fibers_initialized_(false),
      gp_fibers_initialized_(false),
      elementFibers_(0),
      gpFibers_(0),
      elementCylinderCoordinateSystemManager_(),
      gpCylinderCoordinateSystemManagers_(0),
      extensions_(0)
{
  // empty
}

void MAT::Anisotropy::PackAnisotropy(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, numgp_);
  DRT::ParObject::AddtoPack(data, static_cast<int>(element_fibers_initialized_));
  DRT::ParObject::AddtoPack(data, static_cast<int>(gp_fibers_initialized_));
  DRT::ParObject::AddtoPack(data, elementFibers_);
  PackFiberVector<LINALG::Matrix<3, 1>>(data, gpFibers_);

  if (elementCylinderCoordinateSystemManager_)
  {
    DRT::ParObject::AddtoPack(data, static_cast<int>(true));
    elementCylinderCoordinateSystemManager_->Pack(data);
  }
  else
  {
    DRT::ParObject::AddtoPack(data, static_cast<int>(false));
  }

  for (const auto& gpCylinderCoordinateSystemManager : gpCylinderCoordinateSystemManagers_)
  {
    gpCylinderCoordinateSystemManager.Pack(data);
  }
}

void MAT::Anisotropy::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DRT::ParObject::ExtractfromPack(position, data, numgp_);
  element_fibers_initialized_ = static_cast<bool>(DRT::ParObject::ExtractInt(position, data));
  gp_fibers_initialized_ = static_cast<bool>(DRT::ParObject::ExtractInt(position, data));
  DRT::ParObject::ExtractfromPack(position, data, elementFibers_);
  UnpackFiberVector<LINALG::Matrix<3, 1>>(position, data, gpFibers_);

  if (static_cast<bool>(DRT::ParObject::ExtractInt(position, data)))
  {
    elementCylinderCoordinateSystemManager_ = CylinderCoordinateSystemManager();
    elementCylinderCoordinateSystemManager_->Unpack(data, position);
  }
  else
  {
    elementCylinderCoordinateSystemManager_ = boost::none;
  }

  for (auto& gpCylinderCoordinateSystemManager : gpCylinderCoordinateSystemManagers_)
  {
    gpCylinderCoordinateSystemManager.Unpack(data, position);
  }
}

void MAT::Anisotropy::SetNumberOfGaussPoints(int numgp) { numgp_ = numgp; }

void MAT::Anisotropy::ReadAnisotropyFromElement(DRT::INPUT::LineDefinition* lineDefinition)
{
  if (lineDefinition == nullptr)
  {
    // Line definition is not given, so I cannot read anything from the line definition
    return;
  }
  // Read coordinate system

  if (lineDefinition->HaveNamed("RAD") and lineDefinition->HaveNamed("AXI") and
      lineDefinition->HaveNamed("CIR"))
  {
    if (!elementCylinderCoordinateSystemManager_)
    {
      elementCylinderCoordinateSystemManager_ = CylinderCoordinateSystemManager();
    }

    elementCylinderCoordinateSystemManager_->ReadFromElementLineDefinition(lineDefinition);
  }

  // read fibers in FIBERi notation
  // determine number of fibers
  unsigned i = 0;
  while (true)
  {
    if (!lineDefinition->HaveNamed("FIBER" + std::to_string(i + 1)))
    {
      break;
    }
    elementFibers_.resize(i + 1);
    ReadAnisotropyFiber(lineDefinition, "FIBER" + std::to_string(i + 1), elementFibers_[i]);
    i += 1;
  }

  OnElementFibersInitialized();
}

void MAT::Anisotropy::ReadAnisotropyFromParameterList(const Teuchos::ParameterList& params)
{
  // ToDo: Implement GP anisotropy
  if (params.isParameter("fiberholder"))
  {
    const auto& fiberHolder = params.get<DRT::FIBER::NodalFiberHolder>("fiberholder");

    std::vector<LINALG::Matrix<3, 1>> gpfiber1 =
        fiberHolder.GetFiber(DRT::FIBER::FiberType::Fiber1);
    std::vector<LINALG::Matrix<3, 1>> gpfiber2 =
        fiberHolder.GetFiber(DRT::FIBER::FiberType::Fiber2);

    gpFibers_.resize(numgp_);
    for (int gp = 0; gp < numgp_; ++gp)
    {
      gpFibers_[gp].resize(2);
      gpFibers_[gp][0].Update(gpfiber1.at(gp));
      gpFibers_[gp][1].Update(gpfiber2.at(gp));
    }
  }

  OnGPFibersInitialized();
}

const LINALG::Matrix<3, 1>& MAT::Anisotropy::GetElementFiber(unsigned int i) const
{
  if (!element_fibers_initialized_)
  {
    dserror("The element fibers are not yet initialized.");
  }
  if (i >= elementFibers_.size())
  {
    dserror(
        "You requested fiber %d, but only %d fibers are available", i + 1, elementFibers_.size());
  }
  return elementFibers_[i];
}

const std::vector<LINALG::Matrix<3, 1>>& MAT::Anisotropy::GetElementFibers() const
{
  if (!element_fibers_initialized_)
  {
    dserror("The element fibers are not yet initialized.");
  }
  return elementFibers_;
}

const std::vector<std::vector<LINALG::Matrix<3, 1>>>& MAT::Anisotropy::GetGPFibers() const
{
  if (!gp_fibers_initialized_)
  {
    dserror("The Gauss point fibers are not yet initialized.");
  }
  return gpFibers_;
}

const LINALG::Matrix<3, 1>& MAT::Anisotropy::GetGPFiber(unsigned int gp, unsigned int i) const
{
  if (!gp_fibers_initialized_)
  {
    dserror("The GP fibers are not yet initialized.");
  }

  if (gp >= gpFibers_.size())
  {
    dserror("The number of GP is too large. %d instead of maximum allowed %d", gp + 1,
        gpFibers_.size());
  }

  if (i >= gpFibers_[gp].size())
  {
    dserror(
        "You requested fiber %d, but only %d fibers are available", i + 1, elementFibers_.size());
  }
  return gpFibers_[gp][i];
}

void MAT::Anisotropy::RegisterAnisotropyExtension(BaseAnisotropyExtension& extension)
{
  extensions_.emplace_back(Teuchos::rcpFromRef(extension));
  extension.SetAnisotropy(*this);
}

void MAT::Anisotropy::OnElementFibersInitialized()
{
  element_fibers_initialized_ = true;
  for (auto& extension : extensions_)
  {
    extension->OnGlobalElementDataInitialized();
  }

  if (element_fibers_initialized_ and gp_fibers_initialized_)
  {
    for (auto& extension : extensions_)
    {
      extension->OnGlobalDataInitialized();
    }
  }
}

void MAT::Anisotropy::OnGPFibersInitialized()
{
  gp_fibers_initialized_ = true;
  for (auto& extension : extensions_)
  {
    extension->OnGlobalGPDataInitialized();
  }

  if (element_fibers_initialized_ and gp_fibers_initialized_)
  {
    for (auto& extension : extensions_)
    {
      extension->OnGlobalDataInitialized();
    }
  }
}

int MAT::Anisotropy::GetNumberOfGaussPoints() const { return numgp_; }

int MAT::Anisotropy::GetNumberOfElementFibers() const { return elementFibers_.size(); }

int MAT::Anisotropy::GetNumberOfGPFibers() const
{
  if (gpFibers_.empty()) return 0;

  return gpFibers_[0].size();
}

bool MAT::Anisotropy::HasElementCylinderCoordinateSystem() const
{
  return elementCylinderCoordinateSystemManager_.is_initialized();
}

bool MAT::Anisotropy::HasGPCylinderCoordinateSystem() const
{
  return !gpCylinderCoordinateSystemManagers_.empty();
}