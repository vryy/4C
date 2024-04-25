/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of common functionality for anisotropic materials

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mat_anisotropy.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_fiber_nodal_fiber_holder.hpp"
#include "4C_mat_anisotropy_extension.hpp"
#include "4C_mat_anisotropy_utils.hpp"
#include "4C_mat_service.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::Anisotropy::Anisotropy()
    : element_fibers_initialized_(false),
      gp_fibers_initialized_(false),
      element_fibers_(0),
      gp_fibers_(0),
      gp_cylinder_coordinate_system_managers_(0),
      extensions_(0)
{
  // empty
}

void MAT::Anisotropy::PackAnisotropy(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::ParObject::AddtoPack(data, numgp_);
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(element_fibers_initialized_));
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(gp_fibers_initialized_));
  CORE::COMM::ParObject::AddtoPack(data, element_fibers_);
  PackFiberVector<CORE::LINALG::Matrix<3, 1>>(data, gp_fibers_);

  if (element_cylinder_coordinate_system_manager_)
  {
    CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(true));
    element_cylinder_coordinate_system_manager_->Pack(data);
  }
  else
  {
    CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(false));
  }

  for (const auto& gpCylinderCoordinateSystemManager : gp_cylinder_coordinate_system_managers_)
  {
    gpCylinderCoordinateSystemManager.Pack(data);
  }
}

void MAT::Anisotropy::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  CORE::COMM::ParObject::ExtractfromPack(position, data, numgp_);
  element_fibers_initialized_ =
      static_cast<bool>(CORE::COMM::ParObject::ExtractInt(position, data));
  gp_fibers_initialized_ = static_cast<bool>(CORE::COMM::ParObject::ExtractInt(position, data));
  CORE::COMM::ParObject::ExtractfromPack(position, data, element_fibers_);
  UnpackFiberVector<CORE::LINALG::Matrix<3, 1>>(position, data, gp_fibers_);

  if (static_cast<bool>(CORE::COMM::ParObject::ExtractInt(position, data)))
  {
    element_cylinder_coordinate_system_manager_ = CylinderCoordinateSystemManager();
    element_cylinder_coordinate_system_manager_->Unpack(data, position);
  }
  else
  {
    element_cylinder_coordinate_system_manager_ = std::nullopt;
  }

  for (auto& gpCylinderCoordinateSystemManager : gp_cylinder_coordinate_system_managers_)
  {
    gpCylinderCoordinateSystemManager.Unpack(data, position);
  }
}

void MAT::Anisotropy::SetNumberOfGaussPoints(int numgp) { numgp_ = numgp; }

void MAT::Anisotropy::ReadAnisotropyFromElement(INPUT::LineDefinition* lineDefinition)
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
    if (!element_cylinder_coordinate_system_manager_)
    {
      element_cylinder_coordinate_system_manager_ = CylinderCoordinateSystemManager();
    }

    element_cylinder_coordinate_system_manager_->ReadFromElementLineDefinition(lineDefinition);
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
    element_fibers_.resize(i + 1);
    ReadAnisotropyFiber(lineDefinition, "FIBER" + std::to_string(i + 1), element_fibers_[i]);
    i += 1;
  }

  OnElementFibersInitialized();
}

void MAT::Anisotropy::ReadAnisotropyFromParameterList(const Teuchos::ParameterList& params)
{
  if (params.isParameter("fiberholder"))
  {
    const auto& fiberHolder = params.get<DRT::FIBER::NodalFiberHolder>("fiberholder");

    gp_fibers_.resize(numgp_);

    for (const auto& fiber : fiberHolder.GetFibers())
    {
      InsertFibers(fiber);
    }
  }

  OnGPFibersInitialized();
}

void MAT::Anisotropy::InsertFibers(std::vector<CORE::LINALG::Matrix<3, 1>> fiber)
{
  for (unsigned gp = 0; gp < numgp_; ++gp)
  {
    gp_fibers_[gp].emplace_back(fiber[gp]);
  }
}

void MAT::Anisotropy::SetElementFibers(const std::vector<CORE::LINALG::Matrix<3, 1>>& fibers)
{
  element_fibers_ = fibers;

  OnElementFibersInitialized();
}

void MAT::Anisotropy::SetGaussPointFibers(
    const std::vector<std::vector<CORE::LINALG::Matrix<3, 1>>>& fibers)
{
  // check input fibers whether they make sense

  // Check whether the size of the first vector is the number of Gauss points
  if (fibers.size() != numgp_)
  {
    FOUR_C_THROW("The Gauss point fibers don't have the expected size of %d (%d given).", numgp_,
        fibers.size());
  }

  // Check whether every second vector have the same lenghts
  unsigned num_fibs = 1;
  unsigned i = 0;
  for (const auto& gpfibers : fibers)
  {
    if (i == 0)
    {
      num_fibs = gpfibers.size();
    }
    else if (num_fibs != gpfibers.size())
    {
      FOUR_C_THROW(
          "The size of the Gauss point do not match! At every Gauss point, the same amount of "
          "fibers are necessary. Error occured at Gauss point %d. Expected %d fibers, but got %d.",
          i, num_fibs, gpfibers.size());
    }
  }

  gp_fibers_ = fibers;

  OnGPFibersInitialized();
}

const CORE::LINALG::Matrix<3, 1>& MAT::Anisotropy::GetElementFiber(unsigned int i) const
{
  if (!element_fibers_initialized_)
  {
    FOUR_C_THROW("The element fibers are not yet initialized.");
  }
  if (i >= element_fibers_.size())
  {
    FOUR_C_THROW(
        "You requested fiber %d, but only %d fibers are available", i + 1, element_fibers_.size());
  }
  return element_fibers_[i];
}

const std::vector<CORE::LINALG::Matrix<3, 1>>& MAT::Anisotropy::GetElementFibers() const
{
  if (!element_fibers_initialized_)
  {
    FOUR_C_THROW("The element fibers are not yet initialized.");
  }
  return element_fibers_;
}

const std::vector<std::vector<CORE::LINALG::Matrix<3, 1>>>& MAT::Anisotropy::GetGPFibers() const
{
  if (!gp_fibers_initialized_)
  {
    FOUR_C_THROW("The Gauss point fibers are not yet initialized.");
  }
  return gp_fibers_;
}

const CORE::LINALG::Matrix<3, 1>& MAT::Anisotropy::GetGPFiber(unsigned int gp, unsigned int i) const
{
  if (!gp_fibers_initialized_)
  {
    FOUR_C_THROW("The GP fibers are not yet initialized.");
  }

  if (gp >= gp_fibers_.size())
  {
    FOUR_C_THROW("The number of GP is too large. %d instead of maximum allowed %d", gp + 1,
        gp_fibers_.size());
  }

  if (i >= gp_fibers_[gp].size())
  {
    FOUR_C_THROW(
        "You requested fiber %d, but only %d fibers are available", i + 1, element_fibers_.size());
  }
  return gp_fibers_[gp][i];
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

int MAT::Anisotropy::GetNumberOfElementFibers() const { return element_fibers_.size(); }

int MAT::Anisotropy::GetNumberOfGPFibers() const
{
  if (gp_fibers_.empty()) return 0;

  return gp_fibers_[0].size();
}

bool MAT::Anisotropy::HasElementCylinderCoordinateSystem() const
{
  return element_cylinder_coordinate_system_manager_.has_value();
}

bool MAT::Anisotropy::HasGPCylinderCoordinateSystem() const
{
  return !gp_cylinder_coordinate_system_managers_.empty();
}
FOUR_C_NAMESPACE_CLOSE
