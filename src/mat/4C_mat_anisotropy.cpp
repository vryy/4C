/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of common functionality for anisotropic materials

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_mat_anisotropy.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_fem_general_fiber_node_holder.hpp"
#include "4C_mat_anisotropy_extension.hpp"
#include "4C_mat_anisotropy_utils.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Anisotropy::Anisotropy()
    : element_fibers_initialized_(false),
      gp_fibers_initialized_(false),
      element_fibers_(0),
      gp_fibers_(0),
      gp_cylinder_coordinate_system_managers_(0),
      extensions_(0)
{
  // empty
}

void Mat::Anisotropy::pack_anisotropy(Core::Communication::PackBuffer& data) const
{
  Core::Communication::ParObject::add_to_pack(data, numgp_);
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(element_fibers_initialized_));
  Core::Communication::ParObject::add_to_pack(data, static_cast<int>(gp_fibers_initialized_));
  Core::Communication::ParObject::add_to_pack(data, element_fibers_);
  pack_fiber_vector<Core::LinAlg::Matrix<3, 1>>(data, gp_fibers_);

  if (element_cylinder_coordinate_system_manager_)
  {
    Core::Communication::ParObject::add_to_pack(data, static_cast<int>(true));
    element_cylinder_coordinate_system_manager_->pack(data);
  }
  else
  {
    Core::Communication::ParObject::add_to_pack(data, static_cast<int>(false));
  }

  for (const auto& gpCylinderCoordinateSystemManager : gp_cylinder_coordinate_system_managers_)
  {
    gpCylinderCoordinateSystemManager.pack(data);
  }
}

void Mat::Anisotropy::unpack_anisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  Core::Communication::ParObject::extract_from_pack(position, data, numgp_);
  element_fibers_initialized_ =
      static_cast<bool>(Core::Communication::ParObject::extract_int(position, data));
  gp_fibers_initialized_ =
      static_cast<bool>(Core::Communication::ParObject::extract_int(position, data));
  Core::Communication::ParObject::extract_from_pack(position, data, element_fibers_);
  unpack_fiber_vector<Core::LinAlg::Matrix<3, 1>>(position, data, gp_fibers_);

  if (static_cast<bool>(Core::Communication::ParObject::extract_int(position, data)))
  {
    element_cylinder_coordinate_system_manager_ = CylinderCoordinateSystemManager();
    element_cylinder_coordinate_system_manager_->unpack(data, position);
  }
  else
  {
    element_cylinder_coordinate_system_manager_ = std::nullopt;
  }

  for (auto& gpCylinderCoordinateSystemManager : gp_cylinder_coordinate_system_managers_)
  {
    gpCylinderCoordinateSystemManager.unpack(data, position);
  }
}

void Mat::Anisotropy::set_number_of_gauss_points(int numgp) { numgp_ = numgp; }

void Mat::Anisotropy::read_anisotropy_from_element(
    const Core::IO::InputParameterContainer& container)
{
  // Read coordinate system
  if (container.get_if<std::vector<double>>("RAD") != nullptr and
      container.get_if<std::vector<double>>("AXI") != nullptr and
      container.get_if<std::vector<double>>("CIR") != nullptr)
  {
    // read fibers in RAD AXI CIR notation
    if (!element_cylinder_coordinate_system_manager_)
    {
      element_cylinder_coordinate_system_manager_ = CylinderCoordinateSystemManager();
    }

    element_cylinder_coordinate_system_manager_->read_from_element_line_definition(container);
  }

  else
  {
    // read fibers in FIBERi notation and determine number of fibers
    for (int fiber_index = 1;; ++fiber_index)
    {
      std::string fiber_name = "FIBER" + std::to_string(fiber_index);
      if (container.get_if<std::vector<double>>(fiber_name) == nullptr)
      {
        break;
      }

      element_fibers_.resize(fiber_index);
      read_anisotropy_fiber(container, fiber_name, element_fibers_[fiber_index - 1]);
    }
  }
  on_element_fibers_initialized();
}

void Mat::Anisotropy::read_anisotropy_from_parameter_list(const Teuchos::ParameterList& params)
{
  if (params.isParameter("fiberholder"))
  {
    const auto& fiberHolder = params.get<Core::Nodes::NodalFiberHolder>("fiberholder");

    gp_fibers_.resize(numgp_);

    for (const auto& fiber : fiberHolder.get_fibers())
    {
      insert_fibers(fiber);
    }
  }

  on_gp_fibers_initialized();
}

void Mat::Anisotropy::insert_fibers(std::vector<Core::LinAlg::Matrix<3, 1>> fiber)
{
  for (unsigned gp = 0; gp < numgp_; ++gp)
  {
    gp_fibers_[gp].emplace_back(fiber[gp]);
  }
}

void Mat::Anisotropy::set_element_fibers(const std::vector<Core::LinAlg::Matrix<3, 1>>& fibers)
{
  element_fibers_ = fibers;

  on_element_fibers_initialized();
}

void Mat::Anisotropy::set_gauss_point_fibers(
    const std::vector<std::vector<Core::LinAlg::Matrix<3, 1>>>& fibers)
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

  on_gp_fibers_initialized();
}

const Core::LinAlg::Matrix<3, 1>& Mat::Anisotropy::get_element_fiber(unsigned int i) const
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

const std::vector<Core::LinAlg::Matrix<3, 1>>& Mat::Anisotropy::get_element_fibers() const
{
  if (!element_fibers_initialized_)
  {
    FOUR_C_THROW("The element fibers are not yet initialized.");
  }
  return element_fibers_;
}

const std::vector<std::vector<Core::LinAlg::Matrix<3, 1>>>&
Mat::Anisotropy::get_gauss_point_fibers() const
{
  if (!gp_fibers_initialized_)
  {
    FOUR_C_THROW("The Gauss point fibers are not yet initialized.");
  }
  return gp_fibers_;
}

const Core::LinAlg::Matrix<3, 1>& Mat::Anisotropy::get_gauss_point_fiber(
    unsigned int gp, unsigned int i) const
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

void Mat::Anisotropy::register_anisotropy_extension(BaseAnisotropyExtension& extension)
{
  extensions_.emplace_back(Teuchos::rcpFromRef(extension));
  extension.set_anisotropy(*this);
}

void Mat::Anisotropy::on_element_fibers_initialized()
{
  element_fibers_initialized_ = true;
  for (auto& extension : extensions_)
  {
    extension->on_global_element_data_initialized();
  }

  if (element_fibers_initialized_ and gp_fibers_initialized_)
  {
    for (auto& extension : extensions_)
    {
      extension->on_global_data_initialized();
    }
  }
}

void Mat::Anisotropy::on_gp_fibers_initialized()
{
  gp_fibers_initialized_ = true;
  for (auto& extension : extensions_)
  {
    extension->on_global_gp_data_initialized();
  }

  if (element_fibers_initialized_ and gp_fibers_initialized_)
  {
    for (auto& extension : extensions_)
    {
      extension->on_global_data_initialized();
    }
  }
}

int Mat::Anisotropy::get_number_of_gauss_points() const { return numgp_; }

int Mat::Anisotropy::get_number_of_element_fibers() const { return element_fibers_.size(); }

int Mat::Anisotropy::get_number_of_gauss_point_fibers() const
{
  if (gp_fibers_.empty()) return 0;

  return gp_fibers_[0].size();
}

bool Mat::Anisotropy::has_element_cylinder_coordinate_system() const
{
  return element_cylinder_coordinate_system_manager_.has_value();
}

bool Mat::Anisotropy::has_gp_cylinder_coordinate_system() const
{
  return !gp_cylinder_coordinate_system_managers_.empty();
}
FOUR_C_NAMESPACE_CLOSE
