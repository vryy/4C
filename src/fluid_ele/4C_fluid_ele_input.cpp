/*----------------------------------------------------------------------*/
/*! \file

\brief read_element method of the fluid element implementation


\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Fluid::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));

  // set discretization type (setOptimalgaussrule is pushed into element
  // routine)
  set_dis_type(Core::FE::StringToCellType(distype));

  std::string na = container.get<std::string>("NA");
  if (na == "ale" or na == "ALE" or na == "Ale")
  {
    is_ale_ = true;
  }
  else if (na == "euler" or na == "EULER" or na == "Euler")
    is_ale_ = false;
  else
    FOUR_C_THROW("Reading of fluid element failed: Euler/Ale");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
