#include "4C_adapter_structure_scatra_ele.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_membrane_scatra.hpp"
#include "4C_shell7p_ele_scatra.hpp"
#include "4C_solid_scatra_3D_ele.hpp"
#include "4C_truss3_scatra.hpp"
#include "4C_w1_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Inpar::ScaTra::ImplType Adapter::get_sca_tra_impl_type(Core::Elements::Element* ele)
{
  Inpar::ScaTra::ImplType impltype(Inpar::ScaTra::impltype_undefined);

  // the element type name, needed to cast correctly in the following
  const std::string eletypename = ele->element_type().name();

  if (eletypename == "Wall1ScatraType")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::Wall1Scatra*>(ele))->impl_type();
  }
  // shell scatra elements
  else if (eletypename == "Shell7pScatraType")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::Shell7pScatra*>(ele))->impl_type();
  }
  // membrane3 scatra element
  else if (eletypename == "MembraneScatra_tri3Type")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>*>(ele))
                   ->impl_type();
  }
  // membrane6 scatra element
  else if (eletypename == "MembraneScatra_tri6Type")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>*>(ele))
                   ->impl_type();
  }
  // membrane4 scatra element
  else if (eletypename == "MembraneScatra_quad4Type")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>*>(ele))
                   ->impl_type();
  }
  // membrane9 scatra element
  else if (eletypename == "MembraneScatra_quad9Type")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>*>(ele))
                   ->impl_type();
  }
  // truss3 scatra element
  else if (eletypename == "Truss3ScatraType")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::Truss3Scatra*>(ele))->impl_type();
  }
  // SolidScatra element
  else if (eletypename == "SolidScatraType")
  {
    impltype = (dynamic_cast<Discret::ELEMENTS::SolidScatra*>(ele))->impl_type();
  }
  else
  {
    if (!(eletypename == "Bele3Type")) return impltype;
    impltype = Inpar::ScaTra::impltype_no_physics;
  }

  return impltype;
}
FOUR_C_NAMESPACE_CLOSE
