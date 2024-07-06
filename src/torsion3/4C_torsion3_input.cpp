/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional torsion spring element

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_torsion3.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Torsion3::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read type of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);
  set_material(0, Mat::Factory(material_id));

  // read type of bending potential
  std::string buffer;
  linedef->extract_string("BENDINGPOTENTIAL", buffer);

  // bending potential E_bend = 0.5*SPRING*\theta^2
  if (buffer == "quadratic") bendingpotential_ = quadratic;

  // bending potential E_bend = SPRING*(1 - \cos(\theta^2) )
  else if (buffer == "cosine")
    bendingpotential_ = cosine;

  else
    FOUR_C_THROW("Reading of Torsion3 element failed because of unknown potential type!");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
