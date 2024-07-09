/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_truss3.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Truss3::read_element(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material_id = 0;
  linedef->extract_int("MAT", material_id);
  set_material(0, Mat::Factory(material_id));

  linedef->extract_double("CROSS", crosssec_);

  std::string buffer;
  linedef->extract_string("KINEM", buffer);

  if (buffer == "totlag")  // geometrically non-linear with Total Lagrangean approach
    kintype_ = KinematicType::tr3_totlag;
  else if (buffer == "engstr")  // geometrically non-linear approach with engineering strains
    kintype_ = KinematicType::tr3_engstrain;
  else
    FOUR_C_THROW("Reading of Truss3 element failed because of unknown kinematic type!");

  return true;
}

/*------------------------------------------------------------------------*
 | Set cross section area                           (public) mueller 03/12|
 *------------------------------------------------------------------------*/
void Discret::ELEMENTS::Truss3::set_cross_sec(const double& crosssec) { crosssec_ = crosssec; }

FOUR_C_NAMESPACE_CLOSE
