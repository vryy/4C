/*----------------------------------------------------------------------*/
/*! \file
\brief Solid Tet10 Element
\level 2
*----------------------------------------------------------------------*/

#include "so_tet10.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet10::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  SolidMaterial()->Setup(NUMGPT_SOTET10, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::kinem_linear;
    dserror("Reading of SO_TET10 element failed only nonlinear kinematics implemented");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  // geometrically non-linear with Updated Lagrangean approach
  else
    dserror("Reading of SO_TET10 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}
