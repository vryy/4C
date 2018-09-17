/*!----------------------------------------------------------------------**###
\file so_tet4_input.cpp
\brief
\level 1

<pre>
\maintainer Christoph Meier
</pre>

*----------------------------------------------------------------------*/

#include "so_tet4.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet4::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  SolidMaterial()->Setup(NUMGPT_SOTET4, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::kinem_linear;
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else
  {
    dserror("Reading of SO_TET4 element failed KINEM unknown");
  }

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  return true;
}
