/*!----------------------------------------------------------------------
\brief pyramid shaped solid element
\level 2

\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so_pyramid5.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_pyramid5::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  Teuchos::RCP<MAT::Material> mat = Material();

  SolidMaterial()->Setup(NUMGPT_SOP5, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);

  if (buffer == "linear")
  {
    kintype_ = INPAR::STR::kinem_linear;
  }
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
  }
  else
    dserror("Reading SO_PYRAMID5 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    dserror("This element currently does not support the extended update call.");

  // only for linear SVK materials and small strain plastic materials
  bool admissibl_mat = false;
  if ((mat->MaterialType() == INPAR::MAT::m_stvenant) or
      (mat->MaterialType() == INPAR::MAT::m_thermostvenant) or
      (mat->MaterialType() == INPAR::MAT::m_pllinelast) or
      (mat->MaterialType() == INPAR::MAT::m_thermopllinelast) or
      (mat->MaterialType() == INPAR::MAT::m_elpldamage))
    admissibl_mat = true;

  // check for SVK material if geometrically linear
  if ((kintype_ == INPAR::STR::kinem_linear) and (admissibl_mat == false))
    dserror("ERROR: Only linear elasticity (SVK) for geometrically linear pyramid5 element");

  return true;
}
