/*----------------------------------------------------------------------*/
/*! \file

\brief input-related methods of the quadratic NURBS 27 element

\level 2


*----------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_so3_nurbs27.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Discret::ELEMENTS::Nurbs::SoNurbs27::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  const int numgp = 27;
  SolidMaterial()->Setup(numgp, linedef);

  // read possible gaussian points, obsolete for computation
  std::vector<int> ngp;
  linedef->ExtractIntVector("GP", ngp);
  for (int i = 0; i < 3; ++i)
    if (ngp[i] != 3) FOUR_C_THROW("Only version with 3 GP for So_N27 implemented");

  // we expect kintype to be total lagrangian
  kintype_ = Inpar::STR::KinemType::nonlinearTotLag;

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  return true;
}

FOUR_C_NAMESPACE_CLOSE
