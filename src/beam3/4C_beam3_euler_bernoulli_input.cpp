/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear torsionless rod based on a C1 curve

\level 2

*/
/*----------------------------------------------------------------------------*/

#include "4C_beam3_euler_bernoulli.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3eb::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, MAT::Factory(material));

  const auto mat_type = Material()->Parameter()->Type();
  FOUR_C_THROW_UNLESS(
      mat_type == CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper ||
          mat_type == CORE::Materials::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes,
      "The material parameter definition '%s' is not supported by Beam3eb element! "
      "Choose MAT_BeamKirchhoffTorsionFreeElastHyper or "
      "MAT_BeamKirchhoffTorsionFreeElastHyper_ByModes!",
      to_string(mat_type).data());

  return true;
}

FOUR_C_NAMESPACE_CLOSE
