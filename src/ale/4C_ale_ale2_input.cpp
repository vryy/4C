/*----------------------------------------------------------------------------*/
/*! \file

\brief Input of 2D ALE elements

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale2.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool Discret::ELEMENTS::Ale2::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // get gauss rule
  const Core::FE::GaussRule2D gaussrule = get_optimal_gaussrule(Shape());
  const Core::FE::IntegrationPoints2D intpoints(gaussrule);
  const int numgp = intpoints.nquad;

  // get material
  Teuchos::RCP<Core::Mat::Material> mat = Material();
  Teuchos::RCP<Mat::So3Material> so3mat = Teuchos::rcp_dynamic_cast<Mat::So3Material>(mat, true);

  // call material setup
  so3mat->Setup(numgp, linedef);
  return true;
}

FOUR_C_NAMESPACE_CLOSE
