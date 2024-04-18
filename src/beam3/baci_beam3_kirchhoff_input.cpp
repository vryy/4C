/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beam3_kirchhoff.hpp"
#include "baci_discretization_fem_general_largerotations.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_par_parameter.hpp"

FOUR_C_NAMESPACE_OPEN


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3k::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model and cross-sections specs
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  if (Material()->Parameter()->Name() != "MAT_BeamKirchhoffElastHyper" and
      Material()->Parameter()->Name() != "MAT_BeamKirchhoffElastHyper_ByModes")
  {
    FOUR_C_THROW(
        "The material parameter definition '%s' is not supported by Beam3k element! "
        "Choose MAT_BeamKirchhoffElastHyper or MAT_BeamKirchhoffElastHyper_ByModes!",
        Material()->Parameter()->Name().c_str());
  }


  int rotvec = 0;
  linedef->ExtractInt("ROTVEC", rotvec);

  int wk = 0;
  linedef->ExtractInt("WK", wk);

  if (rotvec == 0)
    rotvec_ = false;
  else if (rotvec == 1)
    rotvec_ = true;
  else
    FOUR_C_THROW(
        "The variable ROTVEC can only take on the values 0 (tangent vectors as nodal DoFs) and "
        "1 (rotation vectors as nodal DoFs)!");

  if (wk == 0)
    weakkirchhoff_ = false;
  else if (wk == 1)
  {
    weakkirchhoff_ = true;
#ifdef CONSISTENTSPINSK
    FOUR_C_THROW(
        "The flag CONSISTENTSPINSK is only possible for strong Kirchhoff constraint enforcement "
        "(weakkirchhoff_=false)");
#endif
  }
  else
    FOUR_C_THROW(
        "The variable WK can only take on the values 0 (Kirchhoff constraint enforced in a strong "
        "manner) and "
        "1 (Kirchhoff constraint enforced in a weak manner)!");


  // extract triads at element nodes in reference configuration as rotation vectors and save them as
  // quaternions at each node, respectively
  std::vector<double> nodal_thetas;
  linedef->ExtractDoubleVector("TRIADS", nodal_thetas);
  this->SetUpInitialRotations(nodal_thetas);

  // read whether automatic differentiation via Sacado::Fad package shall be used
  use_fad_ = linedef->HaveNamed("FAD") ? true : false;

  return true;
}

FOUR_C_NAMESPACE_CLOSE
