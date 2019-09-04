/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam3k.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_fem_general/largerotations.H"


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3k::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model and cross-sections specs
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  if (Material()->Parameter()->Name() != "MAT_BeamKirchhoffElastHyper" and
      Material()->Parameter()->Name() != "MAT_BeamKirchhoffElastHyper_ByModes")
  {
    dserror(
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
    dserror(
        "The variable ROTVEC can only take on the values 0 (tangent vectors as nodal DoFs) and "
        "1 (rotation vectors as nodal DoFs)!");

  if (wk == 0)
    weakkirchhoff_ = false;
  else if (wk == 1)
  {
    weakkirchhoff_ = true;
#ifdef CONSISTENTSPINSK
    dserror(
        "The flag CONSISTENTSPINSK is only possible for strong Kirchhoff constraint enforcement "
        "(weakkirchhoff_=false)");
#endif
  }
  else
    dserror(
        "The variable WK can only take on the values 0 (Kirchhoff constraint enforced in a strong "
        "manner) and "
        "1 (Kirchhoff constraint enforced in a weak manner)!");


  // extract triads at element nodes in reference configuration as rotation vectors and save them as
  // quaternions at each node, respectively
  std::vector<double> nodal_thetas;
  linedef->ExtractDoubleVector("TRIADS", nodal_thetas);
  theta0_.resize(BEAM3K_COLLOCATION_POINTS);
  for (int i = 0; i < BEAM3K_COLLOCATION_POINTS; i++)
  {
    for (int j = 0; j < 3; j++) (theta0_[i])(j) = nodal_thetas[3 * i + j];

    // Shift angles by 2PI in case these angles are not in the interval [-PI,PI].
    if (theta0_[i].Norm2() > M_PI)
    {
      LINALG::Matrix<4, 1> Q(true);
      LARGEROTATIONS::angletoquaternion(theta0_[i], Q);
      LARGEROTATIONS::quaterniontoangle(Q, theta0_[i]);
    }
  }

  // read whether automatic differentiation via Sacado::Fad package shall be used
  useFAD_ = linedef->HaveNamed("FAD") ? true : false;

  return true;
}
