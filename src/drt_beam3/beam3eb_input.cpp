/*----------------------------------------------------------------------------*/
/*!
\file beam3eb_input.cpp

\brief three dimensional nonlinear torsionless rod based on a C1 curve

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam3eb.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_fem_general/largerotations.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3eb::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  if (Material()->Parameter()->Name() != "MAT_BeamKirchhoffTorsionFreeElastHyper" and
      Material()->Parameter()->Name() != "MAT_BeamKirchhoffTorsionFreeElastHyper_ByModes")
  {
    dserror(
        "The material parameter definition '%s' is not supported by "
        "Beam3eb element! "
        "Choose MAT_BeamKirchhoffTorsionFreeElastHyper or "
        "MAT_BeamKirchhoffTorsionFreeElastHyper_ByModes!",
        Material()->Parameter()->Name().c_str());
  }

  return true;
}
