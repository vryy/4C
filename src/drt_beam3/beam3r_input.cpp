/*-----------------------------------------------------------------------------------------------*/
/*!

\brief input related methods of 3D nonlinear Reissner beam element

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam3r.H"

#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_fem_general/largerotations.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3r::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  /* the triad field is discretized with Lagrange polynomials of order NumNode()-1;
   * the centerline is either discretized in the same way or with 3rd order Hermite polynomials;
   * we thus make a difference between nnodetriad and nnodecl;
   * assumptions: nnodecl<=nnodetriad
   * first nodes with local ID 0...nnodecl-1 are used for interpolation of centerline AND triad
   * field*/
  const int nnodetriad = NumNode();


  // read number of material model and cross-section specs
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  if (Material()->Parameter()->Name() != "MAT_BeamReissnerElastHyper" and
      Material()->Parameter()->Name() != "MAT_BeamReissnerElastHyper_ByModes")
  {
    dserror(
        "The material parameter definition '%s' is not supported by Beam3r element! "
        "Choose MAT_BeamReissnerElastHyper or MAT_BeamReissnerElastHyper_ByModes!",
        Material()->Parameter()->Name().c_str());
  }


  if (linedef->HaveNamed("HERM2LIN2") or linedef->HaveNamed("HERM2LIN3") or
      linedef->HaveNamed("HERM2LIN4") or linedef->HaveNamed("HERM2LIN5"))
    centerline_hermite_ = true;
  else
    centerline_hermite_ = false;

  // read whether automatic differentiation via Sacado::Fad package shall be used
  useFAD_ = linedef->HaveNamed("FAD") ? true : false;


  // store nodal triads according to input file
  theta0node_.resize(nnodetriad);

  /* Attention! expression "TRIADS" in input file is misleading.
   * The 3 specified values per node define a rotational pseudovector, which
   * parameterizes the orientation of the triad at this node
   * (relative to the global reference coordinate system)*/
  /* extract rotational pseudovectors at element nodes in reference configuration
   *  and save them as quaternions at each node, respectively*/
  std::vector<double> nodal_rotvecs;
  linedef->ExtractDoubleVector("TRIADS", nodal_rotvecs);

  for (int node = 0; node < nnodetriad; node++)
    for (int dim = 0; dim < 3; dim++) theta0node_[node](dim) = nodal_rotvecs[3 * node + dim];

  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3r::SetCenterlineHermite(const bool yesno) { centerline_hermite_ = yesno; }
