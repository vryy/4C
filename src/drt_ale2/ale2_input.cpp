/*----------------------------------------------------------------------------*/
/*!
\file ale2_input.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289 15362
</pre>
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale2.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool DRT::ELEMENTS::Ale2::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // get gauss rule
  const DRT::UTILS::GaussRule2D gaussrule = getOptimalGaussrule(Shape());
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);
  const int numgp = intpoints.nquad;

  // get material
  Teuchos::RCP<MAT::Material> mat = Material();
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(mat, true);

  // call material setup
  so3mat->Setup(numgp, linedef);
  return true;
}
