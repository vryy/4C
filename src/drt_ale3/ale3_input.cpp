/*----------------------------------------------------------------------------*/
/*! \file

\brief Input of 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale3.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool DRT::ELEMENTS::Ale3::ReadElement(
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  linedef->Print(std::cout);

  DiscretizationType shape = StringToDistype(distype);

  std::cout << " distype " << distype << std::endl;

  DRT::UTILS::GaussRule3D gaussrule;

  switch (shape)
  {
    case hex8:
    {
      gaussrule = DRT::UTILS::intrule_hex_8point;
      break;
    }
    case hex20:
    case hex27:
    {
      gaussrule = DRT::UTILS::intrule_hex_27point;
      break;
    }
    case pyramid5:
    {
      gaussrule = DRT::UTILS::intrule_pyramid_8point;
      break;
    }
    case tet4:
    {
      gaussrule = DRT::UTILS::intrule_tet_1point;
      break;
    }
    case tet10:
    {
      gaussrule = DRT::UTILS::intrule_tet_4point;
      break;
    }
    default:
      dserror("Unknown distype %s for ALE3 element", distype.c_str());
      // just set to something to shutup compiler
      gaussrule = DRT::UTILS::intrule3D_undefined;
      break;
  }  // end switch distype

  // set up of materials with GP data (e.g., history variables)
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());

  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  const int numgp = intpoints.nquad;
  so3mat->Setup(numgp, linedef);

  return true;
}
