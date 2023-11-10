/*----------------------------------------------------------------------------*/
/*! \file

\brief Input of 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_ale_ale3.H"
#include "baci_io_linedefinition.H"
#include "baci_mat_so3_material.H"

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

  CORE::FE::CellType shape = CORE::FE::StringToCellType(distype);

  std::cout << " distype " << distype << std::endl;

  CORE::DRT::UTILS::GaussRule3D gaussrule;

  switch (shape)
  {
    case CORE::FE::CellType::hex8:
    {
      gaussrule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
      break;
    }
    case CORE::FE::CellType::hex20:
    case CORE::FE::CellType::hex27:
    {
      gaussrule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
      break;
    }
    case CORE::FE::CellType::pyramid5:
    {
      gaussrule = CORE::DRT::UTILS::GaussRule3D::pyramid_8point;
      break;
    }
    case CORE::FE::CellType::tet4:
    {
      gaussrule = CORE::DRT::UTILS::GaussRule3D::tet_1point;
      break;
    }
    case CORE::FE::CellType::tet10:
    {
      gaussrule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
      break;
    }
    default:
      dserror("Unknown distype %s for ALE3 element", distype.c_str());
      // just set to something to shutup compiler
      gaussrule = CORE::DRT::UTILS::GaussRule3D::undefined;
      break;
  }  // end switch distype

  // set up of materials with GP data (e.g., history variables)
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());

  const CORE::DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  const int numgp = intpoints.nquad;
  so3mat->Setup(numgp, linedef);

  return true;
}
