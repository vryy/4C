/*----------------------------------------------------------------------------*/
/*! \file

\brief Input of 3D ALE element

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_ale_ale3.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_so3_material.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool Discret::ELEMENTS::Ale3::read_element(const std::string& eletype, const std::string& distype,
    const Core::IO::InputParameterContainer& container)
{
  // read number of material model
  int material_id = container.get<int>("MAT");
  set_material(0, Mat::Factory(material_id));

  Core::FE::CellType shape = Core::FE::StringToCellType(distype);

  std::cout << " distype " << distype << std::endl;

  Core::FE::GaussRule3D gaussrule;

  switch (shape)
  {
    case Core::FE::CellType::hex8:
    {
      gaussrule = Core::FE::GaussRule3D::hex_8point;
      break;
    }
    case Core::FE::CellType::hex20:
    case Core::FE::CellType::hex27:
    {
      gaussrule = Core::FE::GaussRule3D::hex_27point;
      break;
    }
    case Core::FE::CellType::pyramid5:
    {
      gaussrule = Core::FE::GaussRule3D::pyramid_8point;
      break;
    }
    case Core::FE::CellType::tet4:
    {
      gaussrule = Core::FE::GaussRule3D::tet_1point;
      break;
    }
    case Core::FE::CellType::tet10:
    {
      gaussrule = Core::FE::GaussRule3D::tet_4point;
      break;
    }
    default:
      FOUR_C_THROW("Unknown distype %s for ALE3 element", distype.c_str());
      // just set to something to shutup compiler
      gaussrule = Core::FE::GaussRule3D::undefined;
      break;
  }  // end switch distype

  // set up of materials with GP data (e.g., history variables)
  Teuchos::RCP<Mat::So3Material> so3mat = Teuchos::rcp_dynamic_cast<Mat::So3Material>(material());

  const Core::FE::IntegrationPoints3D intpoints(gaussrule);
  const int numgp = intpoints.nquad;
  so3mat->setup(numgp, container);

  return true;
}

FOUR_C_NAMESPACE_CLOSE
