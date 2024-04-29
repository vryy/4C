/*----------------------------------------------------------------------*/
/*! \file
\brief input for the artery element

\level 3


*----------------------------------------------------------------------*/

#include "4C_art_net_artery.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_cnst_1d_art.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Artery::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  int ngp;
  linedef->ExtractInt("GP", ngp);

  switch (ngp)
  {
    case 1:
      gaussrule_ = CORE::FE::GaussRule1D::line_1point;
      break;
    case 2:
      gaussrule_ = CORE::FE::GaussRule1D::line_2point;
      break;
    case 3:
      gaussrule_ = CORE::FE::GaussRule1D::line_3point;
      break;
    case 4:
      gaussrule_ = CORE::FE::GaussRule1D::line_4point;
      break;
    case 5:
      gaussrule_ = CORE::FE::GaussRule1D::line_5point;
      break;
    case 6:
      gaussrule_ = CORE::FE::GaussRule1D::line_6point;
      break;
    case 7:
      gaussrule_ = CORE::FE::GaussRule1D::line_7point;
      break;
    case 8:
      gaussrule_ = CORE::FE::GaussRule1D::line_8point;
      break;
    case 9:
      gaussrule_ = CORE::FE::GaussRule1D::line_9point;
      break;
    case 10:
      gaussrule_ = CORE::FE::GaussRule1D::line_10point;
      break;
    default:
      FOUR_C_THROW("Reading of ART element failed: Gaussrule for line not supported!\n");
  }

  // read artery implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = INPAR::ARTDYN::impltype_undefined;
  else if (impltype == "LinExp")
    impltype_ = INPAR::ARTDYN::impltype_lin_exp;
  else if (impltype == "PressureBased")
    impltype_ = INPAR::ARTDYN::impltype_pressure_based;
  else
    FOUR_C_THROW("Invalid implementation type for ARTERY elements!");

  // extract diameter
  double diam = 0.0;
  linedef->ExtractDouble("DIAM", diam);

  // set diameter in material
  SetDiamInMaterial(diam);

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Artery::SetDiamInMaterial(const double diam)
{
  // now the element knows its material, and we can use it to set the diameter
  Teuchos::RCP<MAT::Material> mat = Material();
  if (mat->MaterialType() == CORE::Materials::m_cnst_art)
  {
    MAT::Cnst1dArt* arterymat = dynamic_cast<MAT::Cnst1dArt*>(mat.get());
    arterymat->SetDiam(diam);
    arterymat->SetDiamInitial(diam);
  }
  else
    FOUR_C_THROW("Artery element got unsupported material type %d", mat->MaterialType());
  return;
}

FOUR_C_NAMESPACE_CLOSE
