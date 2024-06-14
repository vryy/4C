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
bool Discret::ELEMENTS::Artery::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  int ngp;
  linedef->extract_int("GP", ngp);

  switch (ngp)
  {
    case 1:
      gaussrule_ = Core::FE::GaussRule1D::line_1point;
      break;
    case 2:
      gaussrule_ = Core::FE::GaussRule1D::line_2point;
      break;
    case 3:
      gaussrule_ = Core::FE::GaussRule1D::line_3point;
      break;
    case 4:
      gaussrule_ = Core::FE::GaussRule1D::line_4point;
      break;
    case 5:
      gaussrule_ = Core::FE::GaussRule1D::line_5point;
      break;
    case 6:
      gaussrule_ = Core::FE::GaussRule1D::line_6point;
      break;
    case 7:
      gaussrule_ = Core::FE::GaussRule1D::line_7point;
      break;
    case 8:
      gaussrule_ = Core::FE::GaussRule1D::line_8point;
      break;
    case 9:
      gaussrule_ = Core::FE::GaussRule1D::line_9point;
      break;
    case 10:
      gaussrule_ = Core::FE::GaussRule1D::line_10point;
      break;
    default:
      FOUR_C_THROW("Reading of ART element failed: Gaussrule for line not supported!\n");
  }

  // read artery implementation type
  std::string impltype;
  linedef->extract_string("TYPE", impltype);

  if (impltype == "Undefined")
    impltype_ = Inpar::ArtDyn::impltype_undefined;
  else if (impltype == "LinExp")
    impltype_ = Inpar::ArtDyn::impltype_lin_exp;
  else if (impltype == "PressureBased")
    impltype_ = Inpar::ArtDyn::impltype_pressure_based;
  else
    FOUR_C_THROW("Invalid implementation type for ARTERY elements!");

  // extract diameter
  double diam = 0.0;
  linedef->extract_double("DIAM", diam);

  // set diameter in material
  SetDiamInMaterial(diam);

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Discret::ELEMENTS::Artery::SetDiamInMaterial(const double diam)
{
  // now the element knows its material, and we can use it to set the diameter
  Teuchos::RCP<Core::Mat::Material> mat = Material();
  if (mat->MaterialType() == Core::Materials::m_cnst_art)
  {
    Mat::Cnst1dArt* arterymat = dynamic_cast<Mat::Cnst1dArt*>(mat.get());
    arterymat->SetDiam(diam);
    arterymat->SetDiamInitial(diam);
  }
  else
    FOUR_C_THROW("Artery element got unsupported material type %d", mat->MaterialType());
  return;
}

FOUR_C_NAMESPACE_CLOSE
