/*---------------------------------------------------------------------*/
/*! \file

\brief Implements an acinus element


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

using namespace CORE::FE;

DRT::ELEMENTS::RedAcinusType DRT::ELEMENTS::RedAcinusType::instance_;

DRT::ELEMENTS::RedAcinusType& DRT::ELEMENTS::RedAcinusType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::RedAcinusType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::RedAcinus* object = new DRT::ELEMENTS::RedAcinus(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::RedAcinusType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_ACINUS")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::RedAcinusType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id, owner));
  return ele;
}


/*--------------------------------------------------------------------  *
 | Read RED_ACINUS element line and add element specific parameters     |
 |                                                             (public) |
 |                                                           roth 10/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinusType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["RED_ACINUS"];

  defs["LINE2"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddNamedDouble("AcinusVolume")
                      .AddNamedDouble("AlveolarDuctVolume")
                      .add_optional_named_double("E1_0")
                      .add_optional_named_double("E1_LIN")
                      .add_optional_named_double("E1_EXP")
                      .add_optional_named_double("TAU")
                      .add_optional_named_double("E1_01")
                      .add_optional_named_double("E1_LIN1")
                      .add_optional_named_double("E1_EXP1")
                      .add_optional_named_double("TAU1")
                      .add_optional_named_double("E1_02")
                      .add_optional_named_double("E1_LIN2")
                      .add_optional_named_double("E1_EXP2")
                      .add_optional_named_double("TAU2")
                      .add_optional_named_double("KAPPA")
                      .add_optional_named_double("BETA")
                      .add_optional_named_double("Area")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(int id, int owner) : CORE::Elements::Element(id, owner) {}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(const DRT::ELEMENTS::RedAcinus& old)
    : CORE::Elements::Element(old),
      elem_type_(old.elem_type_),
      resistance_(old.elem_type_),
      acinus_params_(old.acinus_params_)
{
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAcinus and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::RedAcinus::Clone() const
{
  DRT::ELEMENTS::RedAcinus* newelement = new DRT::ELEMENTS::RedAcinus(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::RedAcinus::Shape() const
{
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);

  AddtoPack(data, elem_type_);
  AddtoPack(data, resistance_);

  AddtoPack(data, acinus_params_.volume_relaxed);
  AddtoPack(data, acinus_params_.alveolar_duct_volume);
  AddtoPack(data, acinus_params_.area);
  AddtoPack(data, acinus_params_.volume_init);
  AddtoPack(data, acinus_params_.generation);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  ExtractfromPack(position, data, elem_type_);
  ExtractfromPack(position, data, resistance_);

  ExtractfromPack(position, data, acinus_params_.volume_relaxed);
  ExtractfromPack(position, data, acinus_params_.alveolar_duct_volume);
  ExtractfromPack(position, data, acinus_params_.area);
  ExtractfromPack(position, data, acinus_params_.volume_init);
  ExtractfromPack(position, data, acinus_params_.generation);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::Print(std::ostream& os) const
{
  os << "RedAcinus ";
  Element::Print(os);

  return;
}

/*-----------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::RedAcinus::element_center_refe_coords()
{
  //  // update element geometry
  CORE::Nodes::Node** nodes = Nodes();

  CORE::LINALG::SerialDenseMatrix mat(num_node(), 3, false);
  for (int i = 0; i < num_node(); ++i)
  {
    const auto& x = nodes[i]->X();
    mat(i, 0) = x[0];
    mat(i, 1) = x[1];
    mat(i, 2) = x[2];
  }

  std::vector<double> centercoords(3, 0);
  for (int i = 0; i < 3; ++i)
  {
    double var = 0;
    for (int j = 0; j < num_node(); ++j)
    {
      var = var + mat(j, i);
    }
    centercoords[i] = var / num_node();
  }

  return centercoords;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::VisNames(std::map<std::string, int>& names)
{
  Teuchos::RCP<CORE::MAT::Material> mat = Material();

  // cast to specific material, because general material does not have VisNames/VisData
  Teuchos::RCP<MAT::Maxwell0dAcinus> mxwll_0d_acin =
      Teuchos::rcp_dynamic_cast<MAT::Maxwell0dAcinus>(Material());
  mxwll_0d_acin->VisNames(names);
}


/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAcinus::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (CORE::Elements::Element::VisData(name, data)) return true;

  // cast to specific material, because general material does not have VisNames/VisData
  Teuchos::RCP<MAT::Maxwell0dAcinus> mxwll_0d_acin =
      Teuchos::rcp_dynamic_cast<MAT::Maxwell0dAcinus>(Material());

  return mxwll_0d_acin->VisData(name, data, this->Id());
}


void DRT::ELEMENTS::RedAcinus::UpdateRelaxedVolume(double newVol)
{
  acinus_params_.volume_relaxed = newVol;
}


const DRT::REDAIRWAYS::AcinusParams& DRT::ELEMENTS::RedAcinus::GetAcinusParams() const
{
  return acinus_params_;
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)              ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::RedAcinus::Lines()
{
  FOUR_C_ASSERT(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
