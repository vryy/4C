/*---------------------------------------------------------------------*/
/*! \file

\brief Implements an acinus element


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_io_linedefinition.hpp"
#include "baci_lib_discret.hpp"
#include "baci_mat_maxwell_0d_acinus.hpp"
#include "baci_red_airways_elementbase.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

using namespace CORE::FE;

DRT::ELEMENTS::RedAcinusType DRT::ELEMENTS::RedAcinusType::instance_;

DRT::ELEMENTS::RedAcinusType& DRT::ELEMENTS::RedAcinusType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::RedAcinusType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::RedAcinus* object = new DRT::ELEMENTS::RedAcinus(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAcinusType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_ACINUS")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAcinusType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::RedAcinus(id, owner));
  return ele;
}


/*--------------------------------------------------------------------  *
 | Read RED_ACINUS element line and add element specific parameters     |
 |                                                             (public) |
 |                                                           roth 10/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinusType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["RED_ACINUS"];

  defs["LINE2"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddNamedDouble("AcinusVolume")
                      .AddNamedDouble("AlveolarDuctVolume")
                      .AddOptionalNamedDouble("E1_0")
                      .AddOptionalNamedDouble("E1_LIN")
                      .AddOptionalNamedDouble("E1_EXP")
                      .AddOptionalNamedDouble("TAU")
                      .AddOptionalNamedDouble("E1_01")
                      .AddOptionalNamedDouble("E1_LIN1")
                      .AddOptionalNamedDouble("E1_EXP1")
                      .AddOptionalNamedDouble("TAU1")
                      .AddOptionalNamedDouble("E1_02")
                      .AddOptionalNamedDouble("E1_LIN2")
                      .AddOptionalNamedDouble("E1_EXP2")
                      .AddOptionalNamedDouble("TAU2")
                      .AddOptionalNamedDouble("KAPPA")
                      .AddOptionalNamedDouble("BETA")
                      .AddOptionalNamedDouble("Area")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(int id, int owner) : DRT::Element(id, owner) {}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAcinus::RedAcinus(const DRT::ELEMENTS::RedAcinus& old)
    : DRT::Element(old),
      elemType_(old.elemType_),
      resistance_(old.elemType_),
      acinusParams_(old.acinusParams_)
{
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAcinus and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::RedAcinus::Clone() const
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
  switch (NumNode())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
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

  AddtoPack(data, elemType_);
  AddtoPack(data, resistance_);

  AddtoPack(data, acinusParams_.volume_relaxed);
  AddtoPack(data, acinusParams_.alveolar_duct_volume);
  AddtoPack(data, acinusParams_.area);
  AddtoPack(data, acinusParams_.volume_init);
  AddtoPack(data, acinusParams_.generation);

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

  ExtractfromPack(position, data, elemType_);
  ExtractfromPack(position, data, resistance_);

  ExtractfromPack(position, data, acinusParams_.volume_relaxed);
  ExtractfromPack(position, data, acinusParams_.alveolar_duct_volume);
  ExtractfromPack(position, data, acinusParams_.area);
  ExtractfromPack(position, data, acinusParams_.volume_init);
  ExtractfromPack(position, data, acinusParams_.generation);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

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
std::vector<double> DRT::ELEMENTS::RedAcinus::ElementCenterRefeCoords()
{
  //  // update element geometry
  DRT::Node** nodes = Nodes();

  CORE::LINALG::SerialDenseMatrix mat(NumNode(), 3, false);
  for (int i = 0; i < NumNode(); ++i)
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
    for (int j = 0; j < NumNode(); ++j)
    {
      var = var + mat(j, i);
    }
    centercoords[i] = var / NumNode();
  }

  return centercoords;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAcinus::VisNames(std::map<std::string, int>& names)
{
  Teuchos::RCP<MAT::Material> mat = Material();

  // cast to specific material, because general material does not have VisNames/VisData
  Teuchos::RCP<MAT::Maxwell_0d_acinus> mxwll_0d_acin =
      Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus>(Material());
  mxwll_0d_acin->VisNames(names);
}


/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAcinus::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  // cast to specific material, because general material does not have VisNames/VisData
  Teuchos::RCP<MAT::Maxwell_0d_acinus> mxwll_0d_acin =
      Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus>(Material());

  return mxwll_0d_acin->VisData(name, data, this->Id());
}


void DRT::ELEMENTS::RedAcinus::UpdateRelaxedVolume(double newVol)
{
  acinusParams_.volume_relaxed = newVol;
}


const DRT::REDAIRWAYS::AcinusParams& DRT::ELEMENTS::RedAcinus::GetAcinusParams() const
{
  return acinusParams_;
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)              ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::RedAcinus::Lines()
{
  dsassert(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

BACI_NAMESPACE_CLOSE
