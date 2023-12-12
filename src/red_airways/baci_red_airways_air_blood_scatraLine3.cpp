/*---------------------------------------------------------------------*/
/*! \file

\brief Incomplete! - Purpose: Implements RedAirBloodScatraLine3 element


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_red_airways_elementbase.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

using namespace CORE::DRT::UTILS;

DRT::ELEMENTS::RedAirBloodScatraLine3Type DRT::ELEMENTS::RedAirBloodScatraLine3Type::instance_;

DRT::ELEMENTS::RedAirBloodScatraLine3Type& DRT::ELEMENTS::RedAirBloodScatraLine3Type::Instance()
{
  return instance_;
}

CORE::COMM::ParObject* DRT::ELEMENTS::RedAirBloodScatraLine3Type::Create(
    const std::vector<char>& data)
{
  DRT::ELEMENTS::RedAirBloodScatraLine3* object = new DRT::ELEMENTS::RedAirBloodScatraLine3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAirBloodScatraLine3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_AIR_BLOOD_SCATRA_LINE3")
  {
    Teuchos::RCP<DRT::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::RedAirBloodScatraLine3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RedAirBloodScatraLine3Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele =
      Teuchos::rcp(new DRT::ELEMENTS::RedAirBloodScatraLine3(id, owner));
  return ele;
}


void DRT::ELEMENTS::RedAirBloodScatraLine3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs =
      definitions["RED_AIR_BLOOD_SCATRA_LINE3"];

  defs["LINE3"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE3", 3)
                      .AddNamedDouble("DiffusionCoefficient")
                      .AddNamedDouble("WallThickness")
                      .AddNamedDouble("PercentageOfDiffusionArea")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 05/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirBloodScatraLine3::RedAirBloodScatraLine3(int id, int owner)
    : DRT::Element(id, owner), data_()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 05/13|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirBloodScatraLine3::RedAirBloodScatraLine3(
    const DRT::ELEMENTS::RedAirBloodScatraLine3& old)
    : DRT::Element(old),
      data_(old.data_),
      elemParams_(old.elemParams_),
      generation_(old.generation_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAirBloodScatraLine3 and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::RedAirBloodScatraLine3::Clone() const
{
  DRT::ELEMENTS::RedAirBloodScatraLine3* newelement =
      new DRT::ELEMENTS::RedAirBloodScatraLine3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::RedAirBloodScatraLine3::Shape() const
{
  switch (NumNode())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      dserror("unexpected number of nodes %d", NumNode());
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirBloodScatraLine3::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);


  std::map<std::string, double>::const_iterator it;

  AddtoPack(data, (int)(elemParams_.size()));
  for (it = elemParams_.begin(); it != elemParams_.end(); it++)
  {
    AddtoPack(data, it->first);
    AddtoPack(data, it->second);
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirBloodScatraLine3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  std::map<std::string, double> it;
  int n = 0;

  ExtractfromPack(position, data, n);

  for (int i = 0; i < n; i++)
  {
    std::string name;
    double val;
    ExtractfromPack(position, data, name);
    ExtractfromPack(position, data, val);
    elemParams_[name] = val;
  }


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 05/13|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirBloodScatraLine3::Print(std::ostream& os) const
{
  os << "RedAirBloodScatraLine3 ";
  Element::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data                     ismail 05/13 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirBloodScatraLine3::VisNames(std::map<std::string, int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirBloodScatraLine3::VisData(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name, data)) return true;

  return false;
}



/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 04/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirBloodScatraLine3::getParams(std::string name, double& var)
{
  std::map<std::string, double>::iterator it;
  it = elemParams_.find(name);
  if (it == elemParams_.end())
  {
    dserror("[%s] is not found with in the element variables", name.c_str());
    exit(1);
  }
  var = elemParams_[name];
}

/*----------------------------------------------------------------------*
 |  Get element parameters (public)                        ismail 03/11 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirBloodScatraLine3::getParams(std::string name, int& var)
{
  //  if (name == "Generation")
  //  {
  //    var = generation_;
  //  }
  //  else
  {
    dserror("[%s] is not found with in the element INT variables", name.c_str());
    exit(1);
  }
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)              ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::RedAirBloodScatraLine3::Lines()
{
  dsassert(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

BACI_NAMESPACE_CLOSE
