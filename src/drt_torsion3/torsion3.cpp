/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional torsion spring element

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "torsion3.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"

DRT::ELEMENTS::Torsion3Type DRT::ELEMENTS::Torsion3Type::instance_;

DRT::ELEMENTS::Torsion3Type& DRT::ELEMENTS::Torsion3Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::Torsion3Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Torsion3* object = new DRT::ELEMENTS::Torsion3(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Torsion3Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TORSION3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Torsion3(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Torsion3Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Torsion3(id, owner));
  return ele;
}


void DRT::ELEMENTS::Torsion3Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
}

void DRT::ELEMENTS::Torsion3Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  DRT::UTILS::ComputeStructure3DNullSpace(dis, ns, x0, numdf, dimns);
}

void DRT::ELEMENTS::Torsion3Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["TORSION3"];

  defs["LINE3"].AddIntVector("LINE3", 3).AddNamedInt("MAT").AddNamedString("BENDINGPOTENTIAL");
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(int id, int owner) : DRT::Element(id, owner), data_() { return; }
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::Torsion3(const DRT::ELEMENTS::Torsion3& old)
    : DRT::Element(old), data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 | Deep copy this instance of Torsion3 and return pointer to it (public)|
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Torsion3::Clone() const
{
  DRT::ELEMENTS::Torsion3* newelement = new DRT::ELEMENTS::Torsion3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Torsion3::~Torsion3() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                              cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Print(std::ostream& os) const { return; }


/*----------------------------------------------------------------------*
 |(public)                                                   cyron 02/10|
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Torsion3::Shape() const { return line3; }


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  Element::Pack(data);
  AddtoPack(data, bendingpotential_);
  AddtoPack(data, data_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           cyron 02/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);
  bendingpotential_ = static_cast<BendingPotential>(ExtractInt(position, data));
  std::vector<char> tmp(0);
  ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             cyron 02/10|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Torsion3::Lines()
{
  std::vector<Teuchos::RCP<Element>> lines(1);
  lines[0] = Teuchos::rcp(this, false);
  return lines;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Torsion3::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Torsion3::ParamsInterfacePtr()
{
  return interface_ptr_;
}
