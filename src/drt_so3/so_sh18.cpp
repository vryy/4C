/*----------------------------------------------------------------------*/
/*! \file

\brief ToDo Add meaningful comment.

\level 1

\maintainer Christoph Meier

*----------------------------------------------------------------------*/

#include "so_sh18.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_mat/so3_material.H"

DRT::ELEMENTS::So_sh18Type DRT::ELEMENTS::So_sh18Type::instance_;

DRT::ELEMENTS::So_sh18Type& DRT::ELEMENTS::So_sh18Type::Instance() { return instance_; }
namespace
{
  const std::string name = DRT::ELEMENTS::So_sh18Type::Instance().Name();
}

DRT::ParObject* DRT::ELEMENTS::So_sh18Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_sh18* object = new DRT::ELEMENTS::So_sh18(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh18Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDSH18")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh18(id, owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh18Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh18(id, owner));
  return ele;
}

void DRT::ELEMENTS::So_sh18Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDSH18"];

  defs["HEX18"]
      .AddIntVector("HEX18", 18)
      .AddNamedInt("MAT")
      .AddNamedString("KINEM")
      .AddNamedString("TSL")
      .AddNamedString("MEL")
      .AddNamedString("CTL")
      .AddNamedString("VOL")
      .AddOptionalNamedDoubleVector("RAD", 3)
      .AddOptionalNamedDoubleVector("AXI", 3)
      .AddOptionalNamedDoubleVector("CIR", 3)
      .AddOptionalNamedDoubleVector("FIBER1", 3)
      .AddOptionalNamedDoubleVector("FIBER2", 3)
      .AddOptionalNamedDoubleVector("FIBER3", 3)
      .AddOptionalNamedDouble("STRENGTH")
      .AddOptionalNamedDouble("HU")
      .AddOptionalNamedDouble("lambda");
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh18::So_sh18(int id, int owner) : So_base(id, owner), So_hex18(id, owner)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh18::So_sh18(const DRT::ELEMENTS::So_sh18& old)
    : So_base(old),
      So_hex18(old),
      dsg_shear_(old.dsg_shear_),
      dsg_membrane_(old.dsg_membrane_),
      dsg_ctl_(old.dsg_ctl_),
      eas_(old.eas_)
{
  SetupDSG();
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh18::Clone() const
{
  DRT::ELEMENTS::So_sh18* newelement = new DRT::ELEMENTS::So_sh18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  So_base::Pack(data);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  // element technology bools
  AddtoPack(data, (int)dsg_shear_);
  AddtoPack(data, (int)dsg_membrane_);
  AddtoPack(data, (int)dsg_ctl_);
  AddtoPack(data, (int)eas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  So_base::Unpack(basedata);

  // detJ_
  ExtractfromPack(position, data, detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOH18, NUMDIM_SOH18>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  // element technology bools
  dsg_shear_ = ExtractInt(position, data);
  dsg_membrane_ = ExtractInt(position, data);
  dsg_ctl_ = ExtractInt(position, data);
  eas_ = ExtractInt(position, data);
  SetupDSG();

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh18::~So_sh18() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                             seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Print(std::ostream& os) const
{
  os << "So_sh18 ";
  Element::Print(os);
  return;
}
