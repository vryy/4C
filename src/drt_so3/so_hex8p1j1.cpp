/*----------------------------------------------------------------------*/
/*! \file
\brief 'Q1P0' element in 8-node hexahedron shape

\level 2

*/
/*----------------------------------------------------------------------*/

#include "so_hex8p1j1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "so_hex8.H"
#include "../drt_lib/drt_linedefinition.H"

DRT::ELEMENTS::So_Hex8P1J1Type DRT::ELEMENTS::So_Hex8P1J1Type::instance_;

DRT::ELEMENTS::So_Hex8P1J1Type& DRT::ELEMENTS::So_Hex8P1J1Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_Hex8P1J1Type::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::So_Hex8P1J1* object = new DRT::ELEMENTS::So_Hex8P1J1(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_Hex8P1J1Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SOLIDH8P1J1")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_Hex8P1J1(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_Hex8P1J1Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_Hex8P1J1(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_Hex8P1J1Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  //   numdf = 3;
  //   dimns = 6;
  //   nv = 3;
}

void DRT::ELEMENTS::So_Hex8P1J1Type::ComputeNullSpace(
    DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
{
  //   DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_Hex8P1J1Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8P1J1"];

  defs["HEX8"].AddIntVector("HEX8", 8).AddNamedInt("MAT").AddNamedString("KINEM");
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                               lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_Hex8P1J1::So_Hex8P1J1(int id, int owner) : DRT::ELEMENTS::So_hex8(id, owner)
{
  K_pu_.PutScalar(0.0);
  K_tu_.PutScalar(0.0);

  R_t_.PutScalar(0.0);
  R_p_.PutScalar(0.0);

  K_tt_ = 0.0;
  K_pt_ = 0.0;

  p_.PutScalar(0.0);
  p_o_.PutScalar(0.0);

  t_.PutScalar(1.0);
  t_o_.PutScalar(1.0);

  m_.PutScalar(0.0);
  for (int i = 0; i < 3; ++i)
  {
    m_(i, 0) = 1.0;
  }

  Identity6_.PutScalar(0.0);
  for (int i = 0; i < 6; ++i)
  {
    Identity6_(i, i) = 1.0;
  }

  I_d_ = Identity6_;
  I_d_.MultiplyNT(-1.0 / 3.0, m_, m_, 1.0);

  I_0_.PutScalar(0.0);

  for (int i = 0; i < 3; ++i)
  {
    I_0_(i, i) = 1.0;
  }
  for (int i = 3; i < 6; ++i)
  {
    I_0_(i, i) = 0.5;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                          lw 12/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_Hex8P1J1::So_Hex8P1J1(const DRT::ELEMENTS::So_Hex8P1J1& old)
    : DRT::ELEMENTS::So_hex8(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_Hex8P1J1::Clone() const
{
  DRT::ELEMENTS::So_Hex8P1J1* newelement = new DRT::ELEMENTS::So_Hex8P1J1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_hex8 Element
  DRT::ELEMENTS::So_hex8::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                              lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::So_hex8::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                               lw 12/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_Hex8P1J1::~So_Hex8P1J1() { return; }


/*----------------------------------------------------------------------*
 |  print this element (public)                                 lw 12/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_Hex8P1J1::Print(std::ostream& os) const
{
  os << "So_Hex8P1J1 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}
