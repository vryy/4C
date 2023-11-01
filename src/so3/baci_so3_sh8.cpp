/*----------------------------------------------------------------------*/
/*! \file

\brief solid shell8 element formulation

\level 1


*----------------------------------------------------------------------*/

#include "baci_so3_sh8.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_so3_hex8.H"
#include "baci_so3_nullspace.H"
#include "baci_so3_utils.H"
#include "baci_utils_exceptions.H"


DRT::ELEMENTS::So_sh8Type DRT::ELEMENTS::So_sh8Type::instance_;

DRT::ELEMENTS::So_sh8Type& DRT::ELEMENTS::So_sh8Type::Instance() { return instance_; }

DRT::ParObject* DRT::ELEMENTS::So_sh8Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So_sh8(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh8Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh8(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh8Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh8(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_sh8Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::So_sh8Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::So_sh8Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["HEX8"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("KINEM")
                     .AddNamedString("EAS")
                     .AddNamedString("ANS")
                     .AddNamedString("THICKDIR")
                     .AddOptionalNamedDoubleVector("RAD", 3)
                     .AddOptionalNamedDoubleVector("AXI", 3)
                     .AddOptionalNamedDoubleVector("CIR", 3)
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .AddOptionalNamedDoubleVector("FIBER2", 3)
                     .AddOptionalNamedDoubleVector("FIBER3", 3)
                     .AddOptionalNamedDouble("STRENGTH")
                     .AddOptionalNamedDouble("GROWTHTRIG")
                     .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::So_sh8(int id, int owner)
    : DRT::ELEMENTS::So_hex8(id, owner),
      thickdir_(undefined),
      anstype_(ansnone),
      nodes_rearranged_(false),
      thickvec_(3, 0.0)
{
  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        DRT::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::So_sh8(const DRT::ELEMENTS::So_sh8& old)
    : DRT::ELEMENTS::So_hex8(old),
      thickdir_(old.thickdir_),
      anstype_(old.anstype_),
      nodes_rearranged_(old.nodes_rearranged_),
      thickvec_(old.thickvec_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh8::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So_sh8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_hex8 Element
  DRT::ELEMENTS::So_hex8::Pack(data);
  // thickdir
  AddtoPack(data, thickdir_);
  AddtoPack(data, thickvec_);
  AddtoPack(data, anstype_);
  AddtoPack(data, nodes_rearranged_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_hex8 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::So_hex8::Unpack(basedata);
  // thickdir
  thickdir_ = static_cast<ThicknessDirection>(ExtractInt(position, data));
  ExtractfromPack(position, data, thickvec_);
  anstype_ = static_cast<ANSType>(ExtractInt(position, data));
  nodes_rearranged_ = ExtractInt(position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8::Print(std::ostream& os) const
{
  os << "So_sh8 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}
