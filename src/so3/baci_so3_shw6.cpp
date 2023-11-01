/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/

#include "baci_so3_shw6.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_globalproblem.H"
#include "baci_so3_nullspace.H"
#include "baci_so3_utils.H"
#include "baci_so3_weg6.H"


DRT::ELEMENTS::So_shw6Type DRT::ELEMENTS::So_shw6Type::instance_;

DRT::ELEMENTS::So_shw6Type& DRT::ELEMENTS::So_shw6Type::Instance() { return instance_; }


DRT::ParObject* DRT::ELEMENTS::So_shw6Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::So_shw6(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_shw6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_shw6(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_shw6Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_shw6(id, owner));
  return ele;
}


void DRT::ELEMENTS::So_shw6Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::So_shw6Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::So_shw6Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

  defs["WEDGE6"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("WEDGE6", 6)
                       .AddNamedInt("MAT")
                       .AddNamedString("KINEM")
                       .AddNamedString("EAS")
                       .AddOptionalTag("OPTORDER")
                       .AddOptionalNamedDoubleVector("RAD", 3)
                       .AddOptionalNamedDoubleVector("AXI", 3)
                       .AddOptionalNamedDoubleVector("CIR", 3)
                       .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_shw6::So_shw6(int id, int owner) : DRT::ELEMENTS::So_weg6(id, owner)
{
  eastype_ = soshw6_easnone;
  neas_ = 0;
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;

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
DRT::ELEMENTS::So_shw6::So_shw6(const DRT::ELEMENTS::So_shw6& old) : DRT::ELEMENTS::So_weg6(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_shw6::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So_shw6(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_shw6::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_weg6 Element
  DRT::ELEMENTS::So_weg6::Pack(data);
  // eastype_
  AddtoPack(data, eastype_);
  // neas_
  AddtoPack(data, neas_);
  // reordering
  AddtoPack(data, optimal_parameterspace_map_);
  AddtoPack(data, nodes_rearranged_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_shw6::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_weg6 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::So_weg6::Unpack(basedata);
  // eastype_
  eastype_ = static_cast<EASType>(ExtractInt(position, data));
  // neas_
  ExtractfromPack(position, data, neas_);
  // reordering
  optimal_parameterspace_map_ = ExtractInt(position, data);
  nodes_rearranged_ = ExtractInt(position, data);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_shw6::Print(std::ostream& os) const
{
  os << "So_shw6 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}
