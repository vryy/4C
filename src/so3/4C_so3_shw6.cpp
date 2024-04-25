/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/

#include "4C_so3_shw6.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_utils.hpp"
#include "4C_so3_weg6.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::SoShw6Type DRT::ELEMENTS::SoShw6Type::instance_;

DRT::ELEMENTS::SoShw6Type& DRT::ELEMENTS::SoShw6Type::Instance() { return instance_; }


CORE::COMM::ParObject* DRT::ELEMENTS::SoShw6Type::Create(const std::vector<char>& data)
{
  auto* object = new DRT::ELEMENTS::SoShw6(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoShw6Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == GetElementTypeString())
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoShw6(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::SoShw6Type::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::SoShw6(id, owner));
  return ele;
}


void DRT::ELEMENTS::SoShw6Type::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::SoShw6Type::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void DRT::ELEMENTS::SoShw6Type::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions[GetElementTypeString()];

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
DRT::ELEMENTS::SoShw6::SoShw6(int id, int owner) : DRT::ELEMENTS::SoWeg6(id, owner)
{
  eastype_ = soshw6_easnone;
  neas_ = 0;
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;

  Teuchos::RCP<const Teuchos::ParameterList> params =
      GLOBAL::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    DRT::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        GLOBAL::Problem::Instance()->StructuralDynamicParams(), GetElementTypeString());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoShw6::SoShw6(const DRT::ELEMENTS::SoShw6& old) : DRT::ELEMENTS::SoWeg6(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoShw6::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::SoShw6(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoShw6::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class So_weg6 Element
  DRT::ELEMENTS::SoWeg6::Pack(data);
  // eastype_
  AddtoPack(data, eastype_);
  // neas_
  AddtoPack(data, neas_);
  // easdata_
  PackEasData(data);
  // reordering
  AddtoPack(data, optimal_parameterspace_map_);
  AddtoPack(data, nodes_rearranged_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoShw6::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class So_weg6 Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  DRT::ELEMENTS::SoWeg6::Unpack(basedata);
  // eastype_
  eastype_ = static_cast<EASType>(ExtractInt(position, data));
  // neas_
  ExtractfromPack(position, data, neas_);
  // easdata_
  UnpackEasData(position, data);
  // reordering
  optimal_parameterspace_map_ = ExtractInt(position, data);
  nodes_rearranged_ = ExtractInt(position, data);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoShw6::Print(std::ostream& os) const
{
  os << "So_shw6 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

FOUR_C_NAMESPACE_CLOSE
