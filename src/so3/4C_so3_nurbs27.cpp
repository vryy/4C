/*----------------------------------------------------------------------*/
/*! \file

\brief implementation of the quadratic NURBS 27 element

\level 2


*----------------------------------------------------------------------*/

#include "4C_so3_nurbs27.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_so3_surface.hpp"
#include "4C_so3_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Discret::ELEMENTS::Nurbs::SoNurbs27Type Discret::ELEMENTS::Nurbs::SoNurbs27Type::instance_;

Discret::ELEMENTS::Nurbs::SoNurbs27Type& Discret::ELEMENTS::Nurbs::SoNurbs27Type::Instance()
{
  return instance_;
}

Core::Communication::ParObject* Discret::ELEMENTS::Nurbs::SoNurbs27Type::Create(
    const std::vector<char>& data)
{
  auto* object = new Discret::ELEMENTS::Nurbs::SoNurbs27(-1, -1);
  object->unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::SoNurbs27Type::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == get_element_type_string())
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Nurbs::SoNurbs27(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::Nurbs::SoNurbs27Type::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Nurbs::SoNurbs27(id, owner));
  return ele;
}


void Discret::ELEMENTS::Nurbs::SoNurbs27Type::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::Nurbs::SoNurbs27Type::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return ComputeSolid3DNullSpace(node, x0);
}

void Discret::ELEMENTS::Nurbs::SoNurbs27Type::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions[get_element_type_string()];

  defs["NURBS27"] = Input::LineDefinition::Builder()
                        .add_int_vector("NURBS27", 27)
                        .add_named_int("MAT")
                        .add_named_int_vector("GP", 3)
                        .build();
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                                       |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::SoNurbs27::SoNurbs27(int id, int owner) : SoBase(id, owner)
{
  invJ_.resize(NUMGPT_SONURBS27, Core::LinAlg::Matrix<NUMDIM_SONURBS27, NUMDIM_SONURBS27>(true));
  detJ_.resize(NUMGPT_SONURBS27, 0.0);
  SetNurbsElement() = true;

  Teuchos::RCP<const Teuchos::ParameterList> params =
      Global::Problem::Instance()->getParameterList();
  if (params != Teuchos::null)
  {
    Discret::ELEMENTS::UTILS::ThrowErrorFDMaterialTangent(
        Global::Problem::Instance()->structural_dynamic_params(), get_element_type_string());
  }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                                  |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Nurbs::SoNurbs27::SoNurbs27(const Discret::ELEMENTS::Nurbs::SoNurbs27& old)
    : SoBase(old), detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (int i = 0; i < (int)invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }
  SetNurbsElement() = true;

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Nurbs::SoNurbs27::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::Nurbs::SoNurbs27(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Nurbs::SoNurbs27::Shape() const
{
  return Core::FE::CellType::nurbs27;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Nurbs::SoNurbs27::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  SoBase::pack(data);

  // detJ_
  add_to_pack(data, detJ_);

  // invJ_
  const auto size = (int)invJ_.size();
  add_to_pack(data, size);
  for (int i = 0; i < size; ++i) add_to_pack(data, invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Nurbs::SoNurbs27::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  SoBase::unpack(basedata);
  // detJ_
  extract_from_pack(position, data, detJ_);
  // invJ_
  int size = 0;
  extract_from_pack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<NUMDIM_SONURBS27, NUMDIM_SONURBS27>(true));
  for (int i = 0; i < size; ++i) extract_from_pack(position, data, invJ_[i]);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                                         |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Nurbs::SoNurbs27::Print(std::ostream& os) const
{
  os << "So_nurbs27 ";
  Element::Print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                                      |
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Nurbs::SoNurbs27::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, SoNurbs27>(
      Core::Communication::buildSurfaces, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                                        |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Nurbs::SoNurbs27::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, SoNurbs27>(
      Core::Communication::buildLines, *this);
}

FOUR_C_NAMESPACE_CLOSE
