/*----------------------------------------------------------------------------*/
/*! \file
\brief A Kirchhoff-Love shell element

\level 1

*/
/*---------------------------------------------------------------------------*/

#include "4C_shell_kl_nurbs.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_so3_nullspace.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
Discret::ELEMENTS::KirchhoffLoveShellNurbsType
    Discret::ELEMENTS::KirchhoffLoveShellNurbsType::instance_;


/**
 *
 */
Discret::ELEMENTS::KirchhoffLoveShellNurbsType&
Discret::ELEMENTS::KirchhoffLoveShellNurbsType::instance()
{
  return instance_;
}


/**
 *
 */
Core::Communication::ParObject* Discret::ELEMENTS::KirchhoffLoveShellNurbsType::create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::KirchhoffLoveShellNurbs* object =
      new Discret::ELEMENTS::KirchhoffLoveShellNurbs(-1, -1);
  object->unpack(data);
  return object;
}


/**
 *
 */
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::KirchhoffLoveShellNurbsType::create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "SHELL_KIRCHHOFF_LOVE_NURBS" and eledistype == "NURBS9")
  {
    return Teuchos::rcp(new Discret::ELEMENTS::KirchhoffLoveShellNurbs(id, owner));
  }
  return Teuchos::null;
}

/**
 *
 */
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::KirchhoffLoveShellNurbsType::create(
    const int id, const int owner)
{
  return Teuchos::rcp(new Discret::ELEMENTS::KirchhoffLoveShellNurbs(id, owner));
}

/**
 *
 */
void Discret::ELEMENTS::KirchhoffLoveShellNurbsType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  FOUR_C_THROW("NodalBlockInformation not implemented");
}

/**
 *
 */
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::KirchhoffLoveShellNurbsType::compute_null_space(
    Core::Nodes::Node& node, const double* x0, int const numdof, int const dimnsp)
{
  FOUR_C_THROW("ComputeNullSpace not implemented");
}

/**
 *
 */
void Discret::ELEMENTS::KirchhoffLoveShellNurbsType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["SHELL_KIRCHHOFF_LOVE_NURBS"];

  defs["NURBS9"] = Input::LineDefinition::Builder()
                       .add_int_vector("NURBS9", 9)
                       .add_named_int("MAT")
                       .add_named_int_vector("GP", 2)
                       .build();
}


/**
 *
 */
Discret::ELEMENTS::KirchhoffLoveShellNurbs::KirchhoffLoveShellNurbs(int id, int owner)
    : Core::Elements::Element(id, owner),
      material_(0),
      gaussrule_({Core::FE::GaussRule1D::undefined, Core::FE::GaussRule1D::undefined})
{
}

/**
 *
 */
Discret::ELEMENTS::KirchhoffLoveShellNurbs::KirchhoffLoveShellNurbs(
    const Discret::ELEMENTS::KirchhoffLoveShellNurbs& old)
    : Core::Elements::Element(old), material_(old.material_), gaussrule_(old.gaussrule_)
{
}

/**
 *
 */
Core::Elements::Element* Discret::ELEMENTS::KirchhoffLoveShellNurbs::clone() const
{
  return new Discret::ELEMENTS::KirchhoffLoveShellNurbs(*this);
}

/**
 *
 */
void Discret::ELEMENTS::KirchhoffLoveShellNurbs::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // add base class Element
  Core::Elements::Element::pack(data);
  // material_
  add_to_pack(data, material_);
  // gaussrule_
  add_to_pack(data, gaussrule_[0]);
  add_to_pack(data, gaussrule_[1]);
}

/**
 *
 */
void Discret::ELEMENTS::KirchhoffLoveShellNurbs::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Core::Elements::Element::unpack(basedata);
  // material_
  extract_from_pack(position, data, material_);
  // gaussrule_
  extract_from_pack(position, data, gaussrule_[0]);
  extract_from_pack(position, data, gaussrule_[1]);
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/**
 *
 */
void Discret::ELEMENTS::KirchhoffLoveShellNurbs::set_params_interface_ptr(
    const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<Solid::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/**
 *
 */
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::KirchhoffLoveShellNurbs::surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
