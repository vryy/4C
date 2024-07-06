/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element

*----------------------------------------------------------------------*/
#include "4C_membrane.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneLine2Type::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new MembraneLine( id, owner ) );
  return Teuchos::null;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::MembraneLine3Type::create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new MembraneLine( id, owner ) );
  return Teuchos::null;
}

/*-----------------------------------------------------------------------------*
 |  constructor (public)                                          fbraeu 06/16 |
 |  id          (in)  this element's global id                                 |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::Membrane<distype>::Membrane(int id, int owner)
    : Core::Elements::Element(id, owner),
      thickness_(0.0),
      cur_thickness_(0),
      planetype_(plane_stress),
      intpoints_(Core::FE::GaussRule2D::tri_3point)
{
  switch (distype)
  {
    case Core::FE::CellType::tri3:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::tri_3point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    case Core::FE::CellType::tri6:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::tri_6point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    case Core::FE::CellType::quad4:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::quad_4point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    case Core::FE::CellType::quad9:
    {
      Core::FE::GaussRule2D gaussrule = Core::FE::GaussRule2D::quad_9point;
      // get gauss integration points
      intpoints_ = Core::FE::IntegrationPoints2D(gaussrule);
      break;
    }
    default:
      FOUR_C_THROW("shape type unknown!\n");
      break;
  }
  cur_thickness_.resize(intpoints_.nquad, thickness_);
  return;
}

/*-----------------------------------------------------------------------------*
 |  copy-constructor (public)                                     fbraeu 06/16 |
 |  id               (in)  this element's global id                            |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::Membrane<distype>::Membrane(const Discret::ELEMENTS::Membrane<distype>& old)
    : Core::Elements::Element(old),
      thickness_(old.thickness_),
      cur_thickness_(old.cur_thickness_),
      planetype_(old.planetype_),
      intpoints_(old.intpoints_)
{
  cur_thickness_.resize(intpoints_.nquad, thickness_);
  return;
}

/*------------------------------------------------------------------------*
 |  Deep copy this instance of Membrane and return pointer to it (public) |
 |                                                           fbraeu 06/16 |
 *------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::Membrane<distype>::clone() const
{
  Discret::ELEMENTS::Membrane<distype>* newelement =
      new Discret::ELEMENTS::Membrane<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Core::FE::CellType Discret::ELEMENTS::Membrane<distype>::shape() const
{
  return distype;
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element                     (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::Membrane<distype>::num_line() const
{
  return Core::FE::getNumberOfElementLines(distype);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Membrane<distype>::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Element::pack(data);

  // thickness_
  add_to_pack(data, thickness_);

  // current thickness_
  add_to_pack(data, cur_thickness_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Membrane<distype>::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::unpack(basedata);
  // thickness_
  extract_from_pack(position, data, thickness_);
  // current thickness_
  extract_from_pack(position, data, cur_thickness_);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


/*----------------------------------------------------------------------*
 |  return solid material (public)                         sfuchs 05/17 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<Mat::So3Material> Discret::ELEMENTS::Membrane<distype>::solid_material(
    int nummat) const
{
  return Teuchos::rcp_dynamic_cast<Mat::So3Material>(
      Core::Elements::Element::material(nummat), true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Membrane<distype>::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface");
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Teuchos::RCP<Core::Elements::ParamsInterface>
Discret::ELEMENTS::Membrane<distype>::params_interface_ptr()
{
  return interface_ptr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Solid::ELEMENTS::ParamsInterface& Discret::ELEMENTS::Membrane<distype>::str_params_interface()
{
  if (not is_params_interface()) FOUR_C_THROW("The interface ptr is not set!");
  return *(Teuchos::rcp_dynamic_cast<Solid::ELEMENTS::ParamsInterface>(interface_ptr_, true));
}

/*----------------------------------------------------------------------*
 |  print this element (public)                            fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::Membrane<distype>::print(std::ostream& os) const
{
  os << "Membrane ";
  os << " discretization type: " << Core::FE::CellTypeToString(distype).c_str();
  Element::print(os);
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                           fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Membrane<distype>::lines()
{
  return Core::Communication::ElementBoundaryFactory<MembraneLine<distype>, Membrane<distype>>(
      Core::Communication::buildLines, *this);
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                        fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Membrane<distype>::surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

template class Discret::ELEMENTS::Membrane<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::Membrane<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::Membrane<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::Membrane<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
