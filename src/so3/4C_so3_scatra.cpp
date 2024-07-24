/*----------------------------------------------------------------------*/
/*! \file

\brief Solid-scatra elements base class

\level 2


*----------------------------------------------------------------------*/


#include "4C_so3_scatra.hpp"

#include "4C_fem_general_element_integration_select.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ssi.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3Scatra<So3Ele, distype>::So3Scatra(int id, int owner)
    : So3Ele(id, owner),
      impltype_(Inpar::ScaTra::impltype_undefined),
      intpoints_(distype == Core::FE::CellType::tet4
                     ? Core::FE::GaussRule3D::tet_1point
                     : (distype == Core::FE::CellType::tet10
                               ? Core::FE::GaussRule3D::tet_4point
                               : Discret::ELEMENTS::DisTypeToOptGaussRule<distype>::rule)),
      numgpt_(intpoints_.nquad),
      xsi_(0.0),
      inv_j_(0.0),
      det_j_(0.0)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3Scatra<So3Ele, distype>::So3Scatra(
    const Discret::ELEMENTS::So3Scatra<So3Ele, distype>& old)
    : So3Ele(old),
      impltype_(old.impltype_),
      intpoints_(old.intpoints_),
      numgpt_(old.numgpt_),
      xsi_(old.xsi_),
      inv_j_(old.inv_j_),
      det_j_(old.det_j_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3Scatra<So3Ele, distype>::clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3Scatra<So3Ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Scatra<So3Ele, distype>::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  So3Ele::add_to_pack(data, type);

  // pack scalar transport impltype
  So3Ele::add_to_pack(data, impltype_);

  // detJ_
  So3Ele::add_to_pack(data, det_j_);

  // invJ_
  auto size = (int)inv_j_.size();
  So3Ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i) So3Ele::add_to_pack(data, inv_j_[i]);

  // xsi_
  size = (int)xsi_.size();
  So3Ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i) So3Ele::add_to_pack(data, xsi_[i]);


  // add base class Element
  So3Ele::pack(data);
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Scatra<So3Ele, distype>::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // extract scalar transport impltype_
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(So3Ele::extract_int(position, data));

  // detJ_
  So3Ele::extract_from_pack(position, data, det_j_);

  // invJ_
  int size = 0;
  So3Ele::extract_from_pack(position, data, size);
  inv_j_.resize(size, Core::LinAlg::Matrix<numdim_, numdim_>(true));
  for (int i = 0; i < size; ++i) So3Ele::extract_from_pack(position, data, inv_j_[i]);

  // xsi_
  size = 0;
  So3Ele::extract_from_pack(position, data, size);
  xsi_.resize(size, Core::LinAlg::Matrix<numdim_, 1>(true));
  for (int i = 0; i < size; ++i) So3Ele::extract_from_pack(position, data, xsi_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  So3Ele::extract_from_pack(position, data, basedata);

  So3Ele::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Scatra<So3Ele, distype>::print(std::ostream& os) const
{
  os << "So3_scatra ";
  os << " discretization type: " << Core::FE::CellTypeToString(distype).c_str();
  Core::Elements::Element::print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  read this element (public)                             schmidt 09/17|
 *----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Scatra<So3Ele, distype>::read_element(const std::string& eletype,
    const std::string& eledistype, const Core::IO::InputParameterContainer& container)
{
  So3Ele::read_element(eletype, eledistype, container);

  // read scalar transport implementation type
  auto impltype = container.get<std::string>("TYPE");

  if (impltype == "Undefined")
    impltype_ = Inpar::ScaTra::impltype_undefined;
  else if (impltype == "AdvReac")
    impltype_ = Inpar::ScaTra::impltype_advreac;
  else if (impltype == "CardMono")
    impltype_ = Inpar::ScaTra::impltype_cardiac_monodomain;
  else if (impltype == "Chemo")
    impltype_ = Inpar::ScaTra::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = Inpar::ScaTra::impltype_chemoreac;
  else if (impltype == "ElchDiffCond")
    impltype_ = Inpar::ScaTra::impltype_elch_diffcond;
  else if (impltype == "ElchElectrode")
    impltype_ = Inpar::ScaTra::impltype_elch_electrode;
  else if (impltype == "Loma")
    impltype_ = Inpar::ScaTra::impltype_loma;
  else if (impltype == "RefConcReac")
    impltype_ = Inpar::ScaTra::impltype_refconcreac;
  else if (impltype == "Std")
    impltype_ = Inpar::ScaTra::impltype_std;
  else
    FOUR_C_THROW("Invalid implementation type for So3_Scatra elements!");

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Core::Nodes::Node** Discret::ELEMENTS::So3Scatra<So3Ele, distype>::nodes()
{
  return So3Ele::nodes();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline Teuchos::RCP<Core::Mat::Material> Discret::ELEMENTS::So3Scatra<So3Ele, distype>::material()
    const
{
  return So3Ele::material();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <class So3Ele, Core::FE::CellType distype>
inline int Discret::ELEMENTS::So3Scatra<So3Ele, distype>::id() const
{
  return So3Ele::id();
}


/*--------------------------------------------------------------------------*
 | set the material  (public)                                 schmidt 10/17 |
 *                                                                          */
template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Scatra<So3Ele, distype>::set_material(
    int matnum, Teuchos::RCP<Core::Mat::Material> mat)
{
  // call base class
  So3Ele::set_material(0, mat);

  // get the scatra structure control parameter list
  const Teuchos::ParameterList& ssicontrol = Global::Problem::instance()->ssi_control_params();

  if ((Teuchos::getIntegralValue<Inpar::SSI::SolutionSchemeOverFields>(ssicontrol, "COUPALGO") ==
          Inpar::SSI::SolutionSchemeOverFields::ssi_Monolithic) and
      (mat->material_type() != Core::Materials::m_multiplicative_split_defgrad_elasthyper))
    FOUR_C_THROW(
        "When you use the 'COUPALGO' 'ssi_Monolithic' from the 'SSI CONTROL' section,"
        " you need to use the material 'MAT_MultiplicativeSplitDefgradElastHyper'! If you want to "
        "use another material feel free to implement it! ;-)");

  // call base class
  So3Ele::set_material(0, mat);

  return;
}


template class Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8, Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex8fbar,
    Core::FE::CellType::hex8>;
template class Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoHex27, Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet4, Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoTet10, Core::FE::CellType::tet10>;
template class Discret::ELEMENTS::So3Scatra<Discret::ELEMENTS::SoWeg6, Core::FE::CellType::wedge6>;

FOUR_C_NAMESPACE_CLOSE
