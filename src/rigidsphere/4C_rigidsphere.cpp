/*----------------------------------------------------------------------------*/
/*! \file

\brief spherical particle element for brownian dynamics

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_rigidsphere.hpp"

#include "4C_beaminteraction_link_pinjointed.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_largerotations.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_browniandyn.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


Discret::ELEMENTS::RigidsphereType Discret::ELEMENTS::RigidsphereType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::RigidsphereType& Discret::ELEMENTS::RigidsphereType::Instance()
{
  return instance_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::Communication::ParObject* Discret::ELEMENTS::RigidsphereType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Rigidsphere* object = new Discret::ELEMENTS::Rigidsphere(-1, -1);
  object->Unpack(data);
  return (object);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::RigidsphereType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RIGIDSPHERE")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Rigidsphere(id, owner));
    return (ele);
  }
  return (Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::RigidsphereType::Create(
    const int id, const int owner)
{
  return (Teuchos::rcp(new Rigidsphere(id, owner)));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RigidsphereType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  nv = 3;
  dimns = 3;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::RigidsphereType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  Core::LinAlg::SerialDenseMatrix nullspace;
  FOUR_C_THROW("method ComputeNullSpace not implemented!");
  return nullspace;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::RigidsphereType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["RIGIDSPHERE"];

  defs["POINT1"] = Input::LineDefinition::Builder()
                       .add_int_vector("POINT1", 1)
                       .add_named_double("RADIUS")
                       .add_named_double("DENSITY")
                       .Build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Rigidsphere::Rigidsphere(int id, int owner)
    : Core::Elements::Element(id, owner), radius_(0.0), rho_(0.0)
{
  mybondstobeams_.clear();
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Rigidsphere::Rigidsphere(const Discret::ELEMENTS::Rigidsphere& old)
    : Core::Elements::Element(old), radius_(old.radius_), rho_(old.rho_)
{
  mybondstobeams_.clear();
  if (old.mybondstobeams_.size())
  {
    for (auto const& iter : old.mybondstobeams_)
    {
      if (iter.second != Teuchos::null)
        mybondstobeams_[iter.first] =
            Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamLinkPinJointed>(iter.second->Clone());
      else
        FOUR_C_THROW("something went wrong, I am sorry. Please go debugging.");
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Rigidsphere and return pointer to it (public) |
 |                                                            meier 05/12 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::Rigidsphere::Clone() const
{
  Discret::ELEMENTS::Rigidsphere* newelement = new Discret::ELEMENTS::Rigidsphere(*this);
  return (newelement);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Rigidsphere::Print(std::ostream& os) const { return; }


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Rigidsphere::Shape() const
{
  return (Core::FE::CellType::point1);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           meier 05/12/
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Rigidsphere::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  Element::Pack(data);

  // add all class variables
  add_to_pack(data, radius_);
  add_to_pack(data, rho_);

  add_to_pack(data, static_cast<int>(mybondstobeams_.size()));
  for (auto const& iter : mybondstobeams_) iter.second->Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Rigidsphere::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  Element::Unpack(basedata);


  // extract all class variables
  extract_from_pack(position, data, radius_);
  extract_from_pack(position, data, rho_);

  int unsigned numbonds = extract_int(position, data);
  for (int unsigned i = 0; i < numbonds; ++i)
  {
    std::vector<char> tmp;
    extract_from_pack(position, data, tmp);
    Teuchos::RCP<Core::Communication::ParObject> object =
        Teuchos::rcp(Core::Communication::Factory(tmp), true);
    Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> link =
        Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamLinkPinJointed>(object);
    if (link == Teuchos::null) FOUR_C_THROW("Received object is not a beam to beam linkage");
    mybondstobeams_[link->Id()] = link;
  }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             meier 02/14|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Rigidsphere::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}


/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 05/12|
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::RigidsphereType::Initialize(Core::FE::Discretization& dis) { return 0; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Rigidsphere::set_params_interface_ptr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<Core::Elements::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Elements::ParamsInterface> Discret::ELEMENTS::Rigidsphere::ParamsInterfacePtr()
{
  return interface_ptr_;
}

FOUR_C_NAMESPACE_CLOSE
