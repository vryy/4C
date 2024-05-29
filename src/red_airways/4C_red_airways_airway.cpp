/*---------------------------------------------------------------------*/
/*! \file

\brief Implements an airway element


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_io_linedefinition.hpp"
#include "4C_io_pstream.hpp"
#include "4C_lib_discret.hpp"
#include "4C_red_airways_elementbase.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

using namespace CORE::FE;

DRT::ELEMENTS::RedAirwayType DRT::ELEMENTS::RedAirwayType::instance_;


DRT::ELEMENTS::RedAirwayType& DRT::ELEMENTS::RedAirwayType::Instance() { return instance_; }


CORE::COMM::ParObject* DRT::ELEMENTS::RedAirwayType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::RedAirway* object = new DRT::ELEMENTS::RedAirway(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::RedAirwayType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RED_AIRWAY")
  {
    Teuchos::RCP<CORE::Elements::Element> ele =
        Teuchos::rcp(new DRT::ELEMENTS::RedAirway(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<CORE::Elements::Element> DRT::ELEMENTS::RedAirwayType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<CORE::Elements::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::RedAirway(id, owner));
  return ele;
}


/*--------------------------------------------------------------------  *
 | Read RED_AIRWAY element line and add element specific parameters     |
 |                                                             (public) |
 |                                                           roth 10/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirwayType::setup_element_definition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["RED_AIRWAY"];

  defs["LINE2"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedString("ElemSolvingType")
                      .AddNamedString("TYPE")
                      .AddNamedString("Resistance")
                      .AddNamedDouble("PowerOfVelocityProfile")
                      .AddNamedDouble("WallElasticity")
                      .AddNamedDouble("PoissonsRatio")
                      .AddNamedDouble("ViscousTs")
                      .AddNamedDouble("ViscousPhaseShift")
                      .AddNamedDouble("WallThickness")
                      .AddNamedDouble("Area")
                      .AddNamedInt("Generation")
                      .add_optional_named_double("AirwayColl")
                      .add_optional_named_double("S_Close")
                      .add_optional_named_double("S_Open")
                      .add_optional_named_double("Pcrit_Open")
                      .add_optional_named_double("Pcrit_Close")
                      .add_optional_named_double("Open_Init")
                      .add_optional_named_double("BranchLength")
                      .Build();
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                           ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirway::RedAirway(int id, int owner) : CORE::Elements::Element(id, owner) {}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      ismail 01/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RedAirway::RedAirway(const DRT::ELEMENTS::RedAirway& old)
    : CORE::Elements::Element(old),
      elem_type_(old.elem_type_),
      resistance_(old.resistance_),
      elemsolving_type_(old.elemsolving_type_),
      airway_params_(old.airway_params_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of RedAirway and return pointer             |
 |  to it                                                      (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
CORE::Elements::Element* DRT::ELEMENTS::RedAirway::Clone() const
{
  DRT::ELEMENTS::RedAirway* newelement = new DRT::ELEMENTS::RedAirway(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::RedAirway::Shape() const
{
  switch (num_node())
  {
    case 2:
      return CORE::FE::CellType::line2;
    case 3:
      return CORE::FE::CellType::line3;
    default:
      FOUR_C_THROW("unexpected number of nodes %d", num_node());
      break;
  }
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirway::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);

  AddtoPack(data, elem_type_);
  AddtoPack(data, resistance_);
  AddtoPack(data, elemsolving_type_);

  AddtoPack(data, airway_params_.power_velocity_profile);
  AddtoPack(data, airway_params_.wall_elasticity);
  AddtoPack(data, airway_params_.poisson_ratio);
  AddtoPack(data, airway_params_.wall_thickness);
  AddtoPack(data, airway_params_.area);
  AddtoPack(data, airway_params_.viscous_Ts);
  AddtoPack(data, airway_params_.viscous_phase_shift);
  AddtoPack(data, airway_params_.branch_length);
  AddtoPack(data, airway_params_.generation);

  AddtoPack(data, airway_params_.airway_coll);
  AddtoPack(data, airway_params_.s_close);
  AddtoPack(data, airway_params_.s_open);
  AddtoPack(data, airway_params_.p_crit_open);
  AddtoPack(data, airway_params_.p_crit_close);
  AddtoPack(data, airway_params_.open_init);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                         ismail 01/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirway::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  ExtractfromPack(position, data, elem_type_);
  ExtractfromPack(position, data, resistance_);
  ExtractfromPack(position, data, elemsolving_type_);

  ExtractfromPack(position, data, airway_params_.power_velocity_profile);
  ExtractfromPack(position, data, airway_params_.wall_elasticity);
  ExtractfromPack(position, data, airway_params_.poisson_ratio);
  ExtractfromPack(position, data, airway_params_.wall_thickness);
  ExtractfromPack(position, data, airway_params_.area);
  ExtractfromPack(position, data, airway_params_.viscous_Ts);
  ExtractfromPack(position, data, airway_params_.viscous_phase_shift);
  ExtractfromPack(position, data, airway_params_.branch_length);
  ExtractfromPack(position, data, airway_params_.generation);

  ExtractfromPack(position, data, airway_params_.airway_coll);
  ExtractfromPack(position, data, airway_params_.s_close);
  ExtractfromPack(position, data, airway_params_.s_open);
  ExtractfromPack(position, data, airway_params_.p_crit_open);
  ExtractfromPack(position, data, airway_params_.p_crit_close);
  ExtractfromPack(position, data, airway_params_.open_init);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                             ismail 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RedAirway::Print(std::ostream& os) const
{
  os << "RedAirway ";
  Element::Print(os);
}


/*----------------------------------------------------------------------*
 |  Return visualization data (public)                     ismail 02/10 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirway::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (CORE::Elements::Element::VisData(name, data)) return true;

  return false;
}


const DRT::REDAIRWAYS::AirwayParams& DRT::ELEMENTS::RedAirway::GetAirwayParams() const
{
  return airway_params_;
}

/*----------------------------------------------------------------------*
 |  get vector of lines              (public)              ismail  02/13|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<CORE::Elements::Element>> DRT::ELEMENTS::RedAirway::Lines()
{
  FOUR_C_ASSERT(NumLine() == 1, "RED_AIRWAY element must have one and only one line");

  return {Teuchos::rcpFromRef(*this)};
}

FOUR_C_NAMESPACE_CLOSE
