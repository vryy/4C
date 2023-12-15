/*----------------------------------------------------------------------------*/
/*! \file

\brief spherical particle element for brownian dynamics

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "baci_rigidsphere.H"

#include "baci_beaminteraction_link_pinjointed.H"
#include "baci_comm_utils_factory.H"
#include "baci_discretization_fem_general_largerotations.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_inpar_browniandyn.H"
#include "baci_inpar_validparameters.H"
#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_structure_new_elements_paramsinterface.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN


DRT::ELEMENTS::RigidsphereType DRT::ELEMENTS::RigidsphereType::instance_;

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::RigidsphereType& DRT::ELEMENTS::RigidsphereType::Instance() { return instance_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::COMM::ParObject* DRT::ELEMENTS::RigidsphereType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Rigidsphere* object = new DRT::ELEMENTS::Rigidsphere(-1, -1);
  object->Unpack(data);
  return (object);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RigidsphereType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "RIGIDSPHERE")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Rigidsphere(id, owner));
    return (ele);
  }
  return (Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Element> DRT::ELEMENTS::RigidsphereType::Create(const int id, const int owner)
{
  return (Teuchos::rcp(new Rigidsphere(id, owner)));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RigidsphereType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = 3;
  nv = 3;
  dimns = 3;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::RigidsphereType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  CORE::LINALG::SerialDenseMatrix nullspace;
  dserror("method ComputeNullSpace not implemented!");
  return nullspace;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::RigidsphereType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, DRT::INPUT::LineDefinition>& defs = definitions["RIGIDSPHERE"];

  defs["POINT1"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("POINT1", 1)
                       .AddNamedDouble("RADIUS")
                       .AddNamedDouble("DENSITY")
                       .Build();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Rigidsphere::Rigidsphere(int id, int owner)
    : DRT::Element(id, owner), radius_(0.0), rho_(0.0)
{
  mybondstobeams_.clear();
}
/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       meier 05/12|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Rigidsphere::Rigidsphere(const DRT::ELEMENTS::Rigidsphere& old)
    : DRT::Element(old), radius_(old.radius_), rho_(old.rho_)
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
        dserror("something went wrong, I am sorry. Please go debugging.");
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Rigidsphere and return pointer to it (public) |
 |                                                            meier 05/12 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Rigidsphere::Clone() const
{
  DRT::ELEMENTS::Rigidsphere* newelement = new DRT::ELEMENTS::Rigidsphere(*this);
  return (newelement);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                              meier 05/12
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::Print(std::ostream& os) const { return; }


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          meier 05/12 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Rigidsphere::Shape() const
{
  return (CORE::FE::CellType::point1);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           meier 05/12/
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  Element::Pack(data);

  // add all class variables
  AddtoPack(data, radius_);
  AddtoPack(data, rho_);

  AddtoPack(data, static_cast<int>(mybondstobeams_.size()));
  for (auto const& iter : mybondstobeams_) iter.second->Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           meier 05/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);


  // extract all class variables
  ExtractfromPack(position, data, radius_);
  ExtractfromPack(position, data, rho_);

  int unsigned numbonds = ExtractInt(position, data);
  for (int unsigned i = 0; i < numbonds; ++i)
  {
    std::vector<char> tmp;
    ExtractfromPack(position, data, tmp);
    Teuchos::RCP<CORE::COMM::ParObject> object = Teuchos::rcp(CORE::COMM::Factory(tmp), true);
    Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> link =
        Teuchos::rcp_dynamic_cast<BEAMINTERACTION::BeamLinkPinJointed>(object);
    if (link == Teuchos::null) dserror("Received object is not a beam to beam linkage");
    mybondstobeams_[link->Id()] = link;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             meier 02/14|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Rigidsphere::Lines()
{
  return {Teuchos::rcpFromRef(*this)};
}


/*----------------------------------------------------------------------*
 |  Initialize (public)                                      meier 05/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::RigidsphereType::Initialize(DRT::Discretization& dis) { return 0; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Rigidsphere::SetParamsInterfacePtr(const Teuchos::ParameterList& p)
{
  if (p.isParameter("interface"))
    interface_ptr_ = Teuchos::rcp_dynamic_cast<STR::ELEMENTS::ParamsInterface>(
        p.get<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface"));
  else
    interface_ptr_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> DRT::ELEMENTS::Rigidsphere::ParamsInterfacePtr()
{
  return interface_ptr_;
}

BACI_NAMESPACE_CLOSE
