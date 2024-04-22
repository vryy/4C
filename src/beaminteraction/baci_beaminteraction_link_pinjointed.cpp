/*----------------------------------------------------------------------*/
/*! \file

\brief connecting beam linked by pin joint

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_beaminteraction_link_pinjointed.hpp"

#include "baci_beaminteraction_link.hpp"
#include "baci_beaminteraction_link_beam3_reissner_line2_pinjointed.hpp"
#include "baci_beaminteraction_link_truss.hpp"
#include "baci_discretization_fem_general_largerotations.hpp"
#include "baci_inpar_beaminteraction.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


BEAMINTERACTION::BeamLinkPinJointedType BEAMINTERACTION::BeamLinkPinJointedType::instance_;


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkPinJointed::BeamLinkPinJointed() : BeamLink() {}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BEAMINTERACTION::BeamLinkPinJointed::BeamLinkPinJointed(
    const BEAMINTERACTION::BeamLinkPinJointed& old)
    : BEAMINTERACTION::BeamLink(old)
{
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Init(int id,
    const std::vector<std::pair<int, int>>& eleids,
    const std::vector<CORE::LINALG::Matrix<3, 1>>& initpos,
    const std::vector<CORE::LINALG::Matrix<3, 3>>& inittriad,
    INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset)
{
  issetup_ = false;

  BeamLink::Init(id, eleids, initpos, inittriad, linkertype, timelinkwasset);

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Setup(const int matnum)
{
  CheckInit();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // add base class Element
  BeamLink::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  BeamLink::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::ResetState(
    std::vector<CORE::LINALG::Matrix<3, 1>>& bspotpos,
    std::vector<CORE::LINALG::Matrix<3, 3>>& bspottriad)
{
  CheckInitSetup();

  BeamLink::ResetState(bspotpos, bspottriad);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> BEAMINTERACTION::BeamLinkPinJointed::Create(
    INPAR::BEAMINTERACTION::JointType type)
{
  if (type == INPAR::BEAMINTERACTION::beam3r_line2_pin)
    return Teuchos::rcp(new BEAMINTERACTION::BeamLinkBeam3rLine2PinJointed());
  else if (type == INPAR::BEAMINTERACTION::truss)
    return Teuchos::rcp(new BEAMINTERACTION::BeamLinkTruss());
  else
    FOUR_C_THROW(
        "instantiation of new BeamLinkPinJointed object failed due to "
        "unknown type of linker");

  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Print(std::ostream& out) const
{
  CheckInit();

  BeamLink::Print(out);

  out << "\n";
}

FOUR_C_NAMESPACE_CLOSE
