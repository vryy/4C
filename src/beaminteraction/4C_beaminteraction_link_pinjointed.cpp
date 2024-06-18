/*----------------------------------------------------------------------*/
/*! \file

\brief connecting beam linked by pin joint

\level 3

*/
/*----------------------------------------------------------------------*/

#include "4C_beaminteraction_link_pinjointed.hpp"

#include "4C_beaminteraction_link.hpp"
#include "4C_beaminteraction_link_beam3_reissner_line2_pinjointed.hpp"
#include "4C_beaminteraction_link_truss.hpp"
#include "4C_fem_general_largerotations.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_utils_exceptions.hpp"

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
    const std::vector<Core::LinAlg::Matrix<3, 1>>& initpos,
    const std::vector<Core::LinAlg::Matrix<3, 3>>& inittriad,
    Inpar::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset)
{
  issetup_ = false;

  BeamLink::Init(id, eleids, initpos, inittriad, linkertype, timelinkwasset);

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::setup(const int matnum)
{
  check_init();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // add base class Element
  BeamLink::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  BeamLink::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamLinkPinJointed::ResetState(
    std::vector<Core::LinAlg::Matrix<3, 1>>& bspotpos,
    std::vector<Core::LinAlg::Matrix<3, 3>>& bspottriad)
{
  check_init_setup();

  BeamLink::ResetState(bspotpos, bspottriad);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamLinkPinJointed> BEAMINTERACTION::BeamLinkPinJointed::Create(
    Inpar::BEAMINTERACTION::JointType type)
{
  if (type == Inpar::BEAMINTERACTION::beam3r_line2_pin)
    return Teuchos::rcp(new BEAMINTERACTION::BeamLinkBeam3rLine2PinJointed());
  else if (type == Inpar::BEAMINTERACTION::truss)
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
  check_init();

  BeamLink::Print(out);

  out << "\n";
}

FOUR_C_NAMESPACE_CLOSE
