/*----------------------------------------------------------------------------*/
/*! \file

\brief one generic (beam-to-?) contact element pair

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "baci_beaminteraction_contact_pair.H"

#include "baci_beam3_base.H"
#include "baci_beaminteraction_beam_to_beam_contact_pair.H"
#include "baci_beaminteraction_beam_to_beam_contact_params.H"
#include "baci_beaminteraction_beam_to_sphere_contact_pair.H"
#include "baci_beaminteraction_conditions.H"
#include "baci_beaminteraction_contact_params.H"
#include "baci_geometry_pair.H"
#include "baci_lib_element.H"
#include "baci_rigidsphere.H"
#include "baci_so3_base.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_RCP.hpp>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactPair::BeamContactPair()
    : isinit_(false),
      issetup_(false),
      geometry_pair_(Teuchos::null),
      params_(Teuchos::null),
      element1_(nullptr),
      element2_(nullptr)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
    std::vector<DRT::Element const*> elements)
{
  issetup_ = false;

  params_ = params_ptr;

  element1_ = elements[0];
  element2_ = elements[1];

  // If a geometry pair is created at this point, call its Init function.
  if (geometry_pair_ != Teuchos::null) geometry_pair_->Init(element1_, element2_);

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::Setup()
{
  CheckInit();

  // If a geometry pair is created by a derived class, call its Setup function.
  if (geometry_pair_ != Teuchos::null) geometry_pair_->Setup();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BEAMINTERACTION::BeamContactPair::Create(
    std::vector<DRT::Element const*> const& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamInteractionConditions>& beam_interaction_conditions_ptr)
{
  // Check the type of the second element.
  const bool other_is_beam = dynamic_cast<DRT::ELEMENTS::Beam3Base const*>(ele_ptrs[1]) != nullptr;
  const bool other_is_solid = dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]) != nullptr ||
                              dynamic_cast<DRT::ELEMENTS::Solid const*>(ele_ptrs[1]) != nullptr;
  const bool other_is_sphere =
      ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance();

  if (other_is_beam or other_is_solid)
  {
    // Beam-to-beam and beam-to-solid pairs are exclusively created by conditions.
    return beam_interaction_conditions_ptr->CreateContactPair(ele_ptrs);
  }
  else if (other_is_sphere)
  {
    // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.
    // numnodalvalues = 1: only positions as primary nodal DoFs ==> Lagrange interpolation
    // numnodalvalues = 2: positions AND tangents ==> Hermite interpolation

    const DRT::ELEMENTS::Beam3Base* beamele1 =
        dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[0]);

    const unsigned int numnodes_centerline = beamele1->NumCenterlineNodes();
    const unsigned int numnodalvalues = beamele1->HermiteCenterlineInterpolation() ? 2 : 1;

    switch (numnodalvalues)
    {
      case 1:
      {
        switch (numnodes_centerline)
        {
          case 2:
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<2, 1>());
          }
          case 3:
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<3, 1>());
          }
          case 4:
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<4, 1>());
          }
          case 5:
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<5, 1>());
          }
          default:
          {
            dserror(
                "%d and %d is no valid template parameter combination for the "
                "number of nodes and number of types of nodal DoFs used for centerline "
                "interpolation!",
                numnodes_centerline, numnodalvalues);
            break;
          }
        }
        break;
      }
      case 2:
      {
        switch (numnodes_centerline)
        {
          case 2:
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<2, 2>());
          }
          default:
          {
            dserror(
                "%d and %d is no valid template parameter combination for the "
                "number of nodes and number of types of nodal DoFs used for centerline "
                "interpolation!",
                numnodes_centerline, numnodalvalues);
            break;
          }
        }
        break;
      }
      default:
      {
        dserror(
            "%d and %d is no valid template parameter combination for the "
            "number of nodes and number of types of nodal DoFs used for centerline "
            "interpolation!",
            numnodes_centerline, numnodalvalues);
        break;
      }
    }
  }
  else
  {
    dserror("Unknown type of second element in creation of beam contact pair.");
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::CheckInit() const
{
  if (not IsInit()) dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}
