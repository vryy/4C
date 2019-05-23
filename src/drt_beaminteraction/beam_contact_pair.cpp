/*----------------------------------------------------------------------------*/
/*!

\brief one generic (beam-to-?) contact element pair

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam_contact_pair.H"

#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>

#include "beam_to_beam_contact_pair.H"
#include "beam_to_solid_volume_meshtying_pair_factory.H"
#include "beam_to_sphere_contact_pair.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_so3/so_base.H"
#include "beam_to_beam_contact_params.H"

#include "beam_contact_params.H"

#include "../drt_geometry_pair/geometry_pair.H"
#include "../drt_geometry_pair/geometry_pair_evaluation_data_global.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactPair::BeamContactPair()
    : isinit_(false),
      issetup_(false),
      geometry_pair_(Teuchos::null),
      params_(Teuchos::null),
      element1_(NULL),
      element2_(NULL)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
    std::vector<DRT::Element const*> elements)
{
  issetup_ = false;

  params_ = params_ptr;

  element1_ = elements[0];
  element2_ = elements[1];

  // Create the geometry pair.
  CreateGeometryPair(geometry_evaluation_data_ptr);

  // If a geometry pair is created by a derived class, call its Init function.
  if (geometry_pair_ != Teuchos::null)
    geometry_pair_->Init(geometry_evaluation_data_ptr, element1_, element2_);

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
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr)
{
  // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.
  // numnodalvalues = 1: only positions as primary nodal DoFs ==> Lagrange interpolation
  // numnodalvalues = 2: positions AND tangents ==> Hermite interpolation

  const DRT::ELEMENTS::Beam3Base* beamele1 =
      dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[0]);

  // at the moment, both elements of a beam contact pair must be of same type Todo
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
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<2, 1>());
          }
          else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]) != NULL)
          {
            dserror(
                "ERROR: Case numnodalvalues=1 not yet implemented for "
                "Beam-to-solid volume meshtying)");
          }
          else if (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[1]) != NULL)
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<2, 1>());
          }
          else
          {
            dserror("Unknown type of element for beam-to-? contact evaluation!");
            break;
          }
        }
        case 3:
        {
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<3, 1>());
          }
          else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]) != NULL)
          {
            dserror(
                "ERROR: Case numnodalvalues=1 not yet implemented for "
                "Beam-to-solid volume meshtying)");
          }
          else if (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[1]) != NULL)
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<3, 1>());
          }
          else
          {
            dserror("Unknown type of element for beam-to-? contact evaluation!");
            break;
          }
        }
        case 4:
        {
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<4, 1>());
          }
          else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]) != NULL)
          {
            dserror(
                "ERROR: Case numnodalvalues=1 not yet implemented for "
                "Beam-to-solid volume meshtying)");
          }
          else if (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[1]) != NULL)
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<4, 1>());
          }
          else
          {
            dserror("Unknown type of element for beam-to-? contact evaluation!");
            break;
          }
        }
        case 5:
        {
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<5, 1>());
          }
          else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]) != NULL)
          {
            dserror(
                "ERROR: Case numnodalvalues=1 not yet implemented for "
                "Beam-to-solid volume meshtying)");
          }
          else if (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[1]) != NULL)
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<5, 1>());
          }
          else
          {
            dserror("Unknown type of element for beam-to-? contact evaluation!");
            break;
          }
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
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactPair<2, 2>());
          }
          else if (dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]) != NULL)
          {
            // Create the beam to solid pair.
            return BeamToSolidVolumeMeshtyingPairFactory(ele_ptrs, params_ptr);
          }
          else if (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[1]) != NULL)
          {
            return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<2, 2>());
          }
          else
          {
            dserror("Unknown type of element for beam-to-? contact evaluation!");
            break;
          }
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
