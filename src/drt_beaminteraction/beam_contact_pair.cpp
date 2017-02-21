/*----------------------------------------------------------------------------*/
/*!
\file beam_contact_pair.cpp

\brief one generic (beam-to-?) contact element pair

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------     */

#include "beam_contact_pair.H"

#include "../drt_lib/drt_element.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>

#include "beam_contact_params.H"
#include "beam_to_beam_contact_pair.H"
#include "beam_to_sphere_contact_pair.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactPair::BeamContactPair() :
    isinit_(false),
    issetup_(false),
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
    const DRT::Element* element1,
    const DRT::Element* element2)
{
  issetup_ = false;

  params_ = params_ptr;

  element1_ = element1;
  element2_ = element2;


  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::Setup()
{
  CheckInit();

  // the flag issetup_ will be set in the derived method!
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamContactPair::Create(
    std::vector< DRT::Element const *> const& ele_ptrs)
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
          if ( ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
            return Teuchos::rcp (new BEAMINTERACTION::BeamToSphereContactPair<2,1>());
          else
            return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamContactPair<2,1>());
        }
        case 3:
        {
          if ( ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
            return Teuchos::rcp (new BEAMINTERACTION::BeamToSphereContactPair<3,1>());
          else
            return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamContactPair<3,1>());
        }
        case 4:
        {
          if ( ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
            return Teuchos::rcp (new BEAMINTERACTION::BeamToSphereContactPair<4,1>());
          else
            return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamContactPair<4,1>());
        }
        case 5:
        {
          if ( ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
            return Teuchos::rcp (new BEAMINTERACTION::BeamToSphereContactPair<5,1>());
          else
            return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamContactPair<5,1>());
        }
        default:
        {
          dserror("%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!", numnodes_centerline, numnodalvalues);
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
          if ( ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance() )
            return Teuchos::rcp (new BEAMINTERACTION::BeamToSphereContactPair<2,2>());
          else
            return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamContactPair<2,2>());
        }
        default:
          dserror("%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!", numnodes_centerline, numnodalvalues);
          break;
      }
      break;
    }
    default:
    {
      dserror("%d and %d is no valid template parameter combination for the "
          "number of nodes and number of types of nodal DoFs used for centerline "
          "interpolation!", numnodes_centerline, numnodalvalues);
      break;
    }
  }

  return Teuchos::null;

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactPair::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup())
    dserror("Call Init() and Setup() first!");
}
