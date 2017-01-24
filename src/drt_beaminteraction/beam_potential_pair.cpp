/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_potential_pair.cpp

\brief one generic (beam-to-?) element pair interacting via potentials

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_potential_pair.H"

#include "../drt_lib/drt_element.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>
#include "../drt_beaminteraction/beam_contact_params.H"
#include "beam_to_beam_potential_pair.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamPotentialPair::BeamPotentialPair() :
    isinit_(false),
    issetup_(false),
    beam_potential_params_(Teuchos::null),
    element1_(NULL),
    element2_(NULL)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> params_ptr,
    const DRT::Element* element1,
    const DRT::Element* element2)
{
  issetup_ = false;

  beam_potential_params_ = params_ptr;

  element1_ = element1;
  element2_ = element2;


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::Setup()
{
  CheckInit();

  // the flag issetup_ will be set in the derived method!
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair>
BEAMINTERACTION::BeamPotentialPair::Create(
    unsigned int numnodes,
    unsigned int numnodalvalues)
{
  // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.
  // numnodalvalues = 1: only positions as primary nodal DoFs ==> Lagrange interpolation
  // numnodalvalues = 2: positions AND tangents ==> Hermite interpolation

  switch (numnodalvalues)
  {
    case 1:
    {
      switch (numnodes)
      {
        case 2:
        {
          return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamPotentialPair<2,1>());
        }
        case 3:
        {
          return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamPotentialPair<3,1>());
        }
        case 4:
        {
          return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamPotentialPair<4,1>());
        }
        case 5:
        {
          return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamPotentialPair<5,1>());
        }
        default:
        {
          dserror("%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!", numnodes, numnodalvalues);
          break;
        }
      }
      break;
    }
    case 2:
    {
      switch (numnodes)
      {
        case 2:
        {
          return Teuchos::rcp (new BEAMINTERACTION::BeamToBeamPotentialPair<2,2>());
        }
        default:
          dserror("%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!", numnodes, numnodalvalues);
          break;
      }
      break;
    }
    default:
    {
      dserror("%d and %d is no valid template parameter combination for the "
          "number of nodes and number of types of nodal DoFs used for centerline "
          "interpolation!", numnodes, numnodalvalues);
      break;
    }
  }

  return Teuchos::null;

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup())
    dserror("Call Init() and Setup() first!");
}
