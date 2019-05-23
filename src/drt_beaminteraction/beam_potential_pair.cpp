/*-----------------------------------------------------------------------------------------------*/
/*!

\brief one generic (beam-to-?) element pair interacting via potentials

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_potential_pair.H"

#include "../drt_lib/drt_element.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_RCP.hpp>
#include "beam_potential_params.H"
#include "beam_to_beam_potential_pair.H"
#include "beam_to_sphere_potential_pair.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamPotentialPair::BeamPotentialPair()
    : isinit_(false),
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
    const DRT::Element* element1, const DRT::Element* element2)
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
Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> BEAMINTERACTION::BeamPotentialPair::Create(
    std::vector<DRT::Element const*> const& ele_ptrs,
    BEAMINTERACTION::BeamPotentialParams const& beam_potential_params)
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
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSpherePotentialPair<2, 1>());
          else
          {
            if (beam_potential_params.UseFAD())
              return Teuchos::rcp(
                  new BEAMINTERACTION::BeamToBeamPotentialPair<2, 1, Sacado::Fad::DFad<double>>());
            else
              return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialPair<2, 1, double>());
          }
        }
        case 3:
        {
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSpherePotentialPair<3, 1>());
          else
          {
            if (beam_potential_params.UseFAD())
              return Teuchos::rcp(
                  new BEAMINTERACTION::BeamToBeamPotentialPair<3, 1, Sacado::Fad::DFad<double>>());
            else
              return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialPair<3, 1, double>());
          }
        }
        case 4:
        {
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSpherePotentialPair<4, 1>());
          else
          {
            if (beam_potential_params.UseFAD())
              return Teuchos::rcp(
                  new BEAMINTERACTION::BeamToBeamPotentialPair<4, 1, Sacado::Fad::DFad<double>>());
            else
              return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialPair<4, 1, double>());
          }
        }
        case 5:
        {
          if (ele_ptrs[1]->ElementType() == DRT::ELEMENTS::RigidsphereType::Instance())
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSpherePotentialPair<5, 1>());
          else
          {
            if (beam_potential_params.UseFAD())
              return Teuchos::rcp(
                  new BEAMINTERACTION::BeamToBeamPotentialPair<5, 1, Sacado::Fad::DFad<double>>());
            else
              return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialPair<5, 1, double>());
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
            return Teuchos::rcp(new BEAMINTERACTION::BeamToSpherePotentialPair<2, 2>());
          else
          {
            if (beam_potential_params.UseFAD())
              return Teuchos::rcp(
                  new BEAMINTERACTION::BeamToBeamPotentialPair<2, 2, Sacado::Fad::DFad<double>>());
            else
              return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPotentialPair<2, 2, double>());
          }
        }
        default:
          dserror(
              "%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!",
              numnodes_centerline, numnodalvalues);
          break;
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

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::CheckInit() const
{
  if (not IsInit()) dserror("Call Init() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::CheckInitSetup() const
{
  if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
DRT::UTILS::GaussRule1D BEAMINTERACTION::BeamPotentialPair::GetGaussRule() const
{
  switch (Params()->NumberGaussPoints())
  {
    case 5:
    {
      return DRT::UTILS::intrule_line_5point;
      break;
    }

    case 10:
    {
      return DRT::UTILS::intrule_line_10point;
      break;
    }

    case 20:
    {
      return DRT::UTILS::intrule_line_20point;
      break;
    }

    case 32:
    {
      return DRT::UTILS::intrule_line_32point;
      break;
    }

    case 50:
    {
      return DRT::UTILS::intrule_line_50point;
      break;
    }

    default:
      dserror("%d Gauss points are not supported yet!", Params()->NumberGaussPoints());
  }

  return DRT::UTILS::intrule1D_undefined;
}
