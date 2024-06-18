/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief one generic (beam-to-?) element pair interacting via potentials

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_beaminteraction_potential_pair.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_beam_potential_pair.hpp"
#include "4C_beaminteraction_beam_to_sphere_potential_pair.hpp"
#include "4C_beaminteraction_potential_params.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_rigidsphere.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamPotentialPair::BeamPotentialPair()
    : isinit_(false),
      issetup_(false),
      beam_potential_params_(Teuchos::null),
      element1_(nullptr),
      element2_(nullptr)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamPotentialParams> params_ptr,
    const Core::Elements::Element* element1, const Core::Elements::Element* element2)
{
  issetup_ = false;

  beam_potential_params_ = params_ptr;

  element1_ = element1;
  element2_ = element2;


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::setup()
{
  check_init();

  // the flag issetup_ will be set in the derived method!
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamPotentialPair> BEAMINTERACTION::BeamPotentialPair::Create(
    std::vector<Core::Elements::Element const*> const& ele_ptrs,
    BEAMINTERACTION::BeamPotentialParams const& beam_potential_params)
{
  // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.
  // numnodalvalues = 1: only positions as primary nodal DoFs ==> Lagrange interpolation
  // numnodalvalues = 2: positions AND tangents ==> Hermite interpolation

  const Discret::ELEMENTS::Beam3Base* beamele1 =
      dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(ele_ptrs[0]);

  // at the moment, both elements of a beam contact pair must be of same type Todo
  const unsigned int numnodes_centerline = beamele1->NumCenterlineNodes();
  const unsigned int numnodalvalues = beamele1->hermite_centerline_interpolation() ? 2 : 1;

  switch (numnodalvalues)
  {
    case 1:
    {
      switch (numnodes_centerline)
      {
        case 2:
        {
          if (ele_ptrs[1]->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
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
          if (ele_ptrs[1]->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
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
          if (ele_ptrs[1]->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
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
          if (ele_ptrs[1]->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
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
          FOUR_C_THROW(
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
          if (ele_ptrs[1]->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
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
          FOUR_C_THROW(
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
      FOUR_C_THROW(
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
void BEAMINTERACTION::BeamPotentialPair::check_init() const
{
  if (not is_init()) FOUR_C_THROW("Call Init() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamPotentialPair::check_init_setup() const
{
  if (not is_init() or not is_setup()) FOUR_C_THROW("Call Init() and setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
Core::FE::GaussRule1D BEAMINTERACTION::BeamPotentialPair::get_gauss_rule() const
{
  switch (Params()->NumberGaussPoints())
  {
    case 5:
    {
      return Core::FE::GaussRule1D::line_5point;
      break;
    }

    case 10:
    {
      return Core::FE::GaussRule1D::line_10point;
      break;
    }

    case 20:
    {
      return Core::FE::GaussRule1D::line_20point;
      break;
    }

    case 32:
    {
      return Core::FE::GaussRule1D::line_32point;
      break;
    }

    case 50:
    {
      return Core::FE::GaussRule1D::line_50point;
      break;
    }

    default:
      FOUR_C_THROW("%d Gauss points are not supported yet!", Params()->NumberGaussPoints());
  }

  return Core::FE::GaussRule1D::undefined;
}

FOUR_C_NAMESPACE_CLOSE
