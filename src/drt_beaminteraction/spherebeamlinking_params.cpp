/*----------------------------------------------------------------------------*/
/*!
\file spherebeamlinking_params.cpp

\brief data container holding all contractile cells input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------------*/

#include "spherebeamlinking_params.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../drt_inpar/inpar_beaminteraction.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/crosslinkermat.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::SphereBeamLinkingParams::SphereBeamLinkingParams()
  : isinit_(false),
    issetup_(false),
    maxnumintegrins_(-1),
    deltatime_(-1.0),
    contractionrate_(-1.0),
    mat_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SphereBeamLinkingParams::Init( STR::TIMINT::BaseDataGlobalState const & gstate )
{
  issetup_ = false;

  const Teuchos::ParameterList& spherebeamlink_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("SPHERE BEAM LINK");

  // number of cells
  maxnumintegrins_ = spherebeamlink_params_list.get<int> ("MAXNUMINTEGRINSPERCELL");

  // time step for stochastic events concering crosslinking
  deltatime_ = spherebeamlink_params_list.get<double> ("TIMESTEP");

  // contraction rate per second in % [ 0, 1 ] of initial length per second
  contractionrate_ = spherebeamlink_params_list.get<double> ("CONTRACTIONRATE");

  // safety check
  // todo: maybe make input of time step obligatory
  if ( deltatime_ < 0.0 )
  {
    deltatime_ = (*gstate.GetDeltaTime())[0];
    if ( gstate.GetMyRank() == 0 )
      std::cout << " Time step " << (*gstate.GetDeltaTime())[0] << " from Structural Dynamic section "
          "used for sphere beam link.\n"
          "Force dependent unbinding of beam-sphere linker is activated for dt > 0" << std::endl;
  }

  // get material for sphere beam linking
  mat_ = Teuchos::rcp_dynamic_cast<MAT::CrosslinkerMat>( MAT::Material::Factory(
      spherebeamlink_params_list.get<int> ("MATINTEGRINS") ) );
  if ( mat_ == Teuchos::null )
    dserror("Invalid material given for beam sphere link. \n");

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SphereBeamLinkingParams::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}
