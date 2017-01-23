/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_generic.cpp

\brief Generic class for all beaminteraction submodel evaluators.

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#include "beaminteraction_submodel_evaluator_generic.H"

#include "../drt_beaminteraction/periodic_boundingbox.H"
#include "../drt_lib/drt_dserror.H"
#include "str_model_evaluator_beaminteraction_datastate.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      discret_ptr_(Teuchos::null),
      bindis_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      beaminteractiondatastate_(Teuchos::null),
      particlehandler_(Teuchos::null),
      periodic_boundingbox_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Init(
    Teuchos::RCP<DRT::Discretization> const& ia_discret,
    Teuchos::RCP<DRT::Discretization> const& bindis,
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> const& gstate,
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> const& ia_gstate_ptr,
    Teuchos::RCP<PARTICLE::ParticleHandler> const& particlehandler,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& periodic_boundingbox)
{
  issetup_ = false;

  discret_ptr_ = ia_discret;
  bindis_ptr_ = bindis;
  gstate_ptr_ = gstate;
  beaminteractiondatastate_ = ia_gstate_ptr;
  particlehandler_ = particlehandler;
  periodic_boundingbox_ = periodic_boundingbox;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Generic::CheckInitSetup() const
{
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Generic::CheckInit() const
{
  if (not IsInit())
    dserror("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Discretization& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Discret()
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::DiscretPtr()
{
  CheckInit();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Discretization& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Discret() const
{
  CheckInit();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Discretization& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscret()
{
  CheckInit();
  return *bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscretPtr()
{
  CheckInit();
  return bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const DRT::Discretization& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscret() const
{
  CheckInit();
  return *bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataGlobalState& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GState()
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GStatePtr()
{
  CheckInit();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TIMINT::BaseDataGlobalState& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GState()
    const
{
  CheckInit();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteractionDataState& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BeamInteractionDataState()
{
  CheckInit();
  return *beaminteractiondatastate_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BeamInteractionDataStatePtr()
{
  CheckInit();
  return beaminteractiondatastate_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::BeamInteractionDataState& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BeamInteractionDataState()
    const
{
  CheckInit();
  return *beaminteractiondatastate_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
PARTICLE::ParticleHandler& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::ParticleHandler()
{
  CheckInit();
  return *particlehandler_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<PARTICLE::ParticleHandler>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::ParticleHandlerPtr()
{
  CheckInit();
  return particlehandler_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const PARTICLE::ParticleHandler& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::ParticleHandler() const
{
  CheckInit();
  return *particlehandler_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::MESHFREE::BoundingBox& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::PeriodicBoundingBox()
{
  CheckInit();
  return *periodic_boundingbox_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<GEO::MESHFREE::BoundingBox>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::PeriodicBoundingBoxPtr()
{
  CheckInit();
  return periodic_boundingbox_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const GEO::MESHFREE::BoundingBox& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::PeriodicBoundingBox() const
{
  CheckInit();
  return *periodic_boundingbox_;
}
