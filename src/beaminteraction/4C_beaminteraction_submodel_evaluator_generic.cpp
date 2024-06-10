/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all beaminteraction submodel evaluators.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_beaminteraction_submodel_evaluator_generic.hpp"

#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_crosslinker_handler.hpp"
#include "4C_beaminteraction_periodic_boundingbox.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      discret_ptr_(Teuchos::null),
      bindis_ptr_(Teuchos::null),
      gstate_ptr_(Teuchos::null),
      gio_ptr_(Teuchos::null),
      beaminteractiondatastate_(Teuchos::null),
      beam_crosslinker_handler_(Teuchos::null),
      binstrategy_(Teuchos::null),
      periodic_boundingbox_(Teuchos::null),
      eletypeextractor_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Init(
    Teuchos::RCP<Core::FE::Discretization> const& ia_discret,
    Teuchos::RCP<Core::FE::Discretization> const& bindis,
    Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> const& gstate,
    Teuchos::RCP<STR::TimeInt::BaseDataIO> const& gio_ptr,
    Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState> const& ia_gstate_ptr,
    Teuchos::RCP<BEAMINTERACTION::BeamCrosslinkerHandler> const& beamcrosslinkerhandler,
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy,
    Teuchos::RCP<Core::Geo::MeshFree::BoundingBox> const& periodic_boundingbox,
    Teuchos::RCP<BEAMINTERACTION::UTILS::MapExtractor> const& eletypeextractor)
{
  issetup_ = false;

  discret_ptr_ = ia_discret;
  bindis_ptr_ = bindis;
  gstate_ptr_ = gstate;
  gio_ptr_ = gio_ptr;
  beaminteractiondatastate_ = ia_gstate_ptr;
  beam_crosslinker_handler_ = beamcrosslinkerhandler;
  binstrategy_ = binstrategy;
  periodic_boundingbox_ = periodic_boundingbox;
  eletypeextractor_ = eletypeextractor;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Generic::check_init_setup() const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::SUBMODELEVALUATOR::Generic::check_init() const
{
  if (not is_init()) FOUR_C_THROW("Call Init() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::Discretization& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Discret()
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::DiscretPtr()
{
  check_init();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::FE::Discretization>
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::DiscretPtr() const
{
  check_init();
  return discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::Discretization const& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::Discret() const
{
  check_init();
  return *discret_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::Discretization& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscret()
{
  check_init();
  return *bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization>& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscretPtr()
{
  check_init();
  return bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::FE::Discretization>
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscretPtr() const
{
  check_init();
  return bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::FE::Discretization const& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinDiscret() const
{
  check_init();
  return *bindis_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::BaseDataGlobalState& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GState()
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GStatePtr()
{
  check_init();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::BaseDataGlobalState const& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GState() const
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::BaseDataIO& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GInOutput()
{
  check_init();
  return *gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::BaseDataIO const& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::GInOutput() const
{
  check_init();
  return *gio_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteractionDataState&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::beam_interaction_data_state()
{
  check_init();
  return *beaminteractiondatastate_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::MODELEVALUATOR::BeamInteractionDataState>&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::beam_interaction_data_state_ptr()
{
  check_init();
  return beaminteractiondatastate_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::BeamInteractionDataState const&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::beam_interaction_data_state() const
{
  check_init();
  return *beaminteractiondatastate_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamCrosslinkerHandler&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::beam_crosslinker_handler()
{
  check_init();
  return *beam_crosslinker_handler_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::BeamCrosslinkerHandler>&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::beam_crosslinker_handler_ptr()
{
  check_init();
  return beam_crosslinker_handler_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy const& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinStrategy() const
{
  check_init();
  return *binstrategy_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BINSTRATEGY::BinningStrategy& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinStrategy()
{
  check_init();
  return *binstrategy_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BINSTRATEGY::BinningStrategy>&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::BinStrategyPtr()
{
  check_init();
  return binstrategy_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamCrosslinkerHandler const&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::beam_crosslinker_handler() const
{
  check_init();
  return *beam_crosslinker_handler_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::MeshFree::BoundingBox& BEAMINTERACTION::SUBMODELEVALUATOR::Generic::PeriodicBoundingBox()
{
  check_init();
  return *periodic_boundingbox_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Geo::MeshFree::BoundingBox>&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::periodic_bounding_box_ptr()
{
  check_init();
  return periodic_boundingbox_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Core::Geo::MeshFree::BoundingBox const&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::PeriodicBoundingBox() const
{
  check_init();
  return *periodic_boundingbox_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::UTILS::MapExtractor&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::EleTypeMapExtractor()
{
  check_init();
  eletypeextractor_->check_for_valid_map_extractor();
  return *eletypeextractor_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<BEAMINTERACTION::UTILS::MapExtractor>&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::ele_type_map_extractor_ptr()
{
  check_init();
  eletypeextractor_->check_for_valid_map_extractor();
  return eletypeextractor_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::UTILS::MapExtractor const&
BEAMINTERACTION::SUBMODELEVALUATOR::Generic::EleTypeMapExtractor() const
{
  check_init();
  eletypeextractor_->check_for_valid_map_extractor();
  return *eletypeextractor_;
}

FOUR_C_NAMESPACE_CLOSE
