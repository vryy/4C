/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for partitioned scalar structure interaction

 \level 2

 *------------------------------------------------------------------------------------------------*/
#include "4C_ssi_partitioned.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_structure_new.hpp"
#include "4C_contact_nitsche_strategy_ssi.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_ssi_str_model_evaluator_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSIPart::SSIPart(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during setup() in a base class.
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart::init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // call setup of base class
  SSI::SSIBase::init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // safety check
  if (ssi_interface_meshtying() and structparams.get<std::string>("PREDICT") != "TangDis")
  {
    FOUR_C_THROW(
        "Must have TangDis predictor for structural field in partitioned scalar-structure "
        "interaction simulations involving scatra-scatra interface coupling! Otherwise, Dirichlet "
        "boundary conditions on master-side degrees of freedom are not transferred to slave-side "
        "degrees of freedom!");
  }

  if (is_sca_tra_manifold()) FOUR_C_THROW("Manifold not implemented for partitioned SSI");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart::setup()
{
  // call setup of base class
  SSI::SSIBase::setup();

  if (SSIInterfaceContact() and !IsRestart())
  {
    setup_contact_strategy();
    ScaTraField()->SetNitscheContactStrategy(nitsche_strategy_ssi());
  }
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIPart::setup_model_evaluator()
{
  // build and register ssi model evaluator
  Teuchos::RCP<Solid::MODELEVALUATOR::Generic> ssi_model_ptr =
      Teuchos::rcp(new Solid::MODELEVALUATOR::PartitionedSSI(Teuchos::rcp(this, false)));
  structure_base_algorithm()->register_model_evaluator("Partitioned Coupling Model", ssi_model_ptr);

  if (is_s2_i_kinetics_with_pseudo_contact()) set_modelevaluator_base_ssi(ssi_model_ptr);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIPart::read_restart(int restart)
{
  // call base class
  SSIBase::read_restart(restart);

  // do ssi contact specific tasks
  if (SSIInterfaceContact())
  {
    setup_contact_strategy();
    ScaTraField()->SetNitscheContactStrategy(nitsche_strategy_ssi());
  }
}

FOUR_C_NAMESPACE_CLOSE
