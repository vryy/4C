/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for partitioned scalar structure interaction

 \level 2

 *------------------------------------------------------------------------------------------------*/
#include "baci_ssi_partitioned.H"

#include "baci_ssi_str_model_evaluator_partitioned.H"

#include "baci_linalg_utils_sparse_algebra_math.H"

#include "baci_adapter_str_structure_new.H"
#include "baci_adapter_scatra_base_algorithm.H"

#include "baci_scatra_timint_implicit.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSIPart::SSIPart(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSIBase(comm, globaltimeparams)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g.
  // redistribution of elements. Only then call the setup to this class. This will call the setup to
  // all classes in the inheritance hierarchy. This way, this class may also override a method that
  // is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string& struct_disname, const std::string& scatra_disname, bool isAle)
{
  // call setup of base class
  SSI::SSIBase::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // safety check
  if (SSIInterfaceMeshtying() and structparams.get<std::string>("PREDICT") != "TangDis")
  {
    dserror(
        "Must have TangDis predictor for structural field in partitioned scalar-structure "
        "interaction simulations involving scatra-scatra interface coupling! Otherwise, Dirichlet "
        "boundary conditions on master-side degrees of freedom are not transferred to slave-side "
        "degrees of freedom!");
  }

  if (IsScaTraManifold()) dserror("Manifold not implemented for partitioned SSI");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIPart::Setup()
{
  // call setup of base class
  SSI::SSIBase::Setup();
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSIPart::SetupModelEvaluator()
{
  // build and register ssi model evaluator
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> ssi_model_ptr =
      Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedSSI(Teuchos::rcp(this, false)));
  StructureBaseAlgorithm()->RegisterModelEvaluator("Partitioned Coupling Model", ssi_model_ptr);

  if (IsS2IKineticsWithPseudoContact()) SetModelevaluatorBaseSSI(ssi_model_ptr);
}
