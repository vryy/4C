/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned.cpp

 \brief base class for partitioned scalar structure interaction

 \level 2

 \maintainer Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

 *------------------------------------------------------------------------------------------------*/
#include "ssi_partitioned.H"

#include "ssi_str_model_evaluator_partitioned.H"

#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part::SSI_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams)
  : SSI_Base(comm, globaltimeparams)
{
  // Keep this constructor empty!
  // First do everything on the more basic objects like the discretizations, like e.g. redistribution of elements.
  // Only then call the setup to this class. This will call the setup to all classes in the inheritance hierarchy.
  // This way, this class may also override a method that is called during Setup() in a base class.
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetupSystem()
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int SSI::SSI_Part::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname,
    bool isAle)
{
  int returnvar=0;

  // call setup of base class
  returnvar =
      SSI::SSI_Base::Init(comm, globaltimeparams,scatraparams,structparams,struct_disname, scatra_disname,isAle);

  // safety check
  if(scatra_->ScaTraField()->S2ICoupling() and structparams.get<std::string>("PREDICT") != "TangDis")
    dserror("Must have TangDis predictor for structural field in partitioned scalar-structure interaction simulations "
            "involving scatra-scatra interface coupling! Otherwise, Dirichlet boundary conditions on master-side degrees "
            "of freedom are not transferred to slave-side degrees of freedom!");

  return returnvar;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::Setup()
{
  // call setup of base class
  SSI::SSI_Base::Setup();

  return;
}

/*---------------------------------------------------------------------------------*
 | set up structural model evaluator for scalar-structure interaction   fang 01/18 |
 *---------------------------------------------------------------------------------*/
void SSI::SSI_Part::SetupModelEvaluator() const
{
  // build and register ssi model evaluator
  Teuchos::RCP<STR::MODELEVALUATOR::Generic> ssi_model_ptr = Teuchos::rcp(new STR::MODELEVALUATOR::PartitionedSSI(Teuchos::rcp(this,false)));
  struct_adapterbase_ptr_->RegisterModelEvaluator("Partitioned Coupling Model",ssi_model_ptr);

  return;
}
