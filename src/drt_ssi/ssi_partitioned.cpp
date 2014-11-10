/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned.cpp

 \brief base class for partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "ssi_partitioned.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part::SSI_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams)
  : SSI_Base(comm, globaltimeparams,scatraparams,structparams)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetupSystem()
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetStructSolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetScatraSolution()
{
  structure_->Discretization()->SetState(1,"temperature",scatra_->ScaTraField()->Phinp());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetVelocityFields()
{
  scatra_->ScaTraField()->SetVelocityField(
      zeros_, //convective vel.
      Teuchos::null, //acceleration
      structure_->Velnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      structure_->Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetMeshDisp()
{
  scatra_->ScaTraField()->ApplyMeshMovement(
      structure_->Dispnp(),
      structure_->Discretization());
}

