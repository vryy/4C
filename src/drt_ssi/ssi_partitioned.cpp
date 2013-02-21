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

SSI::SSI_Part::SSI_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : SSI_Base(comm, timeparams)
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
  structure_->ApplyTemperatures(scatra_->ScaTraField().Phinp());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetVelocityFields()
{

  scatra_->ScaTraField().SetVelocityField(
      zeros_, //convective vel.
      Teuchos::null, //acceleration
      structure_->ExtractVelnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      structure_->Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part::SetMeshDisp()
{
  scatra_->ScaTraField().ApplyMeshMovement(
      structure_->Dispnp(),
      structure_->Discretization());
}

