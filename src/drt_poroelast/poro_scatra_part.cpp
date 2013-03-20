/*----------------------------------------------------------------------*/
/*!
 \file POROELAST.cpp

 \brief  scalar transport in porous media

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/

/*
 *  Implementation of passive scalar transport in porous media.
 *  Done by Miguel Urrecha (miguel.urrecha@upm.es)
*/


#include "poro_scatra_part.H"
#include "poro_base.H"

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_inpar/inpar_scatra.H"
#include "poroelast_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part::PORO_SCATRA_Part(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams):
    PORO_SCATRA_Base(comm,timeparams)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part::SetPoroSolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part::SetScatraSolution()
{
  //porous structure
  poro_->StructureField()->ApplyCouplingState(scatra_->ScaTraField().Phinp(),"temperature",2);

  //porous fluid
  poro_->FluidField()->SetIterLomaFields(scatra_->ScaTraField().Phinp(),
                                        scatra_->ScaTraField().Phin(),
                                        scatra_->ScaTraField().Phidtnp(),
                                        Teuchos::null,
                                        0.0,
                                        0.0,
                                        0.0,
                                        0.0,
                                        scatra_->ScaTraField().Discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part::SetVelocityFields()
{
  scatra_->ScaTraField().SetVelocityField(
      poro_->FluidField()->ConvectiveVel(), //convective vel.
      Teuchos::null, //acceleration
      poro_->FluidField()->Velnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      poro_->FluidField()->Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part::SetMeshDisp()
{
  scatra_->ScaTraField().ApplyMeshMovement(
      poro_->FluidField()->Dispnp(),
      poro_->FluidField()->Discretization());
}

