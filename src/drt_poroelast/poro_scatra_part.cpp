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
void POROELAST::PORO_SCATRA_Part::DoPoroStep()
{
  //1)  solve the step problem. Methods obtained from poroelast->TimeLoop(sdynparams); --> sdynparams
  //			CUIDADO, aqui vuelve a avanzar el paso de tiempo. Hay que corregir eso.
  //2)  Newton-Raphson iteration
  //3)  calculate stresses, strains, energies
  //4)  update all single field solvers
  //5)  write output to screen and files
  poro_-> Solve();
}

void POROELAST::PORO_SCATRA_Part::DoScatraStep()
{
  const Epetra_Comm& comm =
      DRT::Problem::Instance()->GetDis("structure")->Comm();

  if (comm.MyPID() == 0)
  {
    cout
        << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }
  // -------------------------------------------------------------------
  // prepare time step
  // -------------------------------------------------------------------
  scatra_->ScaTraField().PrepareTimeStep();

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField().Solve();

  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  scatra_->ScaTraField().Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField().Output();
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
void POROELAST::PORO_SCATRA_Part::SetVelocityFields()
{
  scatra_->ScaTraField().SetVelocityField(
      poro_->FluidField().ConvectiveVel(), //convective vel.
      Teuchos::null, //acceleration
      poro_->FluidField().Velnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      poro_->FluidField().Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part::SetMeshDisp()
{
  scatra_->ScaTraField().ApplyMeshMovement(
      poro_->FluidField().Dispnp(),
      poro_->FluidField().Discretization());
}

