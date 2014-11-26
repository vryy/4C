/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.cpp

\brief Basis of xfluid-levelset coupling.

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "../drt_scatra/scatra_timint_ost.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_xfluid/xfluid.H"

#include "xfluid_levelset_coupling_algorithm.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false,"scatra",solverparams)
{
  //Note: The ScaTra base algorithm is initialized without a velocity field. This works in this setting, as the level set
  //      field is not propagated. For future use the velocity field should be initialized in this algorithm. However, the
  //      XFluid class will contain enriched nodes. This will have to be dealt with. Inspiration can be found in the combust
  //      algorithm, where this is done.

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::~Algorithm()
{
  return;
}



/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary TPF problem                                    winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("------Stationary-Xfluid-LevelSet-XFEM------  time step ----------------------------------------\n");
  }

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // write Scatra output (fluid output has been already called in FluidField()->Integrate();
  ScaTraField()->Output();

  // run the simulation, calls the xfluid-"integrate()" routine
  FluidField()->Integrate();



//  // solve level set equation
  if (Comm().MyPID()==0)
    std::cout << "/!\\ warning === Level-set field not solved for Fluid_XFEM_LevelSet problems" << std::endl;


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetScaTraValuesInFluid()
{
  // set level set in fluid field

  switch(FluidField()->TimIntScheme())
  {
  case INPAR::FLUID::timeint_stationary:
  case INPAR::FLUID::timeint_one_step_theta:
  {
    // set the new level set value to the fluid
    Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true);

    xfluid->SetLevelSetField(ScaTraField()->Phinp(), ScaTraField()->Discretization());

    // create new state class
    //TODO
    //      xfluid->CreateInitialStateLS();

  }
  break;
  default:
    dserror("Time integration scheme not supported");
    break;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Output()
{
  FluidField()->Output();
  ScaTraField()->Output();

  return;
}


/*----------------------------------------------------------------------*
 | perform result test                                     winter 06/14 |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TestResults()
{
  //Currently no testing for ScaTra values, as nothing is changed from input.
  std::cout << "------------TestResults:--Fluid values are compared.\n";
  std::cout << "--------------------------ScaTra values are not compared as of now.\n";

  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(Comm());


//  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
//  {
//    // DRT::Problem::Instance()->TestAll() is called in level-set field after adding particles
//    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
//  }
//  else
//  {
//    DRT::Problem::Instance()->AddFieldTest(CreateScaTraFieldTest());
//    DRT::Problem::Instance()->TestAll(Comm());
//  }

  return;
}
