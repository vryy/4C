/*----------------------------------------------------------------------*/
/*!
\file str_mlmc.cpp
\brief Multilevel Monte Carlo analysis for structures

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef HAVE_FFTW

/*----------------------------------------------------------------------*/

/* headers */

#include "str_mlmc.H"
#include "../drt_mlmc/mlmc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"


/*======================================================================*/
/* Multi Level Monte Carlo analysis of structures */
void STR::mlmc()
{
  // get input lists
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // set degrees of freedom in the discretization
  if (not actdis->Filled() || not actdis->HaveDofs()) actdis->FillComplete();



  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));

  // input parameters for structural dynamics
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();

  // show default parameters
  if (actdis->Comm().MyPID() == 0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, sdyn);

  // create a solver
  // get the solver number used for structural solver
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  //Teuchos::RCP<LINALG::Solver> solver
  //  = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
   //                                   actdis->Comm(),
   //                                   NULL));
   //                                   //DRT::Problem::Instance()->ErrorFile()->Handle()));
 // actdis->ComputeNullSpaceIfNecessary(solver->Params());
  bool perform_mlmc = Teuchos::getIntegralValue<int>(mlmcp,"MLMC");
  bool param_cont =Teuchos::getIntegralValue<int>(mlmcp,"PARAMETERCONTINUATION");
  if (perform_mlmc==true)
  {
    STR::MLMC mc(actdis,output);
    // Use another integrate function if we do not reset the prestress
    if(!param_cont)
    {
    	mc.Integrate();
    }
    else
    {
    	mc.IntegrateNoReset();
    }
  }
  else
  {
    dserror("Unknown type of Multi Level Monte Carlo Analysis");
  }

  // done
  return;
} // end str_mlmc()


/*----------------------------------------------------------------------*/
#endif
