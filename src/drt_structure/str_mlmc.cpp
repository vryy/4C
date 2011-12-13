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
#ifdef CCADISCRET
#ifdef HAVE_FFTW

/*----------------------------------------------------------------------*/
/* headers */
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif


#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "stru_resulttest.H"
#include "str_mlmc.H"
#include "../drt_mlmc/mlmc.H"
#include "../drt_inpar/inpar_mlmc.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../drt_adapter/adapter_structure_timint.H"

#include "strtimint.H"
#include "strtimint_impl.H"
#include "strtimint_expl.H"
#include "strtimint_genalpha.H"
#include "strtimint_ost.H"
#include "strtimint_gemm.H"
#include "strtimint_ab2.H"

#include "strtimada.H"
#include "strtimada_zienxie.H"
#include "strtimada_joint.H"

#include "strtimint_create.H"
#include "strtimada_create.H"


/*----------------------------------------------------------------------*/
//! General problem data
//!
//! global variable GENPROB genprob is defined in global_control.c
extern GENPROB genprob;

/*======================================================================*/
/* Multi Level Monte Carlo analysis of structures */
void STR::mlmc()
{
  // get input lists
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf, 0);

  // set degrees of freedom in the discretization
  if (not actdis->Filled()) actdis->FillComplete();

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output
    = Teuchos::rcp(new IO::DiscretizationWriter(actdis));

  // input parameters for structural dynamics
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();

  // show default parameters
  if (actdis->Comm().MyPID() == 0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // create a solver
  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->StructSolverParams(),
                                      actdis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());
  bool perform_mlmc = Teuchos::getIntegralValue<int>(mlmcp,"MLMC");
  //cout << perform_mlmc << "see if we come in here" << __LINE__ << __FILE__ <<endl;
  if (perform_mlmc==true)
  {
    STR::MLMC mc(actdis,solver,output);
    mc.Integrate();
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
#endif  // #ifdef CCADISCRET
