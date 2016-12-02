/*----------------------------------------------------------------------*/
/*!
\file str_invanalysis.cpp
\brief Inverse analysis for structures

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "str_invanalysis.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "stru_resulttest.H"
#include "../drt_inv_analysis/inv_analysis.H"
#include "../drt_inv_analysis/gen_inv_analysis.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_comm/comm_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

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



/*======================================================================*/
/* Inverse analysis of structures */
void STR::invanalysis()
{
  // get input lists
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  // access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("structure");

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();
  if (!actdis->HaveDofs()) actdis->FillComplete();

  // context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();

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

  Teuchos::RCP<LINALG::Solver> solver
    = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                      actdis->Comm(),
                                      DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  switch(DRT::INPUT::IntegralValue<INPAR::STR::InvAnalysisType>(iap,"INV_ANALYSIS"))
  {
    case INPAR::STR::inv_lung:
    {
      STR::InvAnalysis ia(actdis,solver,output);
      ia.Integrate();
    }
    break;
    case INPAR::STR::inv_generalized:
    {
      int ngroup = DRT::Problem::Instance()->GetNPGroup()->NumGroups();
      // check whether there is a micro scale which is equivalent to have a subcomm
      Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
      if(subcomm!=Teuchos::null and ngroup>2)
        dserror("Nested parallelism with more than two groups not yet available for inverse multiscale problems");
      STR::GenInvAnalysis ia(actdis,solver,output);
      if (ngroup==1 or (subcomm!=Teuchos::null and ngroup==2)) ia.Integrate();
      else                                                     ia.NPIntegrate();
    }
    break;
    default:
      dserror("Unknown type of inverse analysis");
    break;
  }

  // done
  return;
} // end str_invanalysis()


/*----------------------------------------------------------------------*/
