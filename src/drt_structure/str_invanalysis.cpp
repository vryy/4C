/*----------------------------------------------------------------------*/
/*!
\file str_invanalysis.cpp
\brief Inverse analysis for structures

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "str_invanalysis.H"
#include "../drt_contact/contactstrugenalpha.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "stru_resulttest.H"
#include "../drt_inv_analysis/inv_analysis.H"
#include "../drt_inv_analysis/gen_inv_analysis.H"
#include "../drt_inpar/inpar_invanalysis.H"

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
//! \author m.gee \date 06/01
extern GENPROB genprob;

/*======================================================================*/
/* Inverse analysis of structures */
void STR::invanalysis()
{
  // get input lists
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

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

  // make sure we IMR-like generalised-alpha requested for Multi-magic
  if (DRT::Problem::Instance()->ProblemType() == "struct_multi")
  {
    if (Teuchos::getIntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP") != INPAR::STR::dyna_genalpha)
      dserror("For PROBLEMTYP=struct_multi you have to use DYNAMICTYP=GenAlpha");
    else if (Teuchos::getIntegralValue<INPAR::STR::MidAverageEnum>(sdyn.sublist("GENALPHA"), "GENAVG") != INPAR::STR::midavg_imrlike)
      dserror("For PROBLEMTYP=struct_multi you have to use DYNAMICTYP=GenAlpha with GENAVG=ImrLike");
  }

  switch(Teuchos::getIntegralValue<INPAR::STR::InvAnalysisType>(iap,"INV_ANALYSIS"))
  {
    case INPAR::STR::inv_lung:
    {
      STR::InvAnalysis ia(actdis,solver,output);
      ia.Integrate();
    }
    break;
    case INPAR::STR::inv_generalized:
    {
      STR::GenInvAnalysis ia(actdis,solver,output);
      ia.Integrate();
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
#endif  // #ifdef CCADISCRET
