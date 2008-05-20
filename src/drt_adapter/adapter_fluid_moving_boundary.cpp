/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_moving_boundary.cpp

\brief 

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "adapter_fluid_ale.H"
#include "adapter_fluid_xfem.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidMovingBoundaryBaseAlgorithm::FluidMovingBoundaryBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
                                                              std::string condname)
{
    const Teuchos::ParameterList& list = DRT::Problem::Instance()->ProblemTypeParams();
    const PROBLEM_TYP probtyp = Teuchos::getIntegralValue<PROBLEM_TYP>(list,"PROBLEMTYP");

    // switch between moving domain fluid implementations
    switch (probtyp)
    {
    case prb_fsi:
    case prb_fluid_ale:
    case prb_freesurf:
    {
      std::cout << "using FluidAle as FluidMovingBoundary" << endl;
      fluid_ = Teuchos::rcp(new FluidAle(prbdyn,condname));
      break;
    }
    case prb_fsi_xfem:
    {
      std::cout << "using FluidXFEM as FluidMovingBoundary" << endl;
      fluid_ = Teuchos::rcp(new FluidXFEM(prbdyn,condname));
      break;
    }
    default:
      dserror("fsi type not supported");
      exit(1);
    }
}


#endif
