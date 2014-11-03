/*----------------------------------------------------------------------*/
/*!
\file ad_fld_moving_boundary.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "ad_fld_fluid_ale.H"
#include "ad_fld_fluid_xfem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidMovingBoundaryBaseAlgorithm::FluidMovingBoundaryBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn,
    std::string condname
    )
{
    const PROBLEM_TYP probtyp = DRT::Problem::Instance()->ProblemType();

    // switch between moving domain fluid implementations
    switch (probtyp)
    {
    case prb_fsi:
    case prb_fluid_fluid_ale:
    case prb_fluid_ale:
    case prb_freesurf:
    case prb_fsi_redmodels:
    {
      //std::cout << "using FluidAle as FluidMovingBoundary" << endl;
      fluid_ = Teuchos::rcp(new FluidAle(prbdyn,condname));
      break;
    }
    case prb_fluid_xfem:
    case prb_fsi_xfem:
    case prb_fpsi_xfem:
    case prb_fsi_crack:
    case prb_immersed_fsi:
    {
      //std::cout << "using FluidXFEM as FluidMovingBoundary" << endl;
      fluid_ = Teuchos::rcp(new FluidXFEM(prbdyn,condname));
      break;
    }
    default:
      dserror("fsi type not supported"); break;
    }
}


