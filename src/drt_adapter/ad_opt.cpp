/*!------------------------------------------------------------------------------------------------*
\file ad_opt.cpp

\brief topology optimization adapter

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "ad_opt.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_opti/topopt_optimizer.H"


/// constructor
ADAPTER::TopOptBaseAlgorithm::TopOptBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, ///< problem-dependent parameters
    const int disnum                   ///< fluid field discretization number (default: 0)
)
{
  // setup topology optimization algorithm

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numff,disnum);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) dserror("fluid discretization should be filled before");

//  // -------------------------------------------------------------------
//  // context for output and restart
//  // -------------------------------------------------------------------
//  RCP<IO::DiscretizationWriter> output =
//    rcp(new IO::DiscretizationWriter(actdis));
//  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // create instance of the optimization class (call the constructor)
  // -------------------------------------------------------------------
  optimizer_ = rcp(new TOPOPT::Optimizer(actdis,prbdyn));

  return;

}
