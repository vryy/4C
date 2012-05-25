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
    const int disnum                   ///< optimization field discretization number (default: 0)
)
{
  // setup topology optimization algorithm

  // -------------------------------------------------------------------
  // access the fluid and the optimization discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> optidis = null;
  optidis = DRT::Problem::Instance()->Dis(genprob.numof,disnum);
  RCP<DRT::Discretization> fluiddis = null;
  fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);

  // -------------------------------------------------------------------
  // check degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!optidis->Filled()) dserror("optimization discretization should be filled before");
  if (!fluiddis->Filled()) dserror("fluid discretization should be filled before");

//  // -------------------------------------------------------------------
//  // context for output and restart
//  // -------------------------------------------------------------------
//  RCP<IO::DiscretizationWriter> output =
//    rcp(new IO::DiscretizationWriter(actdis));
//  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // create instance of the optimization class (call the constructor)
  // -------------------------------------------------------------------
  optimizer_ = rcp(new TOPOPT::Optimizer(optidis,fluiddis,prbdyn));

  return;

}
