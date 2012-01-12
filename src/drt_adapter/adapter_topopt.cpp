/*!------------------------------------------------------------------------------------------------*
\file adapter_fluid_topopt.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_opti/topopt_optimizer.H"
#include "adapter_topopt.H"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/// constructor
ADAPTER::TopOptBaseAlgorithm::TopOptBaseAlgorithm(
    const Teuchos::ParameterList& prbdyn, ///< problem-dependent parameters
    const int disnum                   ///< scatra discretization number (default: 0)
)
{
  // setup topology optimization algorithm

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numof,disnum);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

//  // -------------------------------------------------------------------
//  // context for output and restart
//  // -------------------------------------------------------------------
//  RCP<IO::DiscretizationWriter> output =
//    rcp(new IO::DiscretizationWriter(actdis));
//  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------

  const Teuchos::ParameterList& optimizationParams =
    DRT::Problem::Instance()->OptimizationControlParams();

  // create instance of the optimization class (call the constructor)
  optimizer_ = rcp(new TOPOPT::Optimizer(actdis,optimizationParams));

  return;

}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TOPOPT::Optimizer& ADAPTER::TopOptBaseAlgorithm::Optimizer()
{
  return *optimizer_;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::TopOptBaseAlgorithm::~TopOptBaseAlgorithm()
{
}

#endif  // #ifdef CCADISCRET
