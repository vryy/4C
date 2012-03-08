/*!------------------------------------------------------------------------------------------------*
\file adapter_fluid_topopt_algo.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET


#include "adapter_fluid_topopt_algo.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidTopOptCouplingAlgorithm::FluidTopOptCouplingAlgorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn
    )
:  FluidBaseAlgorithm(prbdyn,false),
   TopOptBaseAlgorithm(prbdyn,0),
   TopOptFluidAdjointAlgorithm(prbdyn),
   params_(prbdyn)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidTopOptCouplingAlgorithm::~FluidTopOptCouplingAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidTopOptCouplingAlgorithm::ReadRestart(int step)
{dserror("change");
  FluidField().ReadRestart(step);
  return;
}

#endif  // #ifdef CCADISCRET
