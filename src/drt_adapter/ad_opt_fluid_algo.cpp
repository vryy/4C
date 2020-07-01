/*----------------------------------------------------------------------*/
/*! \file

\brief fluid - topology optimization adapter

\level 3


 *------------------------------------------------------------------------------------------------*/


#include "ad_opt_fluid_algo.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidTopOptCouplingAlgorithm::FluidTopOptCouplingAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn)
    : FluidBaseAlgorithm(prbdyn, DRT::Problem::Instance()->FluidDynamicParams(), "fluid", false),
      TopOptBaseAlgorithm(prbdyn, "opti"),
      TopOptFluidAdjointAlgorithm(prbdyn),
      params_(prbdyn)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidTopOptCouplingAlgorithm::~FluidTopOptCouplingAlgorithm() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidTopOptCouplingAlgorithm::ReadRestart(int step)
{
  dserror("change");
  FluidField()->ReadRestart(step);
  return;
}
