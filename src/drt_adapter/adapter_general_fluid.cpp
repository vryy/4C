
#ifdef CCADISCRET

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include "adapter_general_fluid.H"
#include "adapter_fluid_ale.H"




/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::GeneralFluidBaseAlgorithm::GeneralFluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn,
                                                              std::string condname)
{
  // here we could do some decision what kind of generalized fluid to build
  fluid_ = Teuchos::rcp(new FluidAleAdapter(prbdyn,condname));
}


#endif
