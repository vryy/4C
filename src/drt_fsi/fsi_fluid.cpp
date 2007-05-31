
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_fluid.H"

FSI::Fluid::Fluid(RefCountPtr<DRT::Discretization> dis,
                  LINALG::Solver&       solver,
                  ParameterList&        params,
                  DiscretizationWriter& output)
  : FluidImplicitTimeInt(dis,solver,params,output)
{
}

#endif
#endif
