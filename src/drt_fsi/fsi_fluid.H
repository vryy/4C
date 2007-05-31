
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_FLUID_H
#define FSI_FLUID_H

#include "../drt_lib/fluidimplicitintegration.H"

namespace FSI {

/// Fluid algorithm enhanced for FSI
class Fluid : public FluidImplicitTimeInt
{
public:
  Fluid(RefCountPtr<DRT::Discretization> dis,
        LINALG::Solver&       solver,
        ParameterList&        params,
        DiscretizationWriter& output);

};

}

#endif
#endif
#endif
