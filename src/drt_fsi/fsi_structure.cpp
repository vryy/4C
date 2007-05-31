
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_structure.H"

FSI::Structure::Structure(ParameterList& params,
                          DRT::Discretization& dis,
                          LINALG::Solver& solver,
                          DiscretizationWriter& output)
  : StruGenAlpha(params, dis, solver, output)
{
}

#endif
#endif
