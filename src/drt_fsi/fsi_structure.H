
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_STRUCTURE_H
#define FSI_STRUCTURE_H

#include "../drt_lib/strugenalpha.H"


namespace FSI
{

/// Structural dynamics algorithm enhanced for FSI
class Structure : public StruGenAlpha
{
public:
  Structure(ParameterList& params,
            DRT::Discretization& dis,
            LINALG::Solver& solver,
            DiscretizationWriter& output);
};

}

#endif
#endif
#endif
