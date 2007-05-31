
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifndef FSI_DIRICHLETNEUMANN_H
#define FSI_DIRICHLETNEUMANN_H

#include <Teuchos_RefCountPtr.hpp>

#include "fsi_structure.H"
#include "fsi_fluid.H"
#include "fsi_ale.H"

namespace FSI {

/// Dirichlet Neumann coupled interface system
class DirichletNeumannCoupling
{
public:

  /// outer level FSI time loop
  void timeloop();

private:

  RefCountPtr<Structure> structure_;
  RefCountPtr<Fluid>     fluid_;
  RefCountPtr<Ale>       ale_;
};

}

#endif
#endif
#endif
