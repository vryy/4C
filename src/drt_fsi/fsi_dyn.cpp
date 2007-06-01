
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"

void fsi_ale_drt()
{
  FSI::DirichletNeumannCoupling fsi;

  fsi.Setup();
  fsi.Timeloop();
}

#endif
#endif
