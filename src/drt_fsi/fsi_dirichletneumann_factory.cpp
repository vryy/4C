/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann_factory.cpp

\brief Factory to create the appropriate DirichletNeumann Algorithm

\maintainer Nora Hagmeyer

\level 1
*/
/*----------------------------------------------------------------------*/


#include "fsi_dirichletneumann_factory.H"
#include "fsi_dirichletneumann_disp.H"
#include "fsi_dirichletneumann_volcoupl.H"
#include "fsi_dirichletneumannslideale.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_inpar/inpar_fsi.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FSI::DirichletNeumann> FSI::DirichletNeumannFactory::CreateAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& fsidyn)
{
  INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(
          fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");
  switch (method)
  {
    case INPAR::FSI::DirichletNeumannSlideale:
      switch (DRT::Problem::Instance()->ProblemType())
      {
        case (prb_fsi):
        case (prb_fsi_redmodels):
        case (prb_fsi_lung):
          return Teuchos::rcp(new FSI::DirichletNeumannSlideale(comm));
          break;
        default:
          dserror("Your problem does not work with DirichletNeumann Slide ALE yet!!");
          break;
      }
      break;
    case INPAR::FSI::DirichletNeumannVolCoupl:
      switch (DRT::Problem::Instance()->ProblemType())
      {
        case (prb_fsi):
        case (prb_fsi_redmodels):
        case (prb_fsi_lung):
          return Teuchos::rcp(new FSI::DirichletNeumannVolCoupl(comm));
          break;
        default:
          dserror("Your problem does not work with DirichletNeumann Volume Coupling yet!!");
          break;
      }
      break;
    case INPAR::FSI::DirichletNeumann:
      switch (DRT::Problem::Instance()->ProblemType())
      {
        case (prb_fsi):
        case (prb_fsi_redmodels):
        case (prb_fsi_lung):
        case (prb_fsi_xfem):
          return Teuchos::rcp(new FSI::DirichletNeumannDisp(comm));
          break;
        default:
          dserror("Your problem does not work with DirichletNeumann yet!!");
          break;
      }
      break;
  }
  return Teuchos::null;
}
