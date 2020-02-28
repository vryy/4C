/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate DirichletNeumann Algorithm

\maintainer Nora Hagmeyer

\level 1
*/
/*----------------------------------------------------------------------*/

#include "fsi_dirichletneumann.H"
#include "fsi_dirichletneumann_factory.H"
#include "fsi_dirichletneumann_disp.H"
#include "fsi_dirichletneumann_vel.H"
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
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsipart, "PARTITIONED");
  switch (method)
  {
    case INPAR::FSI::DirichletNeumannSlideale:
      switch (DRT::Problem::Instance()->GetProblemType())
      {
        case (prb_fsi):
        case (prb_fsi_redmodels):
        case (prb_fsi_lung):
          if (DRT::INPUT::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              INPAR::FSI::CoupVarPart::vel)
          {
            dserror(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannSlideale(comm));
          break;
        default:
          dserror("Your problem does not work with DirichletNeumann Slide ALE yet!!");
          break;
      }
      break;
    case INPAR::FSI::DirichletNeumannVolCoupl:
      switch (DRT::Problem::Instance()->GetProblemType())
      {
        case (prb_fsi):
        case (prb_fsi_redmodels):
        case (prb_fsi_lung):
          if (DRT::INPUT::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              INPAR::FSI::CoupVarPart::vel)
          {
            dserror(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannVolCoupl(comm));
          break;
        default:
          dserror("Your problem does not work with DirichletNeumann Volume Coupling yet!!");
          break;
      }
      break;
    case INPAR::FSI::DirichletNeumann:
      switch (DRT::Problem::Instance()->GetProblemType())
      {
        case (prb_fsi):
        case (prb_fsi_redmodels):
        case (prb_fsi_lung):
        case (prb_fsi_xfem):
          if (DRT::INPUT::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              INPAR::FSI::CoupVarPart::vel)
          {
            dserror(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannDisp(comm));
          break;
        case (prb_fbi):
          if (DRT::INPUT::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              INPAR::FSI::CoupVarPart::disp)
          {
            dserror(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface displacements. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannVel(comm));
          break;
        default:
          dserror("Your problem does not work with DirichletNeumann yet!!");
          break;
      }
      break;
  }
  return Teuchos::null;
}
