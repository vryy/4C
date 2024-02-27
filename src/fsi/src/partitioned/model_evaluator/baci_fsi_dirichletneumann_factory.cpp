/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate DirichletNeumann Algorithm


\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_fsi_dirichletneumann_factory.hpp"

#include "baci_fsi_dirichletneumann.hpp"
#include "baci_fsi_dirichletneumann_disp.hpp"
#include "baci_fsi_dirichletneumann_vel.hpp"
#include "baci_fsi_dirichletneumann_volcoupl.hpp"
#include "baci_fsi_dirichletneumannslideale.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_fsi.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FSI::DirichletNeumann> FSI::DirichletNeumannFactory::CreateAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& fsidyn)
{
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  INPAR::FSI::PartitionedCouplingMethod method =
      CORE::UTILS::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsipart, "PARTITIONED");
  switch (method)
  {
    case INPAR::FSI::DirichletNeumannSlideale:
      switch (GLOBAL::Problem::Instance()->GetProblemType())
      {
        case (GLOBAL::ProblemType::fsi):
        case (GLOBAL::ProblemType::fsi_redmodels):
        case (GLOBAL::ProblemType::fsi_lung):
          if (CORE::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
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
      switch (GLOBAL::Problem::Instance()->GetProblemType())
      {
        case (GLOBAL::ProblemType::fsi):
        case (GLOBAL::ProblemType::fsi_redmodels):
        case (GLOBAL::ProblemType::fsi_lung):
          if (CORE::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
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
      switch (GLOBAL::Problem::Instance()->GetProblemType())
      {
        case (GLOBAL::ProblemType::fsi):
        case (GLOBAL::ProblemType::fsi_redmodels):
        case (GLOBAL::ProblemType::fsi_lung):
        case (GLOBAL::ProblemType::fsi_xfem):
          if (CORE::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              INPAR::FSI::CoupVarPart::vel)
          {
            dserror(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannDisp(comm));
          break;
        case (GLOBAL::ProblemType::fbi):
          if (CORE::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
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

BACI_NAMESPACE_CLOSE
