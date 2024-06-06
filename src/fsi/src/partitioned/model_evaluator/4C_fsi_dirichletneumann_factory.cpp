/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate DirichletNeumann Algorithm


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_dirichletneumann_factory.hpp"

#include "4C_fsi_dirichletneumann.hpp"
#include "4C_fsi_dirichletneumann_disp.hpp"
#include "4C_fsi_dirichletneumann_vel.hpp"
#include "4C_fsi_dirichletneumann_volcoupl.hpp"
#include "4C_fsi_dirichletneumannslideale.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FSI::DirichletNeumann> FSI::DirichletNeumannFactory::CreateAlgorithm(
    const Epetra_Comm& comm, const Teuchos::ParameterList& fsidyn)
{
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  Inpar::FSI::PartitionedCouplingMethod method =
      Core::UTILS::IntegralValue<Inpar::FSI::PartitionedCouplingMethod>(fsipart, "PARTITIONED");
  switch (method)
  {
    case Inpar::FSI::DirichletNeumannSlideale:
      switch (Global::Problem::Instance()->GetProblemType())
      {
        case (Core::ProblemType::fsi):
        case (Core::ProblemType::fsi_redmodels):
        case (Core::ProblemType::fsi_lung):
          if (Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              Inpar::FSI::CoupVarPart::vel)
          {
            FOUR_C_THROW(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannSlideale(comm));
          break;
        default:
          FOUR_C_THROW("Your problem does not work with DirichletNeumann Slide ALE yet!!");
          break;
      }
      break;
    case Inpar::FSI::DirichletNeumannVolCoupl:
      switch (Global::Problem::Instance()->GetProblemType())
      {
        case (Core::ProblemType::fsi):
        case (Core::ProblemType::fsi_redmodels):
        case (Core::ProblemType::fsi_lung):
          if (Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              Inpar::FSI::CoupVarPart::vel)
          {
            FOUR_C_THROW(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannVolCoupl(comm));
          break;
        default:
          FOUR_C_THROW("Your problem does not work with DirichletNeumann Volume Coupling yet!!");
          break;
      }
      break;
    case Inpar::FSI::DirichletNeumann:
      switch (Global::Problem::Instance()->GetProblemType())
      {
        case (Core::ProblemType::fsi):
        case (Core::ProblemType::fsi_redmodels):
        case (Core::ProblemType::fsi_lung):
        case (Core::ProblemType::fsi_xfem):
          if (Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              Inpar::FSI::CoupVarPart::vel)
          {
            FOUR_C_THROW(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface velocities. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannDisp(comm));
          break;
        case (Core::ProblemType::fbi):
          if (Core::UTILS::IntegralValue<int>(fsipart, "COUPVARIABLE") ==
              Inpar::FSI::CoupVarPart::disp)
          {
            FOUR_C_THROW(
                "Displacement coupling is not possible in this case! You are not handling any "
                "interface displacements. Check your problem type!");
          }
          else
            return Teuchos::rcp(new FSI::DirichletNeumannVel(comm));
          break;
        default:
          FOUR_C_THROW("Your problem does not work with DirichletNeumann yet!!");
          break;
      }
      break;
  }
  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
