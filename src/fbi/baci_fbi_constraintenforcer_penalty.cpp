/*----------------------------------------------------------------------*/
/*! \file

\brief Implements the constraint enforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction)

\level 2

*----------------------------------------------------------------------*/

#include "baci_fbi_constraintenforcer_penalty.hpp"

#include "baci_adapter_fld_fbi_movingboundary.hpp"
#include "baci_adapter_str_fbiwrapper.hpp"
#include "baci_fbi_adapter_constraintbridge_penalty.hpp"
#include "baci_fbi_beam_to_fluid_meshtying_output_params.hpp"
#include "baci_fbi_beam_to_fluid_meshtying_params.hpp"
#include "baci_fbi_constraintenforcer.hpp"
#include "baci_global_data.hpp"
#include "baci_io_control.hpp"
#include "baci_linalg_mapextractor.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_linalg_sparseoperator.hpp"
#include "baci_linalg_utils_sparse_algebra_create.hpp"

#include <Epetra_Vector.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::Setup(
    Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid)
{
  ADAPTER::FBIConstraintenforcer::Setup(structure, fluid);
  std::ofstream log;
  if ((GetDiscretizations()[1]->Comm().MyPID() == 0) &&
      (Bridge()->GetParams()->GetVisualizationOuputParamsPtr()->GetConstraintViolationOutputFlag()))
  {
    std::string s = GLOBAL::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".penalty");
    log.open(s.c_str(), std::ofstream::out);
    log << "Time \t Step \t ViolationNorm \t FluidViolationNorm \t StructureViolationNorm"
        << std::endl;
    log.close();
  }
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseOperator>
ADAPTER::FBIPenaltyConstraintenforcer::AssembleFluidCouplingMatrix() const
{
  // Get coupling contributions to the fluid stiffness matrix

  return Bridge()->GetCff();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseMatrix>
ADAPTER::FBIPenaltyConstraintenforcer::AssembleStructureCouplingMatrix() const
{
  // For the classical partitioned algorithm we do not have any contributions to the stiffness
  // matrix of the structure field
  return Teuchos::null;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIPenaltyConstraintenforcer::AssembleFluidCouplingResidual()
    const
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(Bridge(), true)
      ->ScalePenaltyFluidContributions();
  // Get the force acting on the fluid field, scale it with -1 to get the
  // correct direction
  Teuchos::RCP<Epetra_Vector> f =
      Teuchos::rcp(new Epetra_Vector((Bridge()->GetFluidCouplingResidual())->Map()));
  f->Update(-1.0, *(Bridge()->GetFluidCouplingResidual()), 0.0);
  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
ADAPTER::FBIPenaltyConstraintenforcer::AssembleStructureCouplingResidual() const
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(Bridge(), true)
      ->ScalePenaltyStructureContributions();
  // Get the force acting on the structure field, scale it with the penalty factor and -1 to get the
  // correct direction
  Teuchos::RCP<Epetra_Vector> f =
      Teuchos::rcp(new Epetra_Vector(Bridge()->GetStructureCouplingResidual()->Map()));
  f->Update(-1.0, *(Bridge()->GetStructureCouplingResidual()), 0.0);

  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::PrepareFluidSolve() { Bridge()->PrepareFluidSolve(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::Output(double time, int step)
{
  PrintViolation(time, step);
}
/*----------------------------------------------------------------------*/

void ADAPTER::FBIPenaltyConstraintenforcer::PrintViolation(double time, int step)
{
  if (Bridge()->GetParams()->GetVisualizationOuputParamsPtr()->GetConstraintViolationOutputFlag())
  {
    double penalty_parameter = Bridge()->GetParams()->GetPenaltyParameter();

    Teuchos::RCP<Epetra_Vector> violation = CORE::LINALG::CreateVector(
        Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(GetFluid(), true)->Velnp()->Map());

    int err =
        Teuchos::rcp_dynamic_cast<const ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
            ->GetCff()
            ->Multiply(false,
                *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(GetFluid(), true)->Velnp()),
                *violation);

    if (err != 0) FOUR_C_THROW(" Matrix vector product threw error code %i ", err);

    err = violation->Update(1.0, *AssembleFluidCouplingResidual(), -1.0);
    if (err != 0) FOUR_C_THROW(" Epetra_Vector update threw error code %i ", err);

    double norm = 0.0, normf = 0.0, norms = 0.0, norm_vel = 0.0;

    GetVelocityPressureSplitter()
        ->ExtractOtherVector(
            Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(GetFluid(), true)->Velnp())
        ->MaxValue(&norm_vel);

    violation->MaxValue(&norm);
    if (norm_vel > 1e-15) normf = norm / norm_vel;

    Teuchos::rcp_dynamic_cast<const ADAPTER::FBIStructureWrapper>(GetStructure(), true)
        ->Velnp()
        ->MaxValue(&norm_vel);
    if (norm_vel > 1e-15) norms = norm / norm_vel;

    std::ofstream log;
    if (GetDiscretizations()[1]->Comm().MyPID() == 0)
    {
      std::string s = GLOBAL::Problem::Instance()->OutputControlFile()->FileName();
      s.append(".penalty");
      log.open(s.c_str(), std::ofstream::app);
      log << time << "\t" << step << "\t" << norm / penalty_parameter << "\t"
          << normf / penalty_parameter << "\t" << norms / penalty_parameter << std::endl;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
