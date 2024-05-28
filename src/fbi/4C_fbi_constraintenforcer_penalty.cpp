/*----------------------------------------------------------------------*/
/*! \file

\brief Implements the constraint enforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction)

\level 2

*----------------------------------------------------------------------*/

#include "4C_fbi_constraintenforcer_penalty.hpp"

#include "4C_adapter_fld_fbi_movingboundary.hpp"
#include "4C_adapter_str_fbiwrapper.hpp"
#include "4C_fbi_adapter_constraintbridge_penalty.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_output_params.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fbi_constraintenforcer.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"

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
  if ((get_discretizations()[1]->Comm().MyPID() == 0) &&
      (bridge()
              ->GetParams()
              ->get_visualization_ouput_params_ptr()
              ->get_constraint_violation_output_flag()))
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
ADAPTER::FBIPenaltyConstraintenforcer::assemble_fluid_coupling_matrix() const
{
  // Get coupling contributions to the fluid stiffness matrix

  return bridge()->GetCff();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseMatrix>
ADAPTER::FBIPenaltyConstraintenforcer::assemble_structure_coupling_matrix() const
{
  // For the classical partitioned algorithm we do not have any contributions to the stiffness
  // matrix of the structure field
  return Teuchos::null;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
ADAPTER::FBIPenaltyConstraintenforcer::assemble_fluid_coupling_residual() const
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(bridge(), true)
      ->scale_penalty_fluid_contributions();
  // Get the force acting on the fluid field, scale it with -1 to get the
  // correct direction
  Teuchos::RCP<Epetra_Vector> f =
      Teuchos::rcp(new Epetra_Vector((bridge()->get_fluid_coupling_residual())->Map()));
  f->Update(-1.0, *(bridge()->get_fluid_coupling_residual()), 0.0);
  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
ADAPTER::FBIPenaltyConstraintenforcer::assemble_structure_coupling_residual() const
{
  Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(bridge(), true)
      ->scale_penalty_structure_contributions();
  // Get the force acting on the structure field, scale it with the penalty factor and -1 to get the
  // correct direction
  Teuchos::RCP<Epetra_Vector> f =
      Teuchos::rcp(new Epetra_Vector(bridge()->get_structure_coupling_residual()->Map()));
  f->Update(-1.0, *(bridge()->get_structure_coupling_residual()), 0.0);

  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::PrepareFluidSolve() { bridge()->PrepareFluidSolve(); }
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::Output(double time, int step)
{
  print_violation(time, step);
}
/*----------------------------------------------------------------------*/

void ADAPTER::FBIPenaltyConstraintenforcer::print_violation(double time, int step)
{
  if (bridge()
          ->GetParams()
          ->get_visualization_ouput_params_ptr()
          ->get_constraint_violation_output_flag())
  {
    double penalty_parameter = bridge()->GetParams()->GetPenaltyParameter();

    Teuchos::RCP<Epetra_Vector> violation = CORE::LINALG::CreateVector(
        Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(get_fluid(), true)->Velnp()->Map());

    int err =
        Teuchos::rcp_dynamic_cast<const ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
            ->GetCff()
            ->Multiply(false,
                *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(get_fluid(), true)->Velnp()),
                *violation);

    if (err != 0) FOUR_C_THROW(" Matrix vector product threw error code %i ", err);

    err = violation->Update(1.0, *assemble_fluid_coupling_residual(), -1.0);
    if (err != 0) FOUR_C_THROW(" Epetra_Vector update threw error code %i ", err);

    double norm = 0.0, normf = 0.0, norms = 0.0, norm_vel = 0.0;

    get_velocity_pressure_splitter()
        ->ExtractOtherVector(
            Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(get_fluid(), true)->Velnp())
        ->MaxValue(&norm_vel);

    violation->MaxValue(&norm);
    if (norm_vel > 1e-15) normf = norm / norm_vel;

    Teuchos::rcp_dynamic_cast<const ADAPTER::FBIStructureWrapper>(GetStructure(), true)
        ->Velnp()
        ->MaxValue(&norm_vel);
    if (norm_vel > 1e-15) norms = norm / norm_vel;

    std::ofstream log;
    if (get_discretizations()[1]->Comm().MyPID() == 0)
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
