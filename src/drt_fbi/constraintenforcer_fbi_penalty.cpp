/*----------------------------------------------------------------------*/
/*! \file

\brief Implements the constraint enforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction)

\level 2

*----------------------------------------------------------------------*/

#include "ad_fbi_constraintbridge_penalty.H"
#include "beam_to_fluid_meshtying_params.H"
#include "beam_to_fluid_meshtying_vtk_output_params.H"
#include "constraintenforcer_fbi_penalty.H"
#include "constraintenforcer_fbi.H"

#include "../drt_adapter/ad_fld_fbi_movingboundary.H"
#include "../drt_adapter/ad_str_fbiwrapper.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::Setup(
    Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid)
{
  ADAPTER::FBIConstraintenforcer::Setup(structure, fluid);
  std::ofstream log;
  if ((GetDiscretizations()[1]->Comm().MyPID() == 0) &&
      (Bridge()->GetParams()->GetVtkOuputParamsPtr()->GetConstraintViolationOutputFlag()))
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".penalty");
    log.open(s.c_str(), std::ofstream::out);
    log << "Time \t Step \t ViolationNorm \t FluidViolationNorm \t StructureViolationNorm"
        << std::endl;
    log.close();
  }
}

Teuchos::RCP<const LINALG::SparseMatrix>
ADAPTER::FBIPenaltyConstraintenforcer::AssembleFluidCouplingMatrix() const
{
  // Get coupling contributions to the fluid stiffness matrix

  return Bridge()->GetCff();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::SparseMatrix>
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
  if (Bridge()->GetParams()->GetVtkOuputParamsPtr()->GetConstraintViolationOutputFlag())
  {
    double penalty_parameter = Bridge()->GetParams()->GetPenaltyParameter();

    Teuchos::RCP<Epetra_Vector> violation = LINALG::CreateVector(
        Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(GetFluid(), true)->Velnp()->Map());

    int err =
        Teuchos::rcp_dynamic_cast<const ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
            ->GetCff()
            ->Multiply(false,
                *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(GetFluid(), true)->Velnp()),
                *violation);

    if (err != 0) dserror(" Matrix vector product threw error code %i ", err);

    err = violation->Update(1.0, *AssembleFluidCouplingResidual(), -1.0);
    if (err != 0) dserror(" Epetra_Vector update threw error code %i ", err);

    double norm, normf, norms;
    double norm_vel;

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
      std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
      s.append(".penalty");
      log.open(s.c_str(), std::ofstream::app);
      log << time << "\t" << step << "\t" << norm / penalty_parameter << "\t"
          << normf / penalty_parameter << "\t" << norms / penalty_parameter << std::endl;
    }
  }
}
