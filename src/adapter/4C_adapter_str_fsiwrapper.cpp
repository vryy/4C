/*----------------------------------------------------------------------*/
/*! \file

\brief Structural adapter for FSI problems containing the interface
       and methods dependent on the interface


\level 1
*/

#include "4C_adapter_str_fsiwrapper.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fsi_str_model_evaluator_partitioned.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_structure_aux.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  bool prestress_is_active(const double currentTime)
  {
    Inpar::Solid::PreStress pstype = Teuchos::getIntegralValue<Inpar::Solid::PreStress>(
        Global::Problem::instance()->structural_dynamic_params(), "PRESTRESS");
    const double pstime =
        Global::Problem::instance()->structural_dynamic_params().get<double>("PRESTRESSTIME");
    return pstype != Inpar::Solid::PreStress::none && currentTime <= pstime + 1.0e-15;
  }
}  // namespace

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::FSIStructureWrapper::FSIStructureWrapper(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure)
{
  // set-up FSI interface
  interface_ = Teuchos::rcp(new Solid::MapExtractor);

  if (Global::Problem::instance()->get_problem_type() != Core::ProblemType::fpsi)
    interface_->setup(*discretization(), *discretization()->dof_row_map());
  else
    interface_->setup(*discretization(), *discretization()->dof_row_map(),
        true);  // create overlapping maps for fpsi problem

  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");
  predictor_ = Core::UTILS::integral_value<int>(fsipart, "PREDICTOR");
}

/*------------------------------------------------------------------------------------*
 * Rebuild FSI interface on structure side                              sudhakar 09/13
 * This is necessary if elements are added/deleted from interface
 *------------------------------------------------------------------------------------*/
void Adapter::FSIStructureWrapper::rebuild_interface()
{
  interface_ = Teuchos::rcp(new Solid::MapExtractor);
  interface_->setup(*discretization(), *discretization()->dof_row_map());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FSIStructureWrapper::use_block_matrix()
{
  StructureWrapper::use_block_matrix(interface_, interface_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FSIStructureWrapper::relaxation_solve(
    Teuchos::RCP<Epetra_Vector> iforce)
{
  apply_interface_forces(iforce);
  fsi_model_evaluator()->set_is_relaxation_solve(true);
  Teuchos::RCP<const Epetra_Vector> idisi =
      fsi_model_evaluator()->solve_relaxation_linear(structure_);
  fsi_model_evaluator()->set_is_relaxation_solve(false);

  // we are just interested in the incremental interface displacements
  return interface_->extract_fsi_cond_vector(idisi);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FSIStructureWrapper::predict_interface_dispnp()
{
  // prestressing business
  Teuchos::RCP<Epetra_Vector> idis;

  switch (predictor_)
  {
    case 1:
    {
      // d(n)
      // respect Dirichlet conditions at the interface (required for pseudo-rigid body)
      if (prestress_is_active(time()))
      {
        idis = Teuchos::rcp(new Epetra_Vector(*interface_->fsi_cond_map(), true));
      }
      else
      {
        idis = interface_->extract_fsi_cond_vector(dispn());
      }
      break;
    }
    case 2:
      // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
      FOUR_C_THROW("interface velocity v(n-1) not available");
      break;
    case 3:
    {
      // d(n)+dt*v(n)
      if (prestress_is_active(time()))
        FOUR_C_THROW("only constant interface predictor useful for prestressing");

      double current_dt = dt();

      idis = interface_->extract_fsi_cond_vector(dispn());
      Teuchos::RCP<Epetra_Vector> ivel = interface_->extract_fsi_cond_vector(veln());

      idis->Update(current_dt, *ivel, 1.0);
      break;
    }
    case 4:
    {
      // d(n)+dt*v(n)+0.5*dt^2*a(n)
      if (prestress_is_active(time()))
        FOUR_C_THROW("only constant interface predictor useful for prestressing");

      double current_dt = dt();

      idis = interface_->extract_fsi_cond_vector(dispn());
      Teuchos::RCP<Epetra_Vector> ivel = interface_->extract_fsi_cond_vector(veln());
      Teuchos::RCP<Epetra_Vector> iacc = interface_->extract_fsi_cond_vector(accn());

      idis->Update(current_dt, *ivel, 0.5 * current_dt * current_dt, *iacc, 1.0);
      break;
    }
    default:
      FOUR_C_THROW(
          "unknown interface displacement predictor '%s'", Global::Problem::instance()
                                                               ->fsi_dynamic_params()
                                                               .sublist("PARTITIONED SOLVER")
                                                               .get<std::string>("PREDICTOR")
                                                               .c_str());
      break;
  }

  return idis;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FSIStructureWrapper::extract_interface_dispn()
{
  FOUR_C_ASSERT(interface_->full_map()->PointSameAs(dispn()->Map()),
      "Full map of map extractor and Dispn() do not match.");

  // prestressing business
  if (prestress_is_active(time()))
  {
    return Teuchos::rcp(new Epetra_Vector(*interface_->fsi_cond_map(), true));
  }
  else
  {
    return interface_->extract_fsi_cond_vector(dispn());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> Adapter::FSIStructureWrapper::extract_interface_dispnp()
{
  FOUR_C_ASSERT(interface_->full_map()->PointSameAs(dispnp()->Map()),
      "Full map of map extractor and Dispnp() do not match.");

  // prestressing business
  if (prestress_is_active(time()))
  {
    if (discretization()->get_comm().MyPID() == 0)
      std::cout << "Applying no displacements to the fluid since we do prestressing" << std::endl;

    return Teuchos::rcp(new Epetra_Vector(*interface_->fsi_cond_map(), true));
  }
  else
  {
    return interface_->extract_fsi_cond_vector(dispnp());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces
void Adapter::FSIStructureWrapper::apply_interface_forces(Teuchos::RCP<Epetra_Vector> iforce)
{
  fsi_model_evaluator()->get_interface_force_np_ptr()->PutScalar(0.0);
  interface_->add_fsi_cond_vector(iforce, fsi_model_evaluator()->get_interface_force_np_ptr());
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Apply interface forces deprecated version ! Remove as soon as possible!
void Adapter::FSIStructureWrapper::apply_interface_forces_temporary_deprecated(
    Teuchos::RCP<Epetra_Vector> iforce)
{
  Teuchos::RCP<Epetra_Vector> fifc = Core::LinAlg::create_vector(*dof_row_map(), true);

  interface_->add_fsi_cond_vector(iforce, fifc);

  set_force_interface(fifc);

  prepare_partition_step();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Solid::ModelEvaluator::PartitionedFSI>
Adapter::FSIStructureWrapper::fsi_model_evaluator()
{
  return fsi_model_evaluator_;
};

FOUR_C_NAMESPACE_CLOSE
