/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for time step size adaptivity within monolithic FSI

\level 2


*/

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_adapter_str_fsi_timint_adaptive.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io_pstream.hpp"
#include "4C_structure_aux.hpp"
#include "4C_structure_timada.hpp"
#include "4C_structure_timint.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/*======================================================================*/
/* constructor */
Adapter::StructureFSITimIntAda::StructureFSITimIntAda(
    Teuchos::RCP<Solid::TimAda> sta, Teuchos::RCP<Structure> sti)
    : FSIStructureWrapper(sti), StructureTimIntAda(sta, sti), str_time_integrator_(sti)
{
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();
  const Teuchos::ParameterList& sada = sdyn.sublist("TIMEADAPTIVITY");

  // type of error norm
  errnorm_ = Core::UTILS::IntegralValue<Inpar::Solid::VectorNorm>(sada, "LOCERRNORM");

  //----------------------------------------------------------------------------
  // Handling of Dirichlet BCs in error estimation
  //----------------------------------------------------------------------------
  // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface.
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(sti->get_dbc_map_extractor()->cond_map());
  intersectionmaps.push_back(Interface()->fsi_cond_map());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::intersect_maps(intersectionmaps);

  numdbcdofs_ = sti->get_dbc_map_extractor()->cond_map()->NumGlobalElements();
  numdbcfsidofs_ = intersectionmap->NumGlobalElements();
  numdbcinnerdofs_ = numdbcdofs_ - numdbcfsidofs_;
}

/*----------------------------------------------------------------------------*/
/* Indicate norms of local discretization error */
void Adapter::StructureFSITimIntAda::IndicateErrorNorms(double& err, double& errcond,
    double& errother, double& errinf, double& errinfcond, double& errinfother)
{
  // call functionality of adaptive structural time integrator
  str_ada()->evaluate_local_error_dis();

  // Indicate has to be done by the FSI algorithm since it depends on the interface
  indicate_errors(err, errcond, errother, errinf, errinfcond, errinfother);

  return;
}

/*----------------------------------------------------------------------------*/
/* Indicate local discretization error */
void Adapter::StructureFSITimIntAda::indicate_errors(double& err, double& errcond, double& errother,
    double& errinf, double& errinfcond, double& errinfother)
{
  // vector with local discretization error for each DOF
  Teuchos::RCP<Epetra_Vector> error = str_ada()->LocErrDis();

  // extract the condition part of the full error vector
  // (i.e. only interface displacement DOFs)
  Teuchos::RCP<Epetra_Vector> errorcond =
      Teuchos::rcp(new Epetra_Vector(*interface_->extract_fsi_cond_vector(error)));

  // in case of structure split: extract the other part of the full error vector
  // (i.e. only interior displacement DOFs)
  Teuchos::RCP<Epetra_Vector> errorother =
      Teuchos::rcp(new Epetra_Vector(*interface_->extract_fsi_cond_vector(error)));

  // calculate L2-norms of different subsets of local discretization error vector
  err = Solid::calculate_vector_norm(errnorm_, error, numdbcdofs_);
  errcond = Solid::calculate_vector_norm(errnorm_, errorcond, numdbcfsidofs_);
  errother = Solid::calculate_vector_norm(errnorm_, errorother, numdbcinnerdofs_);

  // calculate L-inf-norms of different subsets of local discretization error vector
  errinf = Solid::calculate_vector_norm(Inpar::Solid::norm_inf, error);
  errinfcond = Solid::calculate_vector_norm(Inpar::Solid::norm_inf, errorcond);
  errinfother = Solid::calculate_vector_norm(Inpar::Solid::norm_inf, errorother);

  return;
}

/*----------------------------------------------------------------------------*/
/* Do a single step with auxiliary time integration scheme */
void Adapter::StructureFSITimIntAda::time_step_auxiliar() { str_ada()->integrate_step_auxiliar(); }

/*----------------------------------------------------------------------------*/
/* Calculate time step size suggestion */
double Adapter::StructureFSITimIntAda::calculate_dt(const double norm)
{
  return str_ada()->calculate_dt(norm);
}

/*----------------------------------------------------------------------------*/
/* Get time step size of adaptive structural time integrator */
double Adapter::StructureFSITimIntAda::Dt() const { return str_ada()->Dt(); }

/*----------------------------------------------------------------------------*/
/* Get target time \f$t_{n+1}\f$ of current time step */
double Adapter::StructureFSITimIntAda::Time() const { return str_ada()->Time(); }

/*----------------------------------------------------------------------------*/
/* Set new time step size */
void Adapter::StructureFSITimIntAda::set_dt(const double dtnew) { str_ada()->set_dt(dtnew); }

/*----------------------------------------------------------------------------*/
/* Update step size */
void Adapter::StructureFSITimIntAda::UpdateStepSize(const double dtnew)
{
  str_ada()->UpdateStepSize(dtnew);
}

/*----------------------------------------------------------------------------*/
/*  Reset certain quantities to prepare repetition of current time step */
void Adapter::StructureFSITimIntAda::reset_step() { str_ada()->reset_step(); }

FOUR_C_NAMESPACE_CLOSE
