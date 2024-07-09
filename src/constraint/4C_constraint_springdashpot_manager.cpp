/*----------------------------------------------------------------------*/
/*! \file

\brief Methods for spring and dashpot constraints / boundary conditions:

\level 2


*----------------------------------------------------------------------*/

#include "4C_constraint_springdashpot_manager.hpp"

#include "4C_constraint_springdashpot.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_io.hpp"
#include "4C_io_pstream.hpp"  // has to go before io.hpp

#include <iostream>

FOUR_C_NAMESPACE_OPEN

CONSTRAINTS::SpringDashpotManager::SpringDashpotManager(Teuchos::RCP<Core::FE::Discretization> dis)
    : actdisc_(dis), havespringdashpot_(false)
{
  // get all spring dashpot conditions
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> springdashpots;
  actdisc_->get_condition("RobinSpringDashpot", springdashpots);

  // number of spring dashpot conditions
  n_conds_ = (int)springdashpots.size();

  // check if spring dashpots are present
  if (n_conds_)
  {
    // yes! we have spring dashpot conditions!
    havespringdashpot_ = true;

    // new instance of spring dashpot BC with current condition for every spring dashpot condition
    for (int i = 0; i < n_conds_; ++i)
      springs_.push_back(Teuchos::rcp(new SpringDashpot(actdisc_, springdashpots[i])));
  }

  return;
}

void CONSTRAINTS::SpringDashpotManager::stiffness_and_internal_forces(
    Teuchos::RCP<Core::LinAlg::SparseMatrix> stiff, Teuchos::RCP<Epetra_Vector> fint,
    Teuchos::RCP<Epetra_Vector> disn, Teuchos::RCP<Epetra_Vector> veln,
    Teuchos::ParameterList parlist)
{
  // evaluate all spring dashpot conditions
  for (int i = 0; i < n_conds_; ++i)
  {
    springs_[i]->reset_newton();
    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = springs_[i]->get_spring_type();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
      springs_[i]->evaluate_robin(stiff, fint, disn, veln, parlist);
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
      springs_[i]->evaluate_force_stiff(*stiff, *fint, disn, veln, parlist);
  }

  return;
}

void CONSTRAINTS::SpringDashpotManager::update()
{
  // update all spring dashpot conditions for each new time step
  for (int i = 0; i < n_conds_; ++i) springs_[i]->update();

  return;
}

void CONSTRAINTS::SpringDashpotManager::reset_prestress(Teuchos::RCP<Epetra_Vector> dis)
{
  // loop over all spring dashpot conditions and reset them
  for (int i = 0; i < n_conds_; ++i) springs_[i]->reset_prestress(dis);

  return;
}

void CONSTRAINTS::SpringDashpotManager::output(Teuchos::RCP<Core::IO::DiscretizationWriter> output,
    Teuchos::RCP<Core::FE::Discretization> discret, Teuchos::RCP<Epetra_Vector> disp)
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> gap =
      Teuchos::rcp(new Epetra_Vector(*(actdisc_->node_row_map()), true));
  Teuchos::RCP<Epetra_MultiVector> normals =
      Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->node_row_map()), 3, true));
  Teuchos::RCP<Epetra_MultiVector> springstress =
      Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->node_row_map()), 3, true));

  // collect outputs from all spring dashpot conditions
  bool found_cursurfnormal = false;
  for (int i = 0; i < n_conds_; ++i)
  {
    springs_[i]->output_gap_normal(gap, normals, springstress);

    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = springs_[i]->get_spring_type();
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal) found_cursurfnormal = true;
  }

  // write vectors to output
  if (found_cursurfnormal)
  {
    output->write_vector("gap", gap);
    output->write_vector("curnormals", normals);
  }

  // write spring stress if defined in io-flag
  if (Core::UTILS::IntegralValue<bool>(Global::Problem::instance()->io_params(), "OUTPUT_SPRING") ==
      true)
    output->write_vector("springstress", springstress);

  return;
}

void CONSTRAINTS::SpringDashpotManager::output_restart(
    Teuchos::RCP<Core::IO::DiscretizationWriter> output_restart,
    Teuchos::RCP<Core::FE::Discretization> discret, Teuchos::RCP<Epetra_Vector> disp)
{
  // row maps for export
  Teuchos::RCP<Epetra_Vector> springoffsetprestr =
      Teuchos::rcp(new Epetra_Vector(*actdisc_->dof_row_map()));
  Teuchos::RCP<Epetra_MultiVector> springoffsetprestr_old =
      Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->node_row_map()), 3, true));

  // collect outputs from all spring dashpot conditions
  for (int i = 0; i < n_conds_; ++i)
  {
    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = springs_[i]->get_spring_type();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
      springs_[i]->output_prestr_offset(springoffsetprestr);
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
      springs_[i]->output_prestr_offset_old(springoffsetprestr_old);
  }

  // write vector to output for restart
  output_restart->write_vector("springoffsetprestr", springoffsetprestr);
  // write vector to output for restart
  output_restart->write_vector("springoffsetprestr_old", springoffsetprestr_old);

  // normal output as well
  output(output_restart, discret, disp);

  return;
}


/*----------------------------------------------------------------------*
|(public)                                                      mhv 03/15|
|Read restart information                                               |
 *-----------------------------------------------------------------------*/
void CONSTRAINTS::SpringDashpotManager::read_restart(
    Core::IO::DiscretizationReader& reader, const double& time)
{
  Teuchos::RCP<Epetra_Vector> tempvec = Teuchos::rcp(new Epetra_Vector(*actdisc_->dof_row_map()));
  Teuchos::RCP<Epetra_MultiVector> tempvecold =
      Teuchos::rcp(new Epetra_MultiVector(*(actdisc_->node_row_map()), 3, true));

  reader.read_vector(tempvec, "springoffsetprestr");
  reader.read_multi_vector(tempvecold, "springoffsetprestr_old");

  // loop over all spring dashpot conditions and set restart
  for (int i = 0; i < n_conds_; ++i)
  {
    // get spring type from current condition
    const CONSTRAINTS::SpringDashpot::SpringType stype = springs_[i]->get_spring_type();

    if (stype == CONSTRAINTS::SpringDashpot::xyz or
        stype == CONSTRAINTS::SpringDashpot::refsurfnormal)
      springs_[i]->set_restart(tempvec);
    if (stype == CONSTRAINTS::SpringDashpot::cursurfnormal)
      springs_[i]->set_restart_old(tempvecold);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
