/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned coupling algorithm for scatra-thermo interaction

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_sti_partitioned.hpp"

#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------------*
 | constructor                                                         fang 09/17 |
 *--------------------------------------------------------------------------------*/
STI::Partitioned::Partitioned(const Epetra_Comm& comm,  //! communicator
    const Teuchos::ParameterList& stidyn,  //! parameter list for scatra-thermo interaction
    const Teuchos::ParameterList&
        scatradyn,  //! scalar transport parameter list for scatra and thermo fields
    const Teuchos::ParameterList& solverparams_scatra,  //! solver parameter list for scatra field
    const Teuchos::ParameterList& solverparams_thermo   //! solver parameter list for thermo field
    )
    :  // instantiate base class
      Algorithm(comm, stidyn, scatradyn, solverparams_scatra, solverparams_thermo),
      couplingtype_(Teuchos::getIntegralValue<Inpar::STI::CouplingType>(stidyn, "COUPLINGTYPE")),
      omegamax_(stidyn.sublist("PARTITIONED").get<double>("OMEGAMAX"))
{
  // set control parameters for outer coupling iteration loop
  itermax_ = fieldparameters_->sublist("NONLINEAR").get<int>("ITEMAX_OUTER");
  itertol_ = fieldparameters_->sublist("NONLINEAR").get<double>("CONVTOL_OUTER");

  // initialize vectors for outer coupling iteration loop
  switch (couplingtype_)
  {
    case Inpar::STI::CouplingType::oneway_scatratothermo:
    case Inpar::STI::CouplingType::oneway_thermotoscatra:
      // do nothing
      break;

    case Inpar::STI::CouplingType::twoway_scatratothermo:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken_dofsplit:
    case Inpar::STI::CouplingType::twoway_thermotoscatra:
    case Inpar::STI::CouplingType::twoway_thermotoscatra_aitken:
    {
      // initialize increment vectors
      sca_tra_field()->phinp_inc() =
          Core::LinAlg::CreateVector(*sca_tra_field()->discretization()->dof_row_map(), true);
      thermo_field()->phinp_inc() =
          Core::LinAlg::CreateVector(*thermo_field()->discretization()->dof_row_map(), true);

      // initialize old increment vectors
      if (couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo_aitken or
          couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo_aitken_dofsplit)
        sca_tra_field()->phinp_inc_old() =
            Core::LinAlg::CreateVector(*sca_tra_field()->discretization()->dof_row_map(), true);
      else if (couplingtype_ == Inpar::STI::CouplingType::twoway_thermotoscatra_aitken)
        thermo_field()->phinp_inc_old() =
            Core::LinAlg::CreateVector(*thermo_field()->discretization()->dof_row_map(), true);

      // initialize relaxation parameter
      if (couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo)
        sca_tra_field()->omega().resize(1, stidyn.sublist("PARTITIONED").get<double>("OMEGA"));
      else if (couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo_aitken)
        sca_tra_field()->omega().resize(1, 1.);
      else if (couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo_aitken_dofsplit)
        sca_tra_field()->omega().resize(sca_tra_field()->num_dof_per_node(), 1.);
      else if (couplingtype_ == Inpar::STI::CouplingType::twoway_thermotoscatra)
        thermo_field()->omega().resize(1, stidyn.sublist("PARTITIONED").get<double>("OMEGA"));
      else if (couplingtype_ == Inpar::STI::CouplingType::twoway_thermotoscatra_aitken)
        thermo_field()->omega().resize(1, 1.);
      else
        thermo_field()->omega().resize(sca_tra_field()->num_dof_per_node(), 1.);

      break;
    }

    default:
    {
      FOUR_C_THROW("What the hell?!");
      break;
    }
  }

  return;
}  // STI::Partitioned::Partitioned


/*--------------------------------------------------------------------*
 | convergence check for outer iteration loop              fang 09/17 |
 *--------------------------------------------------------------------*/
bool STI::Partitioned::exit_outer_coupling() const
{
  // extract processor ID
  const int mypid = get_comm().MyPID();

  // compute vector norms
  double L2_scatra(0.), L2_scatra_inc(0.), L2_thermo(0.), L2_thermo_inc(0.);
  sca_tra_field()->phinp()->Norm2(&L2_scatra);
  sca_tra_field()->phinp_inc()->Norm2(&L2_scatra_inc);
  thermo_field()->phinp()->Norm2(&L2_thermo);
  thermo_field()->phinp_inc()->Norm2(&L2_thermo_inc);
  if (L2_scatra < 1.e-10) L2_scatra = 1.;
  if (L2_thermo < 1.e-10) L2_thermo = 1.;

  // print convergence status
  if (mypid == 0)
  {
    std::cout << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
    std::cout << "|                       OUTER ITERATION                        |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
    std::cout << "|- step/max -|- tol      [norm] -|- scatra-inc -|- thermo-inc -|" << std::endl;
    std::cout << "|  " << std::setw(3) << iter_ << "/" << std::setw(3) << itermax_ << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
              << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
              << L2_scatra_inc / L2_scatra << "   | " << std::setw(10) << std::setprecision(3)
              << std::scientific << L2_thermo_inc / L2_thermo << "   |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
  }

  // convergence check
  if (L2_scatra_inc / L2_scatra <= itertol_ and L2_thermo_inc / L2_thermo <= itertol_)
  {
    if (mypid == 0)
    {
      std::cout << "|   OUTER ITERATION LOOP CONVERGED AFTER ITERATION " << std::setw(3) << iter_
                << "/" << std::setw(3) << itermax_ << " !   |" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
    }

    return true;
  }

  // throw error in case maximum number of iteration steps is reached without convergence
  else if (iter_ == itermax_)
  {
    if (mypid == 0)
    {
      std::cout << "| >>>> not converged within maximum number of iteration steps! |" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
    }

    FOUR_C_THROW("Outer iteration did not converge within maximum number of iteration steps!");

    return true;
  }

  // proceed with next outer iteration step
  return false;
}  // STI::Partitioned::exit_outer_coupling()


/*--------------------------------------------------------------------------------*
 | evaluate time step using outer coupling iteration                   fang 09/17 |
 *--------------------------------------------------------------------------------*/
void STI::Partitioned::solve()
{
  switch (couplingtype_)
  {
    case Inpar::STI::CouplingType::oneway_scatratothermo:
    case Inpar::STI::CouplingType::oneway_thermotoscatra:
    {
      solve_one_way();
      break;
    }

    case Inpar::STI::CouplingType::twoway_scatratothermo:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken_dofsplit:
    case Inpar::STI::CouplingType::twoway_thermotoscatra:
    case Inpar::STI::CouplingType::twoway_thermotoscatra_aitken:
    {
      solve_two_way();
      break;
    }

    default:
    {
      FOUR_C_THROW("Dude... what's wrong with you?!");
      break;
    }
  }

  return;
}  // STI::Partitioned::Solve()


/*--------------------------------------------------------------------------------*
 | evaluate time step using one-way coupling iteration                 fang 09/17 |
 *--------------------------------------------------------------------------------*/
void STI::Partitioned::solve_one_way() const
{
  switch (couplingtype_)
  {
    case Inpar::STI::CouplingType::oneway_scatratothermo:
    {
      // pass thermo degrees of freedom to scatra discretization
      transfer_thermo_to_scatra(thermo_field()->phiafnp());

      // solve scatra field
      sca_tra_field()->solve();

      // pass scatra degrees of freedom to thermo discretization
      transfer_scatra_to_thermo(sca_tra_field()->phiafnp());

      // solve thermo field
      thermo_field()->solve();

      break;
    }

    case Inpar::STI::CouplingType::oneway_thermotoscatra:
    {
      // pass scatra degrees of freedom to thermo discretization
      transfer_scatra_to_thermo(sca_tra_field()->phiafnp());

      // solve thermo field
      thermo_field()->solve();

      // pass thermo degrees of freedom to scatra discretization
      transfer_thermo_to_scatra(thermo_field()->phiafnp());

      // solve scatra field
      sca_tra_field()->solve();

      break;
    }

    default:
    {
      FOUR_C_THROW("No, no, noooooooo...!");
      break;
    }
  }
}  // STI::Partitioned::solve_one_way()


/*----------------------------------------------------------------------*
 | evaluate time step using two-way coupling iteration       fang 09/17 |
 *----------------------------------------------------------------------*/
void STI::Partitioned::solve_two_way()
{
  // reset number of outer iterations
  iter_ = 0;

  switch (couplingtype_)
  {
    case Inpar::STI::CouplingType::twoway_scatratothermo:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken:
    case Inpar::STI::CouplingType::twoway_scatratothermo_aitken_dofsplit:
    {
      // initialize relaxed scatra state vector
      const Teuchos::RCP<Epetra_Vector> scatra_relaxed =
          Teuchos::rcp(new Epetra_Vector(*sca_tra_field()->phiafnp()));

      // begin outer iteration loop
      while (true)
      {
        // increment iteration number
        iter_++;

        // pass relaxed scatra degrees of freedom to thermo discretization
        transfer_scatra_to_thermo(scatra_relaxed);

        // store current thermo state vector
        if (thermo_field()->phinp_inc()->Update(1., *thermo_field()->phiafnp(), 0.))
          FOUR_C_THROW("Update failed!");

        // solve thermo field
        thermo_field()->solve();

        // compute increment of thermo state vector
        if (thermo_field()->phinp_inc()->Update(1., *thermo_field()->phiafnp(), -1.))
          FOUR_C_THROW("Update failed!");

        // pass thermo degrees of freedom to scatra discretization
        transfer_thermo_to_scatra(thermo_field()->phiafnp());

        // store current scatra state vector
        if (sca_tra_field()->phinp_inc()->Update(1., *scatra_relaxed, 0.))
          FOUR_C_THROW("Update failed!");

        // solve scatra field
        sca_tra_field()->solve();

        // compute increment of scatra state vector
        if (sca_tra_field()->phinp_inc()->Update(1., *sca_tra_field()->phiafnp(), -1.))
          FOUR_C_THROW("Update failed!");

        // convergence check
        if (exit_outer_coupling()) break;

        // perform static relaxation
        if (couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo)
        {
          if (scatra_relaxed->Update(
                  sca_tra_field()->omega()[0], *sca_tra_field()->phinp_inc(), 1.))
            FOUR_C_THROW("Update failed!");
        }

        // perform dynamic relaxation
        else
        {
          // compute difference between current and previous increments of scatra state vector
          Epetra_Vector scatra_inc_diff(*sca_tra_field()->phinp_inc());
          if (scatra_inc_diff.Update(-1., *sca_tra_field()->phinp_inc_old(), 1.))
            FOUR_C_THROW("Update failed!");

          if (couplingtype_ == Inpar::STI::CouplingType::twoway_scatratothermo_aitken)
          {
            // compute L2 norm of difference between current and previous increments of scatra state
            // vector
            double scatra_inc_diff_L2(0.);
            scatra_inc_diff.Norm2(&scatra_inc_diff_L2);

            // compute dot product between increment of scatra state vector and difference between
            // current and previous increments of scatra state vector
            double scatra_inc_dot_scatra_inc_diff(0.);
            if (scatra_inc_diff.Dot(*sca_tra_field()->phinp_inc(), &scatra_inc_dot_scatra_inc_diff))
              FOUR_C_THROW("Couldn't compute dot product!");

            // compute Aitken relaxation factor
            if (iter_ > 1 and scatra_inc_diff_L2 > 1.e-12)
              sca_tra_field()->omega()[0] *=
                  1 - scatra_inc_dot_scatra_inc_diff / (scatra_inc_diff_L2 * scatra_inc_diff_L2);

            // restrict Aitken relaxation factor if necessary
            if (omegamax_ > 0. and sca_tra_field()->omega()[0] > omegamax_)
              sca_tra_field()->omega()[0] = omegamax_;

            // perform Aitken relaxation
            if (scatra_relaxed->Update(
                    sca_tra_field()->omega()[0], *sca_tra_field()->phinp_inc(), 1.))
              FOUR_C_THROW("Update failed!");
          }

          else
          {
            // safety check
            if (sca_tra_field()->splitter() == Teuchos::null)
              FOUR_C_THROW("Map extractor was not initialized!");

            // loop over all degrees of freedom
            for (int idof = 0; idof < sca_tra_field()->splitter()->num_maps(); ++idof)
            {
              // extract subvectors associated with current degree of freedom
              const Teuchos::RCP<const Epetra_Vector> scatra_inc_dof =
                  sca_tra_field()->splitter()->extract_vector(*sca_tra_field()->phinp_inc(), idof);
              const Teuchos::RCP<const Epetra_Vector> scatra_inc_diff_dof =
                  sca_tra_field()->splitter()->extract_vector(scatra_inc_diff, idof);

              // compute L2 norm of difference between current and previous increments of current
              // degree of freedom
              double scatra_inc_diff_L2(0.);
              scatra_inc_diff_dof->Norm2(&scatra_inc_diff_L2);

              // compute dot product between increment of current degree of freedom and difference
              // between current and previous increments of current degree of freedom
              double scatra_inc_dot_scatra_inc_diff(0.);
              if (scatra_inc_diff_dof->Dot(*scatra_inc_dof, &scatra_inc_dot_scatra_inc_diff))
                FOUR_C_THROW("Couldn't compute dot product!");

              // compute Aitken relaxation factor for current degree of freedom
              if (iter_ > 1 and scatra_inc_diff_L2 > 1.e-12)
                sca_tra_field()->omega()[idof] *=
                    1 - scatra_inc_dot_scatra_inc_diff / (scatra_inc_diff_L2 * scatra_inc_diff_L2);

              // restrict Aitken relaxation factor if necessary
              if (omegamax_ > 0. and sca_tra_field()->omega()[idof] > omegamax_)
                sca_tra_field()->omega()[idof] = omegamax_;

              // perform Aitken relaxation for current degree of freedom
              sca_tra_field()->splitter()->add_vector(
                  *scatra_inc_dof, idof, *scatra_relaxed, sca_tra_field()->omega()[idof]);
            }
          }

          // update increment of scatra state vector
          if (sca_tra_field()->phinp_inc_old()->Update(1., *sca_tra_field()->phinp_inc(), 0.))
            FOUR_C_THROW("Update failed!");
        }
      }

      break;
    }

    case Inpar::STI::CouplingType::twoway_thermotoscatra:
    case Inpar::STI::CouplingType::twoway_thermotoscatra_aitken:
    {
      // initialize relaxed thermo state vector
      const Teuchos::RCP<Epetra_Vector> thermo_relaxed =
          Teuchos::rcp(new Epetra_Vector(*thermo_field()->phiafnp()));

      // begin outer iteration loop
      while (true)
      {
        // increment iteration number
        iter_++;

        // pass relaxed thermo degrees of freedom to scatra discretization
        transfer_thermo_to_scatra(thermo_relaxed);

        // store current scatra state vector
        if (sca_tra_field()->phinp_inc()->Update(1., *sca_tra_field()->phiafnp(), 0.))
          FOUR_C_THROW("Update failed!");

        // solve scatra field
        sca_tra_field()->solve();

        // compute increment of scatra state vector
        if (sca_tra_field()->phinp_inc()->Update(1., *sca_tra_field()->phiafnp(), -1.))
          FOUR_C_THROW("Update failed!");

        // pass scatra degrees of freedom to thermo discretization
        transfer_scatra_to_thermo(sca_tra_field()->phiafnp());

        // store current thermo state vector
        if (thermo_field()->phinp_inc()->Update(1., *thermo_relaxed, 0.))
          FOUR_C_THROW("Update failed!");

        // solve thermo field
        thermo_field()->solve();

        // compute increment of thermo state vector
        if (thermo_field()->phinp_inc()->Update(1., *thermo_field()->phiafnp(), -1.))
          FOUR_C_THROW("Update failed!");

        // convergence check
        if (exit_outer_coupling()) break;

        if (couplingtype_ == Inpar::STI::CouplingType::twoway_thermotoscatra_aitken)
        {
          // compute difference between current and previous increments of thermo state vector
          Epetra_Vector thermo_inc_diff(*thermo_field()->phinp_inc());
          if (thermo_inc_diff.Update(-1., *thermo_field()->phinp_inc_old(), 1.))
            FOUR_C_THROW("Update failed!");

          // compute L2 norm of difference between current and previous increments of thermo state
          // vector
          double thermo_inc_diff_L2(0.);
          thermo_inc_diff.Norm2(&thermo_inc_diff_L2);

          // compute dot product between increment of thermo state vector and difference between
          // current and previous increments of thermo state vector
          double thermo_inc_dot_thermo_inc_diff(0.);
          if (thermo_inc_diff.Dot(*thermo_field()->phinp_inc(), &thermo_inc_dot_thermo_inc_diff))
            FOUR_C_THROW("Couldn't compute dot product!");

          // compute Aitken relaxation factor
          if (iter_ > 1 and thermo_inc_diff_L2 > 1.e-12)
            thermo_field()->omega()[0] *=
                1 - thermo_inc_dot_thermo_inc_diff / (thermo_inc_diff_L2 * thermo_inc_diff_L2);

          // restrict Aitken relaxation factor if necessary
          if (omegamax_ > 0. and thermo_field()->omega()[0] > omegamax_)
            thermo_field()->omega()[0] = omegamax_;

          // update increment of thermo state vector
          if (thermo_field()->phinp_inc_old()->Update(1., *thermo_field()->phinp_inc(), 0.))
            FOUR_C_THROW("Update failed!");
        }

        // perform relaxation
        if (thermo_relaxed->Update(thermo_field()->omega()[0], *thermo_field()->phinp_inc(), 1.))
          FOUR_C_THROW("Update failed!");
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Please stop coding a bunch of crap...");
      break;
    }
  }

  return;
}  // STI::Partitioned::solve_two_way()

FOUR_C_NAMESPACE_CLOSE
