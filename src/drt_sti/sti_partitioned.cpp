/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned coupling algorithm for scatra-thermo interaction

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
*/
/*----------------------------------------------------------------------*/
#include "sti_partitioned.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

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
      couplingtype_(Teuchos::getIntegralValue<INPAR::STI::CouplingType>(stidyn, "COUPLINGTYPE")),
      omegamax_(stidyn.sublist("PARTITIONED").get<double>("OMEGAMAX"))
{
  // set control parameters for outer coupling iteration loop
  itermax_ = fieldparameters_->sublist("NONLINEAR").get<int>("ITEMAX_OUTER");
  itertol_ = fieldparameters_->sublist("NONLINEAR").get<double>("CONVTOL_OUTER");

  // initialize vectors for outer coupling iteration loop
  switch (couplingtype_)
  {
    case INPAR::STI::CouplingType::coupling_oneway_scatratothermo:
    case INPAR::STI::CouplingType::coupling_oneway_thermotoscatra:
      // do nothing
      break;

    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo:
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken:
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken_dofsplit:
    case INPAR::STI::CouplingType::coupling_twoway_thermotoscatra:
    case INPAR::STI::CouplingType::coupling_twoway_thermotoscatra_aitken:
    {
      // initialize increment vectors
      scatra_->ScaTraField()->PhinpInc() =
          LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(), true);
      thermo_->ScaTraField()->PhinpInc() =
          LINALG::CreateVector(*thermo_->ScaTraField()->Discretization()->DofRowMap(), true);

      // initialize old increment vectors
      if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken or
          couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken_dofsplit)
        scatra_->ScaTraField()->PhinpIncOld() =
            LINALG::CreateVector(*scatra_->ScaTraField()->Discretization()->DofRowMap(), true);
      else if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_thermotoscatra_aitken)
        thermo_->ScaTraField()->PhinpIncOld() =
            LINALG::CreateVector(*thermo_->ScaTraField()->Discretization()->DofRowMap(), true);

      // initialize relaxation parameter
      if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_scatratothermo)
        scatra_->ScaTraField()->Omega().resize(
            1, stidyn.sublist("PARTITIONED").get<double>("OMEGA"));
      else if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken)
        scatra_->ScaTraField()->Omega().resize(1, 1.);
      else if (couplingtype_ ==
               INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken_dofsplit)
        scatra_->ScaTraField()->Omega().resize(scatra_->ScaTraField()->NumDofPerNode(), 1.);
      else if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_thermotoscatra)
        thermo_->ScaTraField()->Omega().resize(
            1, stidyn.sublist("PARTITIONED").get<double>("OMEGA"));
      else if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_thermotoscatra_aitken)
        thermo_->ScaTraField()->Omega().resize(1, 1.);
      else
        thermo_->ScaTraField()->Omega().resize(scatra_->ScaTraField()->NumDofPerNode(), 1.);

      break;
    }

    default:
    {
      dserror("What the hell?!");
      break;
    }
  }

  return;
}  // STI::Partitioned::Partitioned


/*--------------------------------------------------------------------*
 | convergence check for outer iteration loop              fang 09/17 |
 *--------------------------------------------------------------------*/
bool STI::Partitioned::ExitOuterCoupling() const
{
  // extract processor ID
  const int mypid = Comm().MyPID();

  // compute vector norms
  double L2_scatra(0.), L2_scatra_inc(0.), L2_thermo(0.), L2_thermo_inc(0.);
  scatra_->ScaTraField()->Phinp()->Norm2(&L2_scatra);
  scatra_->ScaTraField()->PhinpInc()->Norm2(&L2_scatra_inc);
  thermo_->ScaTraField()->Phinp()->Norm2(&L2_thermo);
  thermo_->ScaTraField()->PhinpInc()->Norm2(&L2_thermo_inc);
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

    dserror("Outer iteration did not converge within maximum number of iteration steps!");

    return true;
  }

  // proceed with next outer iteration step
  return false;
}  // STI::Partitioned::ExitOuterCoupling()


/*--------------------------------------------------------------------------------*
 | evaluate time step using outer coupling iteration                   fang 09/17 |
 *--------------------------------------------------------------------------------*/
void STI::Partitioned::Solve()
{
  switch (couplingtype_)
  {
    case INPAR::STI::CouplingType::coupling_oneway_scatratothermo:
    case INPAR::STI::CouplingType::coupling_oneway_thermotoscatra:
    {
      SolveOneWay();
      break;
    }

    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo:
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken:
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken_dofsplit:
    case INPAR::STI::CouplingType::coupling_twoway_thermotoscatra:
    case INPAR::STI::CouplingType::coupling_twoway_thermotoscatra_aitken:
    {
      SolveTwoWay();
      break;
    }

    default:
    {
      dserror("Dude... what's wrong with you?!");
      break;
    }
  }

  return;
}  // STI::Partitioned::Solve()


/*--------------------------------------------------------------------------------*
 | evaluate time step using one-way coupling iteration                 fang 09/17 |
 *--------------------------------------------------------------------------------*/
void STI::Partitioned::SolveOneWay() const
{
  switch (couplingtype_)
  {
    case INPAR::STI::CouplingType::coupling_oneway_scatratothermo:
    {
      // pass thermo degrees of freedom to scatra discretization
      TransferThermoToScatra(thermo_->ScaTraField()->Phiafnp());

      // solve scatra field
      scatra_->ScaTraField()->Solve();

      // pass scatra degrees of freedom to thermo discretization
      TransferScatraToThermo(scatra_->ScaTraField()->Phiafnp());

      // solve thermo field
      thermo_->ScaTraField()->Solve();

      break;
    }

    case INPAR::STI::CouplingType::coupling_oneway_thermotoscatra:
    {
      // pass scatra degrees of freedom to thermo discretization
      TransferScatraToThermo(scatra_->ScaTraField()->Phiafnp());

      // solve thermo field
      thermo_->ScaTraField()->Solve();

      // pass thermo degrees of freedom to scatra discretization
      TransferThermoToScatra(thermo_->ScaTraField()->Phiafnp());

      // solve scatra field
      scatra_->ScaTraField()->Solve();

      break;
    }

    default:
    {
      dserror("No, no, noooooooo...!");
      break;
    }
  }
}  // STI::Partitioned::SolveOneWay()


/*----------------------------------------------------------------------*
 | evaluate time step using two-way coupling iteration       fang 09/17 |
 *----------------------------------------------------------------------*/
void STI::Partitioned::SolveTwoWay()
{
  // reset number of outer iterations
  iter_ = 0;

  switch (couplingtype_)
  {
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo:
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken:
    case INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken_dofsplit:
    {
      // initialize relaxed scatra state vector
      const Teuchos::RCP<Epetra_Vector> scatra_relaxed =
          Teuchos::rcp(new Epetra_Vector(*scatra_->ScaTraField()->Phiafnp()));

      // begin outer iteration loop
      while (true)
      {
        // increment iteration number
        iter_++;

        // pass relaxed scatra degrees of freedom to thermo discretization
        TransferScatraToThermo(scatra_relaxed);

        // store current thermo state vector
        if (thermo_->ScaTraField()->PhinpInc()->Update(1., *thermo_->ScaTraField()->Phiafnp(), 0.))
          dserror("Update failed!");
        ;

        // solve thermo field
        thermo_->ScaTraField()->Solve();

        // compute increment of thermo state vector
        if (thermo_->ScaTraField()->PhinpInc()->Update(1., *thermo_->ScaTraField()->Phiafnp(), -1.))
          dserror("Update failed!");

        // pass thermo degrees of freedom to scatra discretization
        TransferThermoToScatra(thermo_->ScaTraField()->Phiafnp());

        // store current scatra state vector
        if (scatra_->ScaTraField()->PhinpInc()->Update(1., *scatra_relaxed, 0.))
          dserror("Update failed!");

        // solve scatra field
        scatra_->ScaTraField()->Solve();

        // compute increment of scatra state vector
        if (scatra_->ScaTraField()->PhinpInc()->Update(1., *scatra_->ScaTraField()->Phiafnp(), -1.))
          dserror("Update failed!");

        // convergence check
        if (ExitOuterCoupling()) break;

        // perform static relaxation
        if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_scatratothermo)
        {
          if (scatra_relaxed->Update(
                  scatra_->ScaTraField()->Omega()[0], *scatra_->ScaTraField()->PhinpInc(), 1.))
            dserror("Update failed!");
        }

        // perform dynamic relaxation
        else
        {
          // compute difference between current and previous increments of scatra state vector
          Epetra_Vector scatra_inc_diff(*scatra_->ScaTraField()->PhinpInc());
          if (scatra_inc_diff.Update(-1., *scatra_->ScaTraField()->PhinpIncOld(), 1.))
            dserror("Update failed!");

          if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_scatratothermo_aitken)
          {
            // compute L2 norm of difference between current and previous increments of scatra state
            // vector
            double scatra_inc_diff_L2(0.);
            scatra_inc_diff.Norm2(&scatra_inc_diff_L2);

            // compute dot product between increment of scatra state vector and difference between
            // current and previous increments of scatra state vector
            double scatra_inc_dot_scatra_inc_diff(0.);
            if (scatra_inc_diff.Dot(
                    *scatra_->ScaTraField()->PhinpInc(), &scatra_inc_dot_scatra_inc_diff))
              dserror("Couldn't compute dot product!");

            // compute Aitken relaxation factor
            if (iter_ > 1 and scatra_inc_diff_L2 > 1.e-12)
              scatra_->ScaTraField()->Omega()[0] *=
                  1 - scatra_inc_dot_scatra_inc_diff / (scatra_inc_diff_L2 * scatra_inc_diff_L2);

            // restrict Aitken relaxation factor if necessary
            if (omegamax_ > 0. and scatra_->ScaTraField()->Omega()[0] > omegamax_)
              scatra_->ScaTraField()->Omega()[0] = omegamax_;

            // perform Aitken relaxation
            if (scatra_relaxed->Update(
                    scatra_->ScaTraField()->Omega()[0], *scatra_->ScaTraField()->PhinpInc(), 1.))
              dserror("Update failed!");
          }

          else
          {
            // safety check
            if (scatra_->ScaTraField()->Splitter() == Teuchos::null)
              dserror("Map extractor was not initialized!");

            // loop over all degrees of freedom
            for (int idof = 0; idof < scatra_->ScaTraField()->Splitter()->NumMaps(); ++idof)
            {
              // extract subvectors associated with current degree of freedom
              const Teuchos::RCP<const Epetra_Vector> scatra_inc_dof =
                  scatra_->ScaTraField()->Splitter()->ExtractVector(
                      *scatra_->ScaTraField()->PhinpInc(), idof);
              const Teuchos::RCP<const Epetra_Vector> scatra_inc_diff_dof =
                  scatra_->ScaTraField()->Splitter()->ExtractVector(scatra_inc_diff, idof);

              // compute L2 norm of difference between current and previous increments of current
              // degree of freedom
              double scatra_inc_diff_L2(0.);
              scatra_inc_diff_dof->Norm2(&scatra_inc_diff_L2);

              // compute dot product between increment of current degree of freedom and difference
              // between current and previous increments of current degree of freedom
              double scatra_inc_dot_scatra_inc_diff(0.);
              if (scatra_inc_diff_dof->Dot(*scatra_inc_dof, &scatra_inc_dot_scatra_inc_diff))
                dserror("Couldn't compute dot product!");

              // compute Aitken relaxation factor for current degree of freedom
              if (iter_ > 1 and scatra_inc_diff_L2 > 1.e-12)
                scatra_->ScaTraField()->Omega()[idof] *=
                    1 - scatra_inc_dot_scatra_inc_diff / (scatra_inc_diff_L2 * scatra_inc_diff_L2);

              // restrict Aitken relaxation factor if necessary
              if (omegamax_ > 0. and scatra_->ScaTraField()->Omega()[idof] > omegamax_)
                scatra_->ScaTraField()->Omega()[idof] = omegamax_;

              // perform Aitken relaxation for current degree of freedom
              scatra_->ScaTraField()->Splitter()->AddVector(
                  *scatra_inc_dof, idof, *scatra_relaxed, scatra_->ScaTraField()->Omega()[idof]);
            }
          }

          // update increment of scatra state vector
          if (scatra_->ScaTraField()->PhinpIncOld()->Update(
                  1., *scatra_->ScaTraField()->PhinpInc(), 0.))
            dserror("Update failed!");
        }
      }

      break;
    }

    case INPAR::STI::CouplingType::coupling_twoway_thermotoscatra:
    case INPAR::STI::CouplingType::coupling_twoway_thermotoscatra_aitken:
    {
      // initialize relaxed thermo state vector
      const Teuchos::RCP<Epetra_Vector> thermo_relaxed =
          Teuchos::rcp(new Epetra_Vector(*thermo_->ScaTraField()->Phiafnp()));

      // begin outer iteration loop
      while (true)
      {
        // increment iteration number
        iter_++;

        // pass relaxed thermo degrees of freedom to scatra discretization
        TransferThermoToScatra(thermo_relaxed);

        // store current scatra state vector
        if (scatra_->ScaTraField()->PhinpInc()->Update(1., *scatra_->ScaTraField()->Phiafnp(), 0.))
          dserror("Update failed!");

        // solve scatra field
        scatra_->ScaTraField()->Solve();

        // compute increment of scatra state vector
        if (scatra_->ScaTraField()->PhinpInc()->Update(1., *scatra_->ScaTraField()->Phiafnp(), -1.))
          dserror("Update failed!");

        // pass scatra degrees of freedom to thermo discretization
        TransferScatraToThermo(scatra_->ScaTraField()->Phiafnp());

        // store current thermo state vector
        if (thermo_->ScaTraField()->PhinpInc()->Update(1., *thermo_relaxed, 0.))
          dserror("Update failed!");

        // solve thermo field
        thermo_->ScaTraField()->Solve();

        // compute increment of thermo state vector
        if (thermo_->ScaTraField()->PhinpInc()->Update(1., *thermo_->ScaTraField()->Phiafnp(), -1.))
          dserror("Update failed!");

        // convergence check
        if (ExitOuterCoupling()) break;

        if (couplingtype_ == INPAR::STI::CouplingType::coupling_twoway_thermotoscatra_aitken)
        {
          // compute difference between current and previous increments of thermo state vector
          Epetra_Vector thermo_inc_diff(*thermo_->ScaTraField()->PhinpInc());
          if (thermo_inc_diff.Update(-1., *thermo_->ScaTraField()->PhinpIncOld(), 1.))
            dserror("Update failed!");

          // compute L2 norm of difference between current and previous increments of thermo state
          // vector
          double thermo_inc_diff_L2(0.);
          thermo_inc_diff.Norm2(&thermo_inc_diff_L2);

          // compute dot product between increment of thermo state vector and difference between
          // current and previous increments of thermo state vector
          double thermo_inc_dot_thermo_inc_diff(0.);
          if (thermo_inc_diff.Dot(
                  *thermo_->ScaTraField()->PhinpInc(), &thermo_inc_dot_thermo_inc_diff))
            dserror("Couldn't compute dot product!");

          // compute Aitken relaxation factor
          if (iter_ > 1 and thermo_inc_diff_L2 > 1.e-12)
            thermo_->ScaTraField()->Omega()[0] *=
                1 - thermo_inc_dot_thermo_inc_diff / (thermo_inc_diff_L2 * thermo_inc_diff_L2);

          // restrict Aitken relaxation factor if necessary
          if (omegamax_ > 0. and thermo_->ScaTraField()->Omega()[0] > omegamax_)
            thermo_->ScaTraField()->Omega()[0] = omegamax_;

          // update increment of thermo state vector
          if (thermo_->ScaTraField()->PhinpIncOld()->Update(
                  1., *thermo_->ScaTraField()->PhinpInc(), 0.))
            dserror("Update failed!");
        }

        // perform relaxation
        if (thermo_relaxed->Update(
                thermo_->ScaTraField()->Omega()[0], *thermo_->ScaTraField()->PhinpInc(), 1.))
          dserror("Update failed!");
      }

      break;
    }

    default:
    {
      dserror("Please stop coding a bunch of crap...");
      break;
    }
  }

  return;
}  // STI::Partitioned::SolveTwoWay()
