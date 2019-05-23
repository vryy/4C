/*----------------------------------------------------------------------*/
/*!
\file elch_algorithm.cpp

\brief Basis of all ELCH algorithms

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251

\level 2
*/
/*----------------------------------------------------------------------*/

#include "../linalg/linalg_mapextractor.H"
#include "../drt_scatra/scatra_timint_elch.H"
#include "elch_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& elchcontrol,
    const Teuchos::ParameterList& scatradyn, const Teuchos::ParameterList& fdyn,
    const Teuchos::ParameterList& solverparams)
    : ScaTraAlgorithm(comm, scatradyn, fdyn, "scatra", solverparams)
{
  // Setup of TurbulenceStatisticManager is performed in the
  // constructor of class ScaTraFluidCouplingAlgorithm

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::~Algorithm() { return; }


/*----------------------------------------------------------------------*
 | Provide information about initial field                   fang 08/14 |
 *----------------------------------------------------------------------*/
void ELCH::Algorithm::PrepareTimeLoop()
{
  // provide information about initial field (do not do for restarts!)
  if (Step() == 0)
  {
    ScaTraField()->OutputProblemSpecific();
    ScaTraField()->OutputTotalAndMeanScalars();

    // compute error for problems with analytical solution (initial field!)
    ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
  }

  return;
}


/*----------------------------------------------------------------------*
 | Print scatra solver type to screen                        fang 08/14 |
 *----------------------------------------------------------------------*/
void ELCH::Algorithm::PrintScaTraSolver()
{
  if (Comm().MyPID() == 0)
    std::cout
        << "\n****************************\n      ELCH SOLVER\n****************************\n";

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool ELCH::Algorithm::ConvergenceCheck(
    int natconvitnum, const int natconvitmax, const double natconvittol)
{
  // convergence check based on the concentration, potential and
  // velocity increment

  //   | concentration increment |_2
  //  -------------------------------- < Tolerance
  //     | concentration_n+1 |_2

  bool stopnonliniter = false;
  Teuchos::RCP<LINALG::MapExtractor> conpotsplitter = ScaTraField()->Splitter();
  // Variables to save different L2 - Norms

  double potincnorm_L2(0.0);
  double potnorm_L2(0.0);
  double velincnorm_L2(0.0);
  double velnorm_L2(0.0);
  double connorm_L2(0.0);
  double conincnorm_L2(0.0);

  // Calculate velocity increment and velocity L2 - Norm
  // velincnp_ = 1.0 * convelnp_ - 1.0 * conveln_

  velincnp_->Update(1.0, *FluidField()->ExtractVelocityPart(FluidField()->Velnp()), -1.0);
  velincnp_->Norm2(&velincnorm_L2);  // Estimation of the L2 - norm save values to both variables
                                     // (velincnorm_L2 and velnorm_L2)
  FluidField()->ExtractVelocityPart(FluidField()->Velnp())->Norm2(&velnorm_L2);

  // Calculate concentration increment and concentration L2 - Norm
  phiincnp_->Update(1.0, *ScaTraField()->Phinp(), -1.0);
  Teuchos::RCP<Epetra_Vector> onlycon = conpotsplitter->ExtractOtherVector(phiincnp_);
  onlycon->Norm2(&conincnorm_L2);
  conpotsplitter->ExtractOtherVector(ScaTraField()->Phinp(), onlycon);
  onlycon->Norm2(&connorm_L2);

  // Calculate potential increment and potential L2 - Norm
  Teuchos::RCP<Epetra_Vector> onlypot = conpotsplitter->ExtractCondVector(phiincnp_);
  onlypot->Norm2(&potincnorm_L2);
  conpotsplitter->ExtractCondVector(ScaTraField()->Phinp(), onlypot);
  onlypot->Norm2(&potnorm_L2);

  // care for the case that there is (almost) zero temperature or velocity
  // (usually not required for temperature)
  if (velnorm_L2 < 1e-6) velnorm_L2 = 1.0;
  if (connorm_L2 < 1e-6) connorm_L2 = 1.0;
  if (potnorm_L2 < 1e-6) potnorm_L2 = 1.0;

  // Print the incremental based convergence check to the screen
  if (natconvitnum != 1)
  {
    if (Comm().MyPID() == 0)
    {
      std::cout << "\n";
      std::cout
          << "*****************************************************************************\n";
      std::cout << "                          OUTER ITERATION STEP\n";
      std::cout
          << "*****************************************************************************\n";
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc ---|-- pot-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   | %10.3E  | %10.3E  |", natconvitnum,
          natconvitmax, natconvittol, conincnorm_L2 / connorm_L2, potincnorm_L2 / potnorm_L2,
          velincnorm_L2 / velnorm_L2);
      printf("\n");
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
    }

    // Converged or Not
    if ((conincnorm_L2 / connorm_L2 <= natconvittol) &&
        (potincnorm_L2 / potnorm_L2 <= natconvittol) &&
        (velincnorm_L2 / velnorm_L2 <= natconvittol))
    // if ((incconnorm_L2/connorm_L2 <= natconvittol))
    {
      stopnonliniter = true;
      if (Comm().MyPID() == 0)
      {
        printf("| Outer Iteration loop converged after iteration %3d/%3d                    |\n",
            natconvitnum, natconvitmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
      }
    }
    else
    {
      if (Comm().MyPID() == 0)
      {
        printf("| Outer Iteration loop is not converged after iteration %3d/%3d             |\n",
            natconvitnum, natconvitmax);
        printf("+---------------------------------------------------------------------------+");
        printf("\n");
        printf("\n");
      }
    }
  }
  else
  {
    // first outer iteration loop: fluid solver has not got the new density yet
    // => minimum two outer iteration loops
    stopnonliniter = false;
    if (Comm().MyPID() == 0)
    {
      std::cout << "\n";
      std::cout
          << "*****************************************************************************\n";
      std::cout << "                          OUTER ITERATION STEP\n";
      std::cout
          << "*****************************************************************************\n";
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
      printf("|- step/max -|- tol      [norm] -|-- con-inc ---|-- pot-inc --|-- vel-inc --|\n");
      printf("|  %3d/%3d   | %10.3E[L_2 ]  |       -      |      -      |      -      |",
          natconvitnum, natconvitmax, natconvittol);
      printf("\n");
      printf("+------------+-------------------+--------------+-------------+-------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next timestep
  // itemax = 1 is also possible for segregated coupling approaches (not fully implicit)
  if (natconvitnum == natconvitmax)
  {
    if (((conincnorm_L2 / connorm_L2 > natconvittol) ||
            (potincnorm_L2 / potnorm_L2 > natconvittol) ||
            (velincnorm_L2 / velnorm_L2 > natconvittol)) or
        (natconvitmax == 1))
    {
      stopnonliniter = true;
      if ((Comm().MyPID() == 0))
      {
        printf("|     >>>>>> not converged in itemax steps!     |\n");
        printf("+-----------------------------------------------+\n");
        printf("\n");
        printf("\n");
      }
    }
  }

  return stopnonliniter;
}
