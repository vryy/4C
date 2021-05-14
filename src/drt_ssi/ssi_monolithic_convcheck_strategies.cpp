/*----------------------------------------------------------------------*/
/*! \file
\brief strategies for Newton-Raphson convergence check for monolithic scalar-structure interaction
problems

To keep the time integrator class for monolithic scalar-structure interaction problems as plain as
possible, the convergence check for the Newton-Raphson iteration has been encapsulated within
separate strategy classes. Every specific convergence check strategy (e.g., for monolithic
scalar-structure interaction problems involving standard scalar transport or electrochemistry)
computes, checks, and outputs different relevant vector norms and is implemented in a subclass
derived from an abstract, purely virtual interface class.

\level 2


 */
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::SSIMono::ConvCheckStrategyBase::ConvCheckStrategyBase(
    const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
    )
    : itermax_(parameters.get<int>("ITEMAX")),
      itertol_(parameters.sublist("MONOLITHIC").get<double>("CONVTOL")),
      restol_(parameters.sublist("MONOLITHIC").get<double>("ABSTOLRES"))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyBase::CheckL2Norm(
    double& incnorm, double& resnorm, double& dofnorm) const
{
  if (std::isnan(incnorm) or std::isnan(resnorm) or std::isnan(dofnorm))
    dserror("Vector norm is not a number!");
  if (std::isinf(incnorm) or std::isinf(resnorm) or std::isinf(dofnorm))
    dserror("Vector norm is infinity!");

  if (dofnorm < 1.e-10) dofnorm = 1.e-10;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyBase::GetAndCheckL2NormStructure(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.MapsSubProblems()
      ->ExtractVector(
          ssi_mono.ssi_vectors_->Increment(), ssi_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&incnorm);

  ssi_mono.MapsSubProblems()
      ->ExtractVector(
          ssi_mono.ssi_vectors_->Residual(), ssi_mono.GetProblemPosition(Subproblem::structure))
      ->Norm2(&resnorm);

  ssi_mono.StructureField()->Dispnp()->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyStd::GetAndCheckL2NormScaTra(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
          ssi_mono.GetProblemPosition(Subproblem::scalar_transport))
      ->Norm2(&incnorm);

  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
          ssi_mono.GetProblemPosition(Subproblem::scalar_transport))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraField()->Phinp()->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyStd::ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const
{
  // initialize exit flag
  bool exit(false);

  double scatraincnorm(0.0), scatraresnorm(0.0), scatradofnorm(0.0), structureincnorm(0.0),
      structureresnorm(0.0), structuredofnorm(0.0);

  GetAndCheckL2NormScaTra(ssi_mono, scatraincnorm, scatraresnorm, scatradofnorm);
  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  // first Newton-Raphson iteration
  if (ssi_mono.IterationCount() == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+"
                << std::endl;
      std::cout << "|- step/max -|- tolerance[norm] -|- scatra-res -|- scatra-inc -|- struct-res "
                   "-|- struct-inc -|"
                << std::endl;

      // print first line of convergence table to screen
      // solution increment not yet available during first Newton-Raphson iteration
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << scatraresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << structureresnorm
                << "   |      --      | "
                << "(       --      , te = " << std::setw(10) << std::setprecision(3)
                << ssi_mono.dtele_ << ")" << std::endl;
    }
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << scatraresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << scatraincnorm / scatradofnorm
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << structureresnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << structureincnorm / structuredofnorm
                << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
                << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
                << std::endl;
    }

    // convergence check
    if (scatraresnorm <= itertol_ and structureresnorm <= itertol_ and
        scatraincnorm / scatradofnorm <= itertol_ and
        structureincnorm / structuredofnorm <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (scatraresnorm < restol_ and structureresnorm < restol_) exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence
  if (ssi_mono.IterationCount() == itermax_ and !exit)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+"
                << std::endl;
      std::cout << "|      Newton-Raphson method has not converged after a maximum number of "
                << std::setw(2) << itermax_ << " iterations!      |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if (exit and ssi_mono.Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+"
              << std::endl;
  }

  return exit;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElch::GetAndCheckL2NormConc(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
              ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
              ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.ScaTraField()->Phinp())
      ->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElch::GetAndCheckL2NormPot(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
              ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
              ssi_mono.GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.ScaTraField()->Phinp())
      ->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElch::ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const
{
  // initialize exit flag
  bool exit(false);

  double concincnorm(0.0), concresnorm(0.0), concdofnorm(0.0), potdofnorm(0.0), potincnorm(0.0),
      potresnorm(0.0), structuredofnorm(0.0), structureresnorm(0.0), structureincnorm(0.0);

  GetAndCheckL2NormConc(ssi_mono, concincnorm, concresnorm, concdofnorm);
  GetAndCheckL2NormPot(ssi_mono, potincnorm, potresnorm, potdofnorm);
  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  // first Newton-Raphson iteration
  if (ssi_mono.IterationCount() == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+"
                << std::endl;
      std::cout << "|- step/max -|- tolerance[norm] -|-- conc-res --|-- conc-inc --|-- pot-res "
                   "---|-- pot-inc ---|- struct-res -|- struct-inc -|"
                << std::endl;

      // print first line of convergence table to screen
      // solution increment not yet available during first Newton-Raphson iteration
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << potresnorm << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureresnorm
                << "   |      --      | "
                << "(       --      , te = " << std::setw(10) << std::setprecision(3)
                << ssi_mono.dtele_ << ")" << std::endl;
    }
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << concincnorm / concdofnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << potresnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific
                << potincnorm / potdofnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << structureresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << structureincnorm / structuredofnorm
                << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
                << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
                << std::endl;
    }

    // convergence check
    if (concresnorm <= itertol_ and potresnorm <= itertol_ and structureresnorm <= itertol_ and
        concincnorm / concdofnorm <= itertol_ and potincnorm / potdofnorm <= itertol_ and
        structureincnorm / structuredofnorm <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (concresnorm < restol_ and potresnorm < restol_ and structureresnorm < restol_) exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence

  if (ssi_mono.IterationCount() == itermax_ and !exit)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+"
                << std::endl;
      std::cout << "|                     Newton-Raphson method has not converged after a maximum "
                   "number of "
                << std::setw(2) << itermax_ << " iterations!                     |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if (exit and ssi_mono.Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+"
              << std::endl;
  }

  return exit;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElch::ExitNewtonRaphsonInitPotCalc(
    const SSI::SSIMono& ssi_mono, const int init_pot_iternum) const
{
  bool converged = false;

  double scatra_pot_dofnorm = 0.0, scatra_pot_resnorm = 0.0, scatra_pot_incnorm = 0.0;

  GetAndCheckL2NormPot(ssi_mono, scatra_pot_incnorm, scatra_pot_resnorm, scatra_pot_dofnorm);

  if (init_pot_iternum == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header
      std::cout << "Calculating initial field for electric potential" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- pot-res ---|-- pot-inc ---|" << std::endl;

      // print only norm of residuals
      std::cout << "|  " << std::setw(3) << init_pot_iternum << "/" << std::setw(3) << itermax_
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
                << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << scatra_pot_resnorm << "   |      --      |" << std::endl;
    }
    if (scatra_pot_resnorm < restol_)
    {
      if (ssi_mono.Comm().MyPID() == 0)
      {
        std::cout << "+------------+-------------------+--------------+--------------+"
                  << std::endl;
      }
      converged = true;
    }
  }
  else
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << init_pot_iternum << "/" << std::setw(3) << itermax_
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
                << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << scatra_pot_resnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << scatra_pot_incnorm / scatra_pot_dofnorm << "   |"
                << std::endl;
    }

    // convergence check
    if ((scatra_pot_resnorm <= itertol_ and scatra_pot_incnorm / scatra_pot_dofnorm <= itertol_) or
        scatra_pot_resnorm < restol_)
    {
      if (ssi_mono.Comm().MyPID() == 0)
      {
        std::cout << "+------------+-------------------+--------------+--------------+" << std::endl
                  << std::endl;
      }
      converged = true;
    }
  }

  // warn if maximum number of iterations is reached without convergence
  if (init_pot_iternum == itermax_)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+--------------------------------------------------------------+" << std::endl;
      std::cout << "|            >>>>>> not converged!                             |" << std::endl;
      std::cout << "+--------------------------------------------------------------+" << std::endl;
    }
    converged = true;
  }

  return converged;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::GetAndCheckL2NormScaTraManifoldConc(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.ssi_vectors_->Increment(), ssi_mono.GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.ssi_vectors_->Residual(), ssi_mono.GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractOtherVector(ssi_mono.ScaTraManifold()->Phinp())
      ->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::GetAndCheckL2NormScaTraManifoldPot(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.ssi_vectors_->Increment(), ssi_mono.GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.MapsSubProblems()->ExtractVector(
          ssi_mono.ssi_vectors_->Residual(), ssi_mono.GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.ScaTraManifold()->Phinp())
      ->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::ExitNewtonRaphson(
    const SSI::SSIMono& ssi_mono) const
{
  // initialize exit flag
  bool exit(false);

  double concincnorm(0.0), concresnorm(0.0), concdofnorm(0.0), manifoldpotdofnorm(0.0),
      manifoldpotincnorm(0.0), manifoldpotresnorm(0.0), manifoldconcdofnorm(0.0),
      manifoldconcincnorm(0.0), manifoldconcresnorm(0.0), potdofnorm(0.0), potincnorm(0.0),
      potresnorm(0.0), structuredofnorm(0.0), structureresnorm(0.0), structureincnorm(0.0);

  GetAndCheckL2NormConc(ssi_mono, concincnorm, concresnorm, concdofnorm);
  GetAndCheckL2NormPot(ssi_mono, potincnorm, potresnorm, potdofnorm);

  GetAndCheckL2NormScaTraManifoldPot(
      ssi_mono, manifoldpotincnorm, manifoldpotresnorm, manifoldpotdofnorm);
  GetAndCheckL2NormScaTraManifoldConc(
      ssi_mono, manifoldconcincnorm, manifoldconcresnorm, manifoldconcdofnorm);

  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  // first Newton-Raphson iteration
  if (ssi_mono.IterationCount() == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+--------------+--------------+----"
                   "----------+--------------+"
                << std::endl;
      std::cout
          << "+------------+-------------------+--                        scatra                  "
             "       --+--        structure        --+--                        manifold          "
             "             --+"
          << std::endl;
      std::cout
          << "|- step/max -|- tolerance[norm] -|-- conc-res --|-- conc-inc --|-- pot-res  --|-- "
             "pot-inc  --|--   res    --|--   inc    --|-- conc-res --|-- conc-inc --|-- "
             "pot-res  --|-- pot-inc  --|"
          << std::endl;

      // print first line of convergence table to screen
      // solution increment not yet available during first Newton-Raphson iteration
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << potresnorm << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << structureresnorm
                << "   |      --      | " << std::setw(10) << std::setprecision(3)
                << std::scientific << manifoldconcresnorm << "   |      --      |  "
                << std::scientific << manifoldpotresnorm
                << "   |      --      | (       --      , te = " << std::setw(10)
                << std::setprecision(3) << ssi_mono.dtele_ << ")" << std::endl;
    }
  }

  // subsequent Newton-Raphson iterations
  else
  {
    // print current line of convergence table to screen
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << concresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << concincnorm / concdofnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << potresnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific
                << potincnorm / potdofnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << structureresnorm << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific << structureincnorm / structuredofnorm
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << manifoldconcresnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << manifoldconcincnorm / manifoldconcdofnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << manifoldpotresnorm
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << manifoldpotincnorm / manifoldpotdofnorm << "   | (ts = " << std::setw(10)
                << std::setprecision(3) << ssi_mono.dtsolve_ << ", te = " << std::setw(10)
                << std::setprecision(3) << ssi_mono.dtele_ << ")" << std::endl;
    }

    // convergence check
    const bool scatra_conv = concresnorm <= itertol_ and potresnorm <= itertol_ and
                             concincnorm / concdofnorm <= itertol_ and
                             potincnorm / potdofnorm <= itertol_;
    const bool struct_conv =
        structureresnorm <= itertol_ and structureincnorm / structuredofnorm <= itertol_;

    const bool manifold_conv = manifoldconcresnorm <= itertol_ and
                               manifoldpotresnorm <= itertol_ and
                               manifoldpotincnorm / manifoldpotdofnorm <= itertol_ and
                               manifoldconcincnorm / manifoldconcdofnorm <= itertol_;

    // exit Newton-Raphson iteration upon convergence
    if (scatra_conv and struct_conv and manifold_conv) exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (concresnorm < restol_ and potresnorm < restol_ and structureresnorm < restol_ and
      manifoldpotresnorm < restol_ and manifoldconcresnorm < restol_)
    exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence

  if (ssi_mono.IterationCount() == itermax_ and !exit)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+--------------+--------------+----"
                   "----------+--------------+"
                << std::endl;
      std::cout
          << "|              >> Newton-Raphson method has not converged after a maximum "
             "number of "
          << itermax_
          << " iterations! <<                                                                      "
             "             |"
          << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if (exit and ssi_mono.Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+--------------+--------------+--------"
                 "------+--------------+"
              << std::endl;
  }

  return exit;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::ExitNewtonRaphsonInitPotCalc(
    const SSI::SSIMono& ssi_mono, const int init_pot_iternum) const
{
  bool converged = false;

  double scatra_pot_dofnorm = 0.0, manifold_pot_dofnorm = 0.0, scatra_pot_resnorm = 0.0,
         manifold_pot_resnorm = 0.0, scatra_pot_incnorm = 0.0, manifold_pot_incnorm = 0.0;

  GetAndCheckL2NormScaTraManifoldPot(
      ssi_mono, manifold_pot_incnorm, manifold_pot_resnorm, manifold_pot_dofnorm);
  GetAndCheckL2NormPot(ssi_mono, scatra_pot_incnorm, scatra_pot_resnorm, scatra_pot_dofnorm);

  if (init_pot_iternum == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header
      std::cout << "Calculating initial field for electric potential" << std::endl;
      std::cout
          << "+------------+-------------------+--------------+--------------+--------------+--"
             "------------+"
          << std::endl;
      std::cout << "+------------+-------------------+--          scatra         --|--             "
                   "manifold    --|"
                << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- pot-res  --|-- pot-inc  --|-- pot-res  "
                   "--|-- pot-inc  --|"
                << std::endl;

      // print only norm of residuals
      std::cout << "|  " << std::setw(3) << init_pot_iternum << "/" << std::setw(3) << itermax_
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
                << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << scatra_pot_resnorm << "   |      --      | " << std::setw(10)
                << std::setprecision(3) << std::scientific << manifold_pot_resnorm
                << "   |      --      | " << std::endl;
    }
    if (scatra_pot_resnorm < restol_ and manifold_pot_resnorm < restol_)
    {
      if (ssi_mono.Comm().MyPID() == 0)
      {
        std::cout << "+------------+-------------------+--------------+--------------+-----------"
                     "---+--------------+"
                  << std::endl;
      }
      converged = true;
    }
  }
  else
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "|  " << std::setw(3) << init_pot_iternum << "/" << std::setw(3) << itermax_
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific << itertol_
                << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                << scatra_pot_resnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << scatra_pot_incnorm / scatra_pot_dofnorm << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << manifold_pot_resnorm
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << manifold_pot_incnorm / manifold_pot_dofnorm << "   | " << std::endl;
    }

    // convergence check
    if ((scatra_pot_resnorm <= itertol_ and scatra_pot_incnorm / scatra_pot_dofnorm <= itertol_ and
            manifold_pot_resnorm <= itertol_ and
            manifold_pot_incnorm / manifold_pot_dofnorm <= itertol_) or
        (scatra_pot_resnorm < restol_ and manifold_pot_resnorm < restol_))
    {
      if (ssi_mono.Comm().MyPID() == 0)
      {
        std::cout << "+------------+-------------------+--------------+--------------+-----------"
                     "---+--------------+"
                  << std::endl
                  << std::endl;
      }
      converged = true;
    }
  }

  // warn if maximum number of iterations is reached without convergence
  if (init_pot_iternum == itermax_)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+--------------------------------------------------------------------------"
                   "------------------+"
                << std::endl;
      std::cout << "|            >>>>>> not converged!                                           "
                   "                 |"
                << std::endl;
      std::cout << "+--------------------------------------------------------------------------"
                   "------------------+"
                << std::endl;
    }
    converged = true;
  }

  return converged;
}