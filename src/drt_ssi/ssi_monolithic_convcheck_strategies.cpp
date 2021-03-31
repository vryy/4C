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
void SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::GetAndCheckL2NormScaTraManifold(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  //! compute L2 norm of manifold state vector
  ssi_mono.ScaTraManifold()->Phinp()->Norm2(&dofnorm);

  //! compute L2 norm of manifold increment vector
  ssi_mono.MapsSubProblems()
      ->ExtractVector(
          ssi_mono.ssi_vectors_->Increment(), ssi_mono.GetProblemPosition(Subproblem::manifold))
      ->Norm2(&incnorm);

  //! compute L2 norm of manifold residual vector
  ssi_mono.MapsSubProblems()
      ->ExtractVector(
          ssi_mono.ssi_vectors_->Residual(), ssi_mono.GetProblemPosition(Subproblem::manifold))
      ->Norm2(&resnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::ExitNewtonRaphson(
    const SSI::SSIMono& ssi_mono) const
{
  // initialize exit flag
  bool exit(false);

  double concincnorm(0.0), concresnorm(0.0), concdofnorm(0.0), manifolddofnorm(0.0),
      manifoldincnorm(0.0), manifoldresnorm(0.0), potdofnorm(0.0), potincnorm(0.0), potresnorm(0.0),
      structuredofnorm(0.0), structureresnorm(0.0), structureincnorm(0.0);

  GetAndCheckL2NormConc(ssi_mono, concincnorm, concresnorm, concdofnorm);
  GetAndCheckL2NormScaTraManifold(ssi_mono, manifoldincnorm, manifoldresnorm, manifolddofnorm);
  GetAndCheckL2NormPot(ssi_mono, potincnorm, potresnorm, potdofnorm);
  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  // first Newton-Raphson iteration
  if (ssi_mono.IterationCount() == 1)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+--------------+--------------+"
                << std::endl;
      std::cout << "|- step/max -|- tolerance[norm] -|-- conc-res --|-- conc-inc --|-- pot-res "
                   "---|-- pot-inc ---|- struct-res -|- struct-inc -|-  surf-res  -|-  surf-inc  -|"
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
                << std::scientific << manifoldresnorm << "   |      --      | "
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
                << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << manifoldresnorm << "   | " << std::setw(10) << std::setprecision(3)
                << std::scientific << manifoldincnorm / manifolddofnorm
                << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
                << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
                << std::endl;
    }

    // convergence check
    if (concresnorm <= itertol_ and potresnorm <= itertol_ and structureresnorm <= itertol_ and
        manifoldresnorm <= itertol_ and concincnorm / concdofnorm <= itertol_ and
        potincnorm / potdofnorm <= itertol_ and structureincnorm / structuredofnorm <= itertol_ and
        manifoldincnorm / manifolddofnorm <= itertol_)
      // exit Newton-Raphson iteration upon convergence
      exit = true;
  }

  // exit Newton-Raphson iteration when residuals are small enough to prevent unnecessary additional
  // solver calls
  if (concresnorm < restol_ and potresnorm < restol_ and structureresnorm < restol_ and
      manifoldresnorm < restol_)
    exit = true;

  // print warning to screen if maximum number of Newton-Raphson iterations is reached without
  // convergence

  if (ssi_mono.IterationCount() == itermax_ and !exit)
  {
    if (ssi_mono.Comm().MyPID() == 0)
    {
      std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                   "--------------+--------------+--------------+--------------+--------------+"
                << std::endl;
      std::cout << "|                     Newton-Raphson method has not converged after a maximum "
                   "number of "
                << std::setw(2) << itermax_
                << " iterations!                                                   |" << std::endl;
    }

    // proceed to next time step
    exit = true;
  }

  // print finish line of convergence table to screen
  if (exit and ssi_mono.Comm().MyPID() == 0)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+--------------+--------------+"
              << std::endl;
  }

  return exit;
}