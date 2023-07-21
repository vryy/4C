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
#include "baci_ssi_monolithic_convcheck_strategies.H"
#include "baci_ssi_utils.H"

#include "baci_adapter_str_ssiwrapper.H"
#include "baci_adapter_scatra_base_algorithm.H"

#include "baci_scatra_timint_implicit.H"

#include "baci_linalg_mapextractor.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::SSIMono::ConvCheckStrategyBase::ConvCheckStrategyBase(
    const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
    )
    : itermax_(parameters.get<int>("ITEMAX")),
      itertol_(parameters.sublist("MONOLITHIC").get<double>("CONVTOL")),
      non_converged_steps_(),
      restol_(parameters.sublist("MONOLITHIC").get<double>("ABSTOLRES"))
{
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyBase::ExitNewtonRaphson(const SSI::SSIMono& ssi_mono)
{
  const auto norms = ComputeNorms(ssi_mono);
  const bool converged = CheckConvergence(ssi_mono, norms);
  const bool exit = ComputeExit(ssi_mono, converged);

  PrintNewtonIterationInformation(ssi_mono, converged, exit, norms);

  if (exit and !converged) non_converged_steps_.insert(ssi_mono.Step());

  return exit;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyBase::ComputeExit(
    const SSI::SSIMono& ssi_mono, const bool converged) const
{
  return (converged or ssi_mono.IterationCount() == itermax_);
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
      ->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::structure))
      ->Norm2(&incnorm);

  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::structure))
      ->Norm2(&resnorm);

  ssi_mono.StructureField()->Dispnp()->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyBase::PrintNonConvergedSteps(const int pid) const
{
  if (pid == 0 and not non_converged_steps_.empty())
  {
    std::cout << std::endl << "Non converged time steps: ";
    for (int step : non_converged_steps_)
    {
      std::cout << step;
      if (step != (*non_converged_steps_.end())) std::cout << ", ";
    }
    std::cout << std::endl << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyStd::GetAndCheckL2NormScaTra(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport))
      ->Norm2(&incnorm);

  ssi_mono.MapsSubProblems()
      ->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
          UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraField()->Phinp()->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
std::map<SSI::L2norm, double> SSI::SSIMono::ConvCheckStrategyStd::ComputeNorms(
    const SSI::SSIMono& ssi_mono) const
{
  double scatraincnorm = 0.0, scatraresnorm = 0.0, scatradofnorm = 0.0, structureincnorm = 0.0,
         structureresnorm = 0.0, structuredofnorm = 0.0;

  GetAndCheckL2NormScaTra(ssi_mono, scatraincnorm, scatraresnorm, scatradofnorm);
  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  return {{SSI::L2norm::scatraincnorm, scatraincnorm}, {SSI::L2norm::scatraresnorm, scatraresnorm},
      {SSI::L2norm::scatradofnorm, scatradofnorm},
      {SSI::L2norm::structureincnorm, structureincnorm},
      {SSI::L2norm::structureresnorm, structureresnorm},
      {SSI::L2norm::structuredofnorm, structuredofnorm}};
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyStd::CheckConvergence(
    const SSI::SSIMono& ssi_mono, const std::map<L2norm, double>& norms) const
{
  if (ssi_mono.IterationCount() > 1 and norms.at(L2norm::scatraresnorm) <= itertol_ and
      norms.at(L2norm::structureresnorm) <= itertol_ and
      norms.at(L2norm::scatraincnorm) / norms.at(L2norm::scatradofnorm) <= itertol_ and
      norms.at(L2norm::structureincnorm) / norms.at(L2norm::structuredofnorm) <= itertol_)
    return true;

  if (norms.at(L2norm::scatraresnorm) < restol_ and norms.at(L2norm::structureresnorm) < restol_)
    return true;

  return false;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyStd::PrintNewtonIterationInformation(
    const SSI::SSIMono& ssi_mono, const bool converged, const bool exit,
    const std::map<L2norm, double>& norms) const
{
  if (ssi_mono.Comm().MyPID() != 0) return;

  if (ssi_mono.IterationCount() == 1)
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
              << std::scientific << norms.at(L2norm::scatraresnorm) << "   |      --      | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(L2norm::structureresnorm) << "   |      --      | "
              << "(       --      , te = " << std::setw(10) << std::setprecision(3)
              << ssi_mono.dtele_ << ")" << std::endl;
  }
  else
  {
    std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
              << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
              << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
              << std::scientific << norms.at(L2norm::scatraresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(L2norm::scatraincnorm) / norms.at(L2norm::scatradofnorm) << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(L2norm::structureresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(L2norm::structureincnorm) / norms.at(L2norm::structuredofnorm)
              << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
              << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
              << std::endl;
  }

  if (exit and !converged)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                 "--------------+"
              << std::endl;
    std::cout << "|      Newton-Raphson method has not converged after a maximum number of "
              << std::setw(2) << itermax_ << " iterations!      |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+"
              << std::endl;
  }

  if (exit and converged)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+"
              << std::endl;
  }
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
              UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractOtherVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
              UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)))
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
              UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
              UTILS::SSIMaps::GetProblemPosition(Subproblem::scalar_transport)))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraField()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.ScaTraField()->Phinp())
      ->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
std::map<SSI::L2norm, double> SSI::SSIMono::ConvCheckStrategyElch::ComputeNorms(
    const SSI::SSIMono& ssi_mono) const
{
  double concincnorm = 0.0, concresnorm = 0.0, concdofnorm = 0.0, potdofnorm = 0.0,
         potincnorm = 0.0, potresnorm = 0.0, structuredofnorm = 0.0, structureresnorm = 0.0,
         structureincnorm = 0.0;

  GetAndCheckL2NormConc(ssi_mono, concincnorm, concresnorm, concdofnorm);
  GetAndCheckL2NormPot(ssi_mono, potincnorm, potresnorm, potdofnorm);
  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  return {
      {SSI::L2norm::concincnorm, concincnorm},
      {SSI::L2norm::concresnorm, concresnorm},
      {SSI::L2norm::concdofnorm, concdofnorm},
      {SSI::L2norm::potdofnorm, potdofnorm},
      {SSI::L2norm::potincnorm, potincnorm},
      {SSI::L2norm::potresnorm, potresnorm},
      {SSI::L2norm::structuredofnorm, structuredofnorm},
      {SSI::L2norm::structureresnorm, structureresnorm},
      {SSI::L2norm::structureincnorm, structureincnorm},
  };
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElch::CheckConvergence(
    const SSI::SSIMono& ssi_mono, const std::map<L2norm, double>& norms) const
{
  if (ssi_mono.IterationCount() > 1 and norms.at(SSI::L2norm::concresnorm) <= itertol_ and
      norms.at(SSI::L2norm::potresnorm) <= itertol_ and
      norms.at(SSI::L2norm::structureresnorm) <= itertol_ and
      norms.at(SSI::L2norm::concincnorm) / norms.at(SSI::L2norm::concdofnorm) <= itertol_ and
      norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) <= itertol_ and
      norms.at(SSI::L2norm::structureincnorm) / norms.at(SSI::L2norm::structuredofnorm) <= itertol_)
    return true;

  if (norms.at(SSI::L2norm::concresnorm) < restol_ and
      norms.at(SSI::L2norm::potresnorm) < restol_ and
      norms.at(SSI::L2norm::structureresnorm) < restol_)
    return true;

  return false;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElch::PrintNewtonIterationInformation(
    const SSI::SSIMono& ssi_mono, const bool converged, const bool exit,
    const std::map<L2norm, double>& norms) const
{
  if (ssi_mono.Comm().MyPID() != 0) return;

  if (ssi_mono.IterationCount() == 1)
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
              << std::scientific << norms.at(SSI::L2norm::concresnorm) << "   |      --      | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::potresnorm) << "   |      --      | " << std::setw(10)
              << std::setprecision(3) << std::scientific << norms.at(SSI::L2norm::structureresnorm)
              << "   |      --      | "
              << "(       --      , te = " << std::setw(10) << std::setprecision(3)
              << ssi_mono.dtele_ << ")" << std::endl;
  }
  else
  {
    std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
              << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
              << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
              << std::scientific << norms.at(SSI::L2norm::concresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::concincnorm) / norms.at(SSI::L2norm::concdofnorm) << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::potresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::structureresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::structureincnorm) / norms.at(SSI::L2norm::structuredofnorm)
              << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
              << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
              << std::endl;
  }

  if (exit and !converged)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+"
                 "--------------+--------------+--------------+"
              << std::endl;
    std::cout << "|                     Newton-Raphson method has not converged after a maximum "
                 "number of "
              << std::setw(2) << itermax_ << " iterations!                     |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+"
              << std::endl;
  }

  if (exit and converged)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+"
              << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElch::ExitNewtonRaphsonInitPotCalc(const SSI::SSIMono& ssi_mono)
{
  const auto norms = ComputeNorms(ssi_mono);
  const bool converged =
      (ssi_mono.IterationCount() > 1 and
          (norms.at(SSI::L2norm::potresnorm) <= itertol_ and
              norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) <= itertol_)) or
      norms.at(SSI::L2norm::potresnorm) < restol_;
  const bool exit = ComputeExit(ssi_mono, converged);
  if (exit and !converged) non_converged_steps_.insert(ssi_mono.Step());

  if (ssi_mono.Comm().MyPID() == 0)
  {
    if (ssi_mono.IterationCount() == 1)
    {
      // print header
      std::cout << "Calculating initial field for electric potential" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- pot-res ---|-- pot-inc ---|" << std::endl;

      // print only norm of residuals
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << norms.at(SSI::L2norm::potresnorm) << "   |      --      |"
                << std::endl;
    }

    else
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << norms.at(SSI::L2norm::potresnorm) << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific
                << norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) << "   |"
                << std::endl;
    }

    if (exit and !converged)
    {
      std::cout << "+--------------------------------------------------------------+" << std::endl;
      std::cout << "|            >>>>>> not converged!                             |" << std::endl;
      std::cout << "+--------------------------------------------------------------+" << std::endl;
    }
    if (exit and converged)
    {
      std::cout << "+------------+-------------------+--------------+--------------+" << std::endl
                << std::endl;
    }
  }

  return exit;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::GetAndCheckL2NormScaTraManifoldConc(
    const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const
{
  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractOtherVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
              UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractOtherVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
              UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)))
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
      ->ExtractCondVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Increment(),
              UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&incnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractCondVector(
          ssi_mono.MapsSubProblems()->ExtractVector(ssi_mono.ssi_vectors_->Residual(),
              UTILS::SSIMaps::GetProblemPosition(Subproblem::manifold)))
      ->Norm2(&resnorm);

  ssi_mono.ScaTraManifold()
      ->Splitter()
      ->ExtractCondVector(ssi_mono.ScaTraManifold()->Phinp())
      ->Norm2(&dofnorm);

  CheckL2Norm(incnorm, resnorm, dofnorm);
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
std::map<SSI::L2norm, double> SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::ComputeNorms(
    const SSI::SSIMono& ssi_mono) const
{
  double concincnorm = 0.0, concresnorm = 0.0, concdofnorm = 0.0, manifoldpotdofnorm = 0.0,
         manifoldpotincnorm = 0.0, manifoldpotresnorm = 0.0, manifoldconcdofnorm = 0.0,
         manifoldconcincnorm = 0.0, manifoldconcresnorm = 0.0, potdofnorm = 0.0, potincnorm = 0.0,
         potresnorm = 0.0, structuredofnorm = 0.0, structureresnorm = 0.0, structureincnorm = 0.0;

  GetAndCheckL2NormConc(ssi_mono, concincnorm, concresnorm, concdofnorm);
  GetAndCheckL2NormPot(ssi_mono, potincnorm, potresnorm, potdofnorm);

  GetAndCheckL2NormScaTraManifoldPot(
      ssi_mono, manifoldpotincnorm, manifoldpotresnorm, manifoldpotdofnorm);
  GetAndCheckL2NormScaTraManifoldConc(
      ssi_mono, manifoldconcincnorm, manifoldconcresnorm, manifoldconcdofnorm);

  GetAndCheckL2NormStructure(ssi_mono, structureincnorm, structureresnorm, structuredofnorm);

  return {
      {SSI::L2norm::concincnorm, concincnorm},
      {SSI::L2norm::concresnorm, concresnorm},
      {SSI::L2norm::concdofnorm, concdofnorm},
      {SSI::L2norm::manifoldpotdofnorm, manifoldpotdofnorm},
      {SSI::L2norm::manifoldpotincnorm, manifoldpotincnorm},
      {SSI::L2norm::manifoldpotresnorm, manifoldpotresnorm},
      {SSI::L2norm::manifoldconcdofnorm, manifoldconcdofnorm},
      {SSI::L2norm::manifoldconcincnorm, manifoldconcincnorm},
      {SSI::L2norm::manifoldconcresnorm, manifoldconcresnorm},
      {SSI::L2norm::potdofnorm, potdofnorm},
      {SSI::L2norm::potincnorm, potincnorm},
      {SSI::L2norm::potresnorm, potresnorm},
      {SSI::L2norm::structuredofnorm, structuredofnorm},
      {SSI::L2norm::structureresnorm, structureresnorm},
      {SSI::L2norm::structureincnorm, structureincnorm},
  };
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::CheckConvergence(
    const SSI::SSIMono& ssi_mono, const std::map<L2norm, double>& norms) const
{
  if (ssi_mono.IterationCount() > 1 and norms.at(SSI::L2norm::concresnorm) <= itertol_ and
      norms.at(SSI::L2norm::potresnorm) <= itertol_ and
      norms.at(SSI::L2norm::structureresnorm) <= itertol_ and
      norms.at(SSI::L2norm::concincnorm) / norms.at(SSI::L2norm::concdofnorm) <= itertol_ and
      norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) <= itertol_ and
      norms.at(SSI::L2norm::structureincnorm) / norms.at(SSI::L2norm::structuredofnorm) <=
          itertol_ and
      norms.at(SSI::L2norm::manifoldconcresnorm) <= itertol_ and
      norms.at(SSI::L2norm::manifoldpotresnorm) <= itertol_ and
      norms.at(SSI::L2norm::manifoldpotincnorm) / norms.at(SSI::L2norm::manifoldpotdofnorm) <=
          itertol_ and
      norms.at(SSI::L2norm::manifoldconcincnorm) / norms.at(SSI::L2norm::manifoldconcdofnorm) <=
          itertol_)
    return true;

  if (norms.at(SSI::L2norm::concresnorm) < restol_ and
      norms.at(SSI::L2norm::potresnorm) < restol_ and
      norms.at(SSI::L2norm::structureresnorm) < restol_ and
      norms.at(SSI::L2norm::manifoldconcresnorm) < restol_ and
      norms.at(SSI::L2norm::manifoldpotresnorm) < restol_)
    return true;

  return false;
}

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::PrintNewtonIterationInformation(
    const SSI::SSIMono& ssi_mono, const bool converged, const bool exit,
    const std::map<L2norm, double>& norms) const
{
  if (ssi_mono.Comm().MyPID() != 0) return;

  if (ssi_mono.IterationCount() == 1)
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
              << std::scientific << norms.at(SSI::L2norm::concresnorm) << "   |      --      | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::potresnorm) << "   |      --      | " << std::setw(10)
              << std::setprecision(3) << std::scientific << norms.at(SSI::L2norm::structureresnorm)
              << "   |      --      | " << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::manifoldconcresnorm) << "   |      --      |  "
              << std::scientific << norms.at(SSI::L2norm::manifoldpotresnorm)
              << "   |      --      | (       --      , te = " << std::setw(10)
              << std::setprecision(3) << ssi_mono.dtele_ << ")" << std::endl;
  }
  else
  {
    std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
              << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
              << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
              << std::scientific << norms.at(SSI::L2norm::concresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::concincnorm) / norms.at(SSI::L2norm::concdofnorm) << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::potresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::structureresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::structureincnorm) / norms.at(SSI::L2norm::structuredofnorm)
              << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::manifoldconcresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::manifoldconcincnorm) /
                     norms.at(SSI::L2norm::manifoldconcdofnorm)
              << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::manifoldpotresnorm) << "   | " << std::setw(10)
              << std::setprecision(3) << std::scientific
              << norms.at(SSI::L2norm::manifoldpotincnorm) /
                     norms.at(SSI::L2norm::manifoldpotdofnorm)
              << "   | (ts = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtsolve_
              << ", te = " << std::setw(10) << std::setprecision(3) << ssi_mono.dtele_ << ")"
              << std::endl;
  }

  if (exit and !converged)
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
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+--------------+--------------+--------"
                 "------+--------------+"
              << std::endl;
  }

  if (exit and converged)
  {
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--"
                 "------------+--------------+--------------+--------------+--------------+--------"
                 "------+--------------+"
              << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool SSI::SSIMono::ConvCheckStrategyElchScaTraManifold::ExitNewtonRaphsonInitPotCalc(
    const SSI::SSIMono& ssi_mono)
{
  const auto norms = ComputeNorms(ssi_mono);
  const bool converged =
      (ssi_mono.IterationCount() > 1 and norms.at(SSI::L2norm::potresnorm) <= itertol_ and
          norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) <= itertol_ and
          norms.at(SSI::L2norm::manifoldpotresnorm) <= itertol_ and
          norms.at(SSI::L2norm::manifoldpotincnorm) / norms.at(SSI::L2norm::manifoldpotdofnorm) <=
              itertol_) or
      (norms.at(SSI::L2norm::potresnorm) < restol_ and
          norms.at(SSI::L2norm::manifoldpotresnorm) < restol_);
  const bool exit = ComputeExit(ssi_mono, converged);
  if (exit and !converged) non_converged_steps_.insert(ssi_mono.Step());

  if (ssi_mono.Comm().MyPID() == 0)
  {
    if (ssi_mono.IterationCount() == 1)
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
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << norms.at(SSI::L2norm::potresnorm) << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific
                << norms.at(SSI::L2norm::manifoldpotresnorm) << "   |      --      | " << std::endl;
    }
    else
    {
      std::cout << "|  " << std::setw(3) << ssi_mono.IterationCount() << "/" << std::setw(3)
                << itermax_ << "   | " << std::setw(10) << std::setprecision(3) << std::scientific
                << itertol_ << "[L_2 ]  | " << std::setw(10) << std::setprecision(3)
                << std::scientific << norms.at(SSI::L2norm::potresnorm) << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific
                << norms.at(SSI::L2norm::potincnorm) / norms.at(SSI::L2norm::potdofnorm) << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific
                << norms.at(SSI::L2norm::manifoldpotresnorm) << "   | " << std::setw(10)
                << std::setprecision(3) << std::scientific
                << norms.at(SSI::L2norm::manifoldpotincnorm) /
                       norms.at(SSI::L2norm::manifoldpotdofnorm)
                << "   | " << std::endl;
    }

    // warn if maximum number of iterations is reached without convergence
    if (exit and !converged)
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
    if (exit and converged)
    {
      std::cout << "+------------+-------------------+--------------+--------------+-----------"
                   "---+--------------+"
                << std::endl
                << std::endl;
    }
  }

  return exit;
}