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
#ifndef BACI_SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H
#define BACI_SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H

#include "baci_config.H"

#include "baci_ssi_monolithic.H"

BACI_NAMESPACE_OPEN

namespace SSI
{
  enum class L2norm
  {
    concdofnorm,
    concincnorm,
    concresnorm,
    manifoldpotdofnorm,
    manifoldpotincnorm,
    manifoldpotresnorm,
    manifoldconcdofnorm,
    manifoldconcincnorm,
    manifoldconcresnorm,
    potdofnorm,
    potincnorm,
    potresnorm,
    scatradofnorm,
    scatraincnorm,
    scatraresnorm,
    structuredofnorm,
    structureincnorm,
    structureresnorm
  };

  class SSIMono::ConvCheckStrategyBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~ConvCheckStrategyBase() = default;

    ConvCheckStrategyBase(const Teuchos::ParameterList& global_time_parameters);

    //! Check, if Newton-Raphson has converged and print residuals and increments to screen
    bool ExitNewtonRaphson(const SSI::SSIMono& ssi_mono);

    //! Check, if Newton-Raphson of initial potential calculation has converged and print residuals
    //! and increments to screen
    virtual bool ExitNewtonRaphsonInitPotCalc(const SSI::SSIMono& ssi_mono) = 0;

    //! print all time steps that have not been converged
    void PrintNonConvergedSteps(int pid) const;

   protected:
    //! get L2 norms from structure field and check if they are reasonable
    //!
    //! \param ssi_mono  ssi time integration
    //! \param incnorm   (out) L2 norm of increment
    //! \param resnorm   (out) L2 norm of residual
    //! \param dofnorm   (out) L2 norm of displacement
    void GetAndCheckL2NormStructure(
        const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const;

    //! check, if L2 norm is inf or nan. For dofnom: check if it is numerical zero
    //!
    //! \param incnorm  (out) L2 norm of increment
    //! \param resnorm  (out) L2 norm of residual
    //! \param dofnorm  (out) L2 norm of state
    void CheckL2Norm(double& incnorm, double& resnorm, double& dofnorm) const;

    //! decide, if Newton loop should be exited, if converged or maximum number of steps are reached
    //!
    //! \param ssi_mono   ssi time integration
    //! \param converged  convergence of Newton loop
    //! \return  decision on exit
    bool ComputeExit(const SSI::SSIMono& ssi_mono, bool converged) const;

    //! maximum number of Newton-Raphson iteration steps
    const int itermax_;

    //! relative tolerance for Newton-Raphson iteration
    const double itertol_;

    //! time steps that have not converged
    std::set<int> non_converged_steps_;

    //! absolute tolerance for residual vectors
    const double restol_;

   private:
    //! decide, if Newton loop has converged
    //!
    //! \param ssi_mono  ssi time integration
    //! \param norms     L2 norms of residual, increment, and state
    //! \return  decision on convergence
    virtual bool CheckConvergence(
        const SSI::SSIMono& ssi_mono, const std::map<SSI::L2norm, double>& norms) const = 0;

    //! compute L2 norms and fill into a map
    //!
    //! \param ssi_mono   ssi time integration
    //! \return  map with string identifier of norm as key and its value
    virtual std::map<SSI::L2norm, double> ComputeNorms(const SSI::SSIMono& ssi_mono) const = 0;

    //! print convergence table to screen
    //!
    //! \param ssi_mono   ssi time integration
    //! \param converged  convergence of Newton loop?
    //! \param exit       exit of Newton loop?
    //! \param norms      computed L2 norms
    virtual void PrintNewtonIterationInformation(const SSI::SSIMono& ssi_mono, bool converged,
        bool exit, const std::map<L2norm, double>& norms) const = 0;
  };

  class SSIMono::ConvCheckStrategyStd : public SSIMono::ConvCheckStrategyBase
  {
   public:
    ConvCheckStrategyStd(const Teuchos::ParameterList& global_time_parameters)
        : ConvCheckStrategyBase(global_time_parameters){};

    bool ExitNewtonRaphsonInitPotCalc(const SSI::SSIMono& ssi_mono) override
    {
      dserror("Caluclation of initial potential only for Elch");
      return {};
    }

   protected:
    //! get L2 norms from scatra field and check if they are reasonable
    void GetAndCheckL2NormScaTra(
        const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const;

   private:
    bool CheckConvergence(
        const SSI::SSIMono& ssi_mono, const std::map<L2norm, double>& norms) const override;

    std::map<SSI::L2norm, double> ComputeNorms(const SSI::SSIMono& ssi_mono) const override;

    void PrintNewtonIterationInformation(const SSI::SSIMono& ssi_mono, const bool converged,
        const bool exit, const std::map<L2norm, double>& norms) const override;
  };


  class SSIMono::ConvCheckStrategyElch : public SSIMono::ConvCheckStrategyBase
  {
   public:
    ConvCheckStrategyElch(const Teuchos::ParameterList& global_time_parameters)
        : ConvCheckStrategyBase(global_time_parameters){};

    bool ExitNewtonRaphsonInitPotCalc(const SSI::SSIMono& ssi_mono) override;

   protected:
    //! get L2 norms from concentration field and check if they are reasonable
    void GetAndCheckL2NormConc(
        const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const;

    //! get L2 norms from potential field and check if they are reasonable
    void GetAndCheckL2NormPot(
        const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const;

   private:
    bool CheckConvergence(
        const SSI::SSIMono& ssi_mono, const std::map<L2norm, double>& norms) const override;

    std::map<SSI::L2norm, double> ComputeNorms(const SSI::SSIMono& ssi_mono) const override;

    void PrintNewtonIterationInformation(const SSI::SSIMono& ssi_mono, const bool converged,
        const bool exit, const std::map<L2norm, double>& norms) const override;
  };

  class SSIMono::ConvCheckStrategyElchScaTraManifold : public SSIMono::ConvCheckStrategyElch
  {
   public:
    ConvCheckStrategyElchScaTraManifold(const Teuchos::ParameterList& global_time_parameters)
        : ConvCheckStrategyElch(global_time_parameters){};

    bool ExitNewtonRaphsonInitPotCalc(const SSI::SSIMono& ssi_mono) override;

   protected:
    //! get L2 norms from concentration field and check if they are reasonable
    void GetAndCheckL2NormScaTraManifoldConc(
        const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const;

    //! get L2 norms from potential field and check if they are reasonable
    void GetAndCheckL2NormScaTraManifoldPot(
        const SSI::SSIMono& ssi_mono, double& incnorm, double& resnorm, double& dofnorm) const;

   private:
    bool CheckConvergence(
        const SSI::SSIMono& ssi_mono, const std::map<L2norm, double>& norms) const override;

    std::map<SSI::L2norm, double> ComputeNorms(const SSI::SSIMono& ssi_mono) const override;

    void PrintNewtonIterationInformation(const SSI::SSIMono& ssi_mono, const bool converged,
        const bool exit, const std::map<L2norm, double>& norms) const override;
  };
}  // namespace SSI
BACI_NAMESPACE_CLOSE

#endif  // SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H
