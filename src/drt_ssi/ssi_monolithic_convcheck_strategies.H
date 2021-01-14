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
#ifndef SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H
#define SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H

#include "ssi_monolithic.H"

namespace SSI
{
  /*!
  \brief abstract strategy for Newton-Raphson convergence check for monolithic scalar-structure
  interaction problems

  */

  class SSIMono::ConvCheckStrategyBase
  {
   public:
    ConvCheckStrategyBase(const Teuchos::ParameterList& global_time_parameters);

    virtual ~ConvCheckStrategyBase() = default;

    /*!
     * \brief Check, if Newton-Raphson has converged and print residuals and increments to screen
     *
     * @param[in] ssi_mono     time integration algorithm for monolithic SSI
     * */
    virtual bool ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const = 0;

   protected:
    //! maximum number of Newton-Raphson iteration steps
    const int itermax_;

    //! relative tolerance for Newton-Raphson iteration
    const double itertol_;

    //! absolute tolerance for residual vectors
    const double restol_;

   private:
    //! copy constructor
    ConvCheckStrategyBase(const ConvCheckStrategyBase& old);
  };  // class ConvCheckStrategyBase


  /*!
  \brief specific strategy for Newton-Raphson convergence check for monolithic scalar-structure
  interaction problems involving standard scalar transport

  */

  class SSIMono::ConvCheckStrategyStd : public SSIMono::ConvCheckStrategyBase
  {
   public:
    ConvCheckStrategyStd(const Teuchos::ParameterList& global_time_parameters)
        : ConvCheckStrategyBase(global_time_parameters){};

    ~ConvCheckStrategyStd() override = default;

    bool ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const override;

   private:
    //! copy constructor
    ConvCheckStrategyStd(const ConvCheckStrategyStd& old);
  };  // class ConvCheckStrategyStd


  /*!
  \brief specific strategy for Newton-Raphson convergence check for monolithic scalar-structure
  interaction problems involving electrochemistry

  */

  class SSIMono::ConvCheckStrategyElch : public SSIMono::ConvCheckStrategyBase
  {
   public:
    ConvCheckStrategyElch(const Teuchos::ParameterList& global_time_parameters)
        : ConvCheckStrategyBase(global_time_parameters){};

    ~ConvCheckStrategyElch() override = default;

    bool ExitNewtonRaphson(const SSI::SSIMono& ssi_mono) const override;

   private:
    //! copy constructor
    ConvCheckStrategyElch(const ConvCheckStrategyElch& old);
  };  // class ConvCheckStrategyElch
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H
