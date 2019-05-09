/*----------------------------------------------------------------------*/
/*!
\brief strategies for Newton-Raphson convergence check for monolithic scalar-structure interaction
problems

To keep the time integrator class for monolithic scalar-structure interaction problems as plain as
possible, the convergence check for the Newton-Raphson iteration has been encapsulated within
separate strategy classes. Every specific convergence check strategy (e.g., for monolithic
scalar-structure interaction problems involving standard scalar transport or electrochemistry)
computes, checks, and outputs different relevant vector norms and is implemented in a subclass
derived from an abstract, purely virtual interface class.

\level 2

\maintainer Christoph Schmidt

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

  \author fang (fang@lnm.mw.tum.de)
  */

  class SSI_Mono::ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyBase(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
    );

    //! destructor
    virtual ~ConvCheckStrategyBase() { return; };

    //! check termination criterion for Newton-Raphson iteration
    virtual bool ExitNewtonRaphson(const SSI::SSI_Mono&
            timint  //!< time integrator for monolithic scalar-structure interaction
        ) const = 0;

   protected:
    //! maximum number of Newton-Raphson iteration steps
    const unsigned itermax_;

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

  \author fang (fang@lnm.mw.tum.de)
  */

  class SSI_Mono::ConvCheckStrategyStd : public SSI_Mono::ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyStd(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };

    //! destructor
    virtual ~ConvCheckStrategyStd() { return; };

    //! check termination criterion for Newton-Raphson iteration
    bool ExitNewtonRaphson(const SSI::SSI_Mono&
            timint  //!< time integrator for monolithic scalar-structure interaction
        ) const;

   private:
    //! copy constructor
    ConvCheckStrategyStd(const ConvCheckStrategyStd& old);
  };  // class ConvCheckStrategyStd


  /*!
  \brief specific strategy for Newton-Raphson convergence check for monolithic scalar-structure
  interaction problems involving electrochemistry

  \author fang (fang@lnm.mw.tum.de)
  */

  class SSI_Mono::ConvCheckStrategyElch : public SSI_Mono::ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyElch(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };

    //! destructor
    virtual ~ConvCheckStrategyElch() { return; };

    //! check termination criterion for Newton-Raphson iteration
    bool ExitNewtonRaphson(const SSI::SSI_Mono&
            timint  //!< time integrator for monolithic scalar-structure interaction
        ) const;

   private:
    //! copy constructor
    ConvCheckStrategyElch(const ConvCheckStrategyElch& old);
  };  // class ConvCheckStrategyElch
}  // namespace SSI
#endif  // #ifndef SSI_MONOLITHIC_CONVCHECK_STRATEGIES_H
