/*----------------------------------------------------------------------*/
/*! \file

\brief strategies for Newton-Raphson convergence check for scalar transport problems

\level 1


 */
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_TIMINT_CONVCHECK_STRATEGIES_HPP
#define FOUR_C_SCATRA_TIMINT_CONVCHECK_STRATEGIES_HPP

#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"

// forward declaration
namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

namespace SCATRA
{
  // forward declaration
  class ScaTraTimIntImpl;

  /*!
  \brief abstract strategy for Newton-Raphson convergence check for scalar transport problems

   To keep the scalar transport time integrator class and derived classes as plain as possible,
   the convergence check for the Newton-Raphson iteration has been encapsulated within separate
   strategy classes. Every specific convergence check strategy (e.g., for standard scalar transport
   or electrochemistry problems with or without scatra-scatra interface coupling involving Lagrange
   multipliers) computes, checks, and outputs different relevant vector norms and is implemented
   in a subclass derived from an abstract, purely virtual interface class.
  */

  class ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyBase(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
    );

    //! destructor
    virtual ~ConvCheckStrategyBase() = default;

    //! perform convergence check for Newton-Raphson iteration
    virtual bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const = 0;

    //! perform convergence check for outer iteration in partitioned coupling schemes
    virtual bool AbortOuterIter(
        const ScaTraTimIntImpl& scatratimint  //!< scalar transport time integrator
    ) const
    {
      dserror("Not implemented in base class!");
      return false;
    }

   protected:
    //! maximum number of Newton-Raphson iteration steps
    const int itmax_;

    //! relative tolerance for Newton-Raphson scheme
    const double ittol_;

    //! absolute tolerance for residuals
    const double abstolres_;

    //! maximum number of outer iteration steps in partitioned coupling schemes
    const int itmax_outer_;

    //! relative tolerance for outer iteration in partitioned coupling schemes
    const double ittol_outer_;

   private:
    //! copy constructor
    ConvCheckStrategyBase(const ConvCheckStrategyBase& old);
  };  // class ConvCheckStrategyBase

  /*!
  \brief specific strategy for Newton-Raphson convergence check for standard scalar transport
  problems without scatra-scatra interface coupling involving Lagrange multipliers

  */

  class ConvCheckStrategyStd : public ConvCheckStrategyBase
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


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

    //! perform convergence check for outer iteration in partitioned coupling schemes
    bool AbortOuterIter(const ScaTraTimIntImpl& scatratimint  //!< scalar transport time integrator
    ) const override;

   private:
    //! copy constructor
    ConvCheckStrategyStd(const ConvCheckStrategyStd& old);
  };  // class ConvCheckStrategyStd

  /*!
  \brief specific strategy for Newton-Raphson convergence check for micro scale in multi-scale
  simulations of scalar transport problems

  */

  class ConvCheckStrategyStdMicroScale : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyStdMicroScale(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

   private:
    //! copy constructor
    ConvCheckStrategyStdMicroScale(const ConvCheckStrategyStdMicroScale& old);
  };  // class ConvCheckStrategyStdMicroScale

  /*!
  \brief specific strategy for Newton-Raphson convergence check for electrochemistry problems
  without scatra-scatra interface coupling involving Lagrange multipliers

  */

  class ConvCheckStrategyStdElch : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyStdElch(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

   private:
    //! copy constructor
    ConvCheckStrategyStdElch(const ConvCheckStrategyStdElch& old);
  };  // class ConvCheckStrategyStdElch

  /*!
  \brief specific strategy for Newton-Raphson convergence check for standard scalar transport
  problems with scatra-scatra interface coupling involving Lagrange multipliers

  */

  class ConvCheckStrategyS2ILM : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyS2ILM(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

   private:
    //! copy constructor
    ConvCheckStrategyS2ILM(const ConvCheckStrategyS2ILM& old);
  };  // class ConvCheckStrategyS2ILM

  /*!
  \brief specific strategy for Newton-Raphson convergence check for electrochemistry problems with
  scatra-scatra interface coupling involving Lagrange multipliers

  */

  class ConvCheckStrategyS2ILMElch : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyS2ILMElch(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

   private:
    //! copy constructor
    ConvCheckStrategyS2ILMElch(const ConvCheckStrategyS2ILMElch& old);
  };  // class ConvCheckStrategyS2ILMElch

  /*!
  \brief specific strategy for Newton-Raphson convergence check for macro scale in multi-scale
  simulations of electrochemistry problems

  */

  class ConvCheckStrategyStdMacroScaleElch : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyStdMacroScaleElch(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

    //! perform convergence check for outer iteration in partitioned coupling schemes
    bool AbortOuterIter(const ScaTraTimIntImpl& scatratimint  //!< scalar transport time integrator
    ) const override;

   private:
    //! copy constructor
    ConvCheckStrategyStdMacroScaleElch(const ConvCheckStrategyStdMacroScaleElch& old);
  };  // class ConvCheckStrategyStdMacroScaleElch

  /*!
  \brief specific strategy for Newton-Raphson convergence check for Multiphase_Poroelasticity_ScaTra
  problems
         || delta phi|| / || phi || < tol    and    || res || / sqrt(totalnumdof) < tol

  */

  class ConvCheckStrategyPoroMultiphaseScatra : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyPoroMultiphaseScatra(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

   private:
    //! copy constructorexit
    ConvCheckStrategyPoroMultiphaseScatra(const ConvCheckStrategyPoroMultiphaseScatra& old);
  };  // class ConvCheckStrategyPoroMultiphaseScatra
  /*!
  \brief specific strategy for Newton-Raphson convergence check for Multiphase_Poroelasticity_ScaTra
  problems with artery coupling
         || delta phi|| / || phi || < tol    and    || res || / sqrt(totalnumdof) < tol
         where phi are both artery_scatra and scatra DOFs and res is the combined residual

  */

  class ConvCheckStrategyPoroMultiphaseScatraArtMeshTying : public ConvCheckStrategyBase
  {
   public:
    //! constructor
    ConvCheckStrategyPoroMultiphaseScatraArtMeshTying(
        const Teuchos::ParameterList& parameters  //!< parameter list for Newton-Raphson iteration
        )
        : ConvCheckStrategyBase(parameters)
    {
      return;
    };


    //! perform convergence check for Newton-Raphson iteration
    bool AbortNonlinIter(
        const ScaTraTimIntImpl& scatratimint,  //!< scalar transport time integrator
        double& actresidual                    //!< return maximum current residual value
    ) const override;

   private:
    //! copy constructorexit
    ConvCheckStrategyPoroMultiphaseScatraArtMeshTying(
        const ConvCheckStrategyPoroMultiphaseScatra& old);
  };  // class ConvCheckStrategyPoroMultiphaseScatra
}  // namespace SCATRA
BACI_NAMESPACE_CLOSE

#endif
