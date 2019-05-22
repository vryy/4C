/*!----------------------------------------------------------------------

\brief Fluid-fluid meshtying strategy for electrochemistry problems

\level 2

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_ELCH_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_ELCH_H

#include "../drt_scatra/scatra_timint_elch.H"

#include "scatra_timint_meshtying_strategy_fluid.H"

namespace SCATRA
{
  /*!
  \brief Fluid-fluid meshtying strategy for electrochemistry problems

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the fluid-fluid meshtying strategy for electrochemistry
  problems.

  */

  class MeshtyingStrategyFluidElch : public MeshtyingStrategyFluid
  {
   public:
    //! constructor
    explicit MeshtyingStrategyFluidElch(SCATRA::ScaTraTimIntElch* elchtimint);

    //! destructor
    virtual ~MeshtyingStrategyFluidElch() { return; };

    //! initialize meshtying objects
    void InitMeshtying();

    //! initialize system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;

   private:
    //! copy constructor
    MeshtyingStrategyFluidElch(const MeshtyingStrategyFluidElch& old);

    //! return pointer to elch time integrator after cast
    SCATRA::ScaTraTimIntElch* ElchTimInt() const
    {
      return dynamic_cast<SCATRA::ScaTraTimIntElch*>(scatratimint_);
    };

    //! instantiate strategy for Newton-Raphson convergence check
    void InitConvCheckStrategy();
  };  // class MeshtyingStrategyFluidElch
}  // namespace SCATRA
#endif  // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_ELCH_H
