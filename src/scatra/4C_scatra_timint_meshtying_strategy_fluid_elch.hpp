/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid-fluid meshtying strategy for electrochemistry problems

\level 2


*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_ELCH_HPP
#define FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_FLUID_ELCH_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_elch.hpp"
#include "4C_scatra_timint_meshtying_strategy_fluid.hpp"

FOUR_C_NAMESPACE_OPEN

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


    //! initialize meshtying objects
    void InitMeshtying() override;

    //! setup meshtying objects
    void setup_meshtying() override;

   private:
    //! copy constructor
    MeshtyingStrategyFluidElch(const MeshtyingStrategyFluidElch& old);

    //! return pointer to elch time integrator after cast
    SCATRA::ScaTraTimIntElch* elch_tim_int() const
    {
      return dynamic_cast<SCATRA::ScaTraTimIntElch*>(scatratimint_);
    };

    //! instantiate strategy for Newton-Raphson convergence check
    void init_conv_check_strategy() override;
  };  // class MeshtyingStrategyFluidElch
}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif
