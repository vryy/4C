/*----------------------------------------------------------------------*/
/*! \file

\brief Standard solution strategy for electrochemistry problems (without meshtying)


\level 2

*----------------------------------------------------------------------*/
#ifndef FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_STD_ELCH_HPP
#define FOUR_C_SCATRA_TIMINT_MESHTYING_STRATEGY_STD_ELCH_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_elch.hpp"
#include "4C_scatra_timint_meshtying_strategy_std.hpp"

FOUR_C_NAMESPACE_OPEN

namespace SCATRA
{
  /*!
  \brief Standard solution strategy for electrochemistry problems (without meshtying)

  To keep the scalar transport time integrator class and derived classes as plain as possible,
  several algorithmic parts have been encapsulated within separate meshtying strategy classes.
  These algorithmic parts include initializing the system matrix and other relevant objects,
  computing meshtying residual terms and their linearizations, and solving the resulting
  linear system of equations. By introducing a hierarchy of strategies for these algorithmic
  parts, a bunch of unhandy if-else selections within the time integrator classes themselves
  can be circumvented. This class contains the standard solution strategy for electrochemistry
  problems without meshtying.

  */

  class MeshtyingStrategyStdElch : public MeshtyingStrategyStd
  {
   public:
    //! constructor
    explicit MeshtyingStrategyStdElch(SCATRA::ScaTraTimIntElch* elchtimint);


    bool system_matrix_initialization_needed() const override { return true; }

    Teuchos::RCP<CORE::LINALG::SparseOperator> InitSystemMatrix() const override;
    //@}

   private:
    //! copy constructor
    MeshtyingStrategyStdElch(const MeshtyingStrategyStdElch& old);

    //! return pointer to elch time integrator after cast
    SCATRA::ScaTraTimIntElch* ElchTimInt() const
    {
      return dynamic_cast<SCATRA::ScaTraTimIntElch*>(scatratimint_);
    };

    //! instantiate strategy for Newton-Raphson convergence check
    void init_conv_check_strategy() override;
  };  // class MeshtyingStrategyStdElch
}  // namespace SCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
