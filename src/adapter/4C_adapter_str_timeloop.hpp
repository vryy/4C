/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STR_TIMELOOP_HPP
#define FOUR_C_ADAPTER_STR_TIMELOOP_HPP

#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /*! \brief Time loop for stuctural simulations
   *
   *  This is a wrapper for the structural time integration which gives
   *  fine-grained access into the time loop by various pre- and post-operators.
   *
   *  To perform such pre- and post-operations, just derive from this class and
   *  overload the respective pre-/post-operator.
   *
   *  Implementations of pre-/post-operators in this class have to remain empty!
   */
  class StructureTimeLoop : public StructureWrapper
  {
   public:
    /// constructor
    explicit StructureTimeLoop(Teuchos::RCP<Structure> structure) : StructureWrapper(structure) {}

    /// actual time loop
    int integrate() override;

    /// wrapper for things that should be done before prepare_time_step is called
    void pre_predict() override{};

    /// wrapper for things that should be done before solving the nonlinear iterations
    void pre_solve() override{};

    /// wrapper for things that should be done before updating
    void pre_update() override{};

    /// wrapper for things that should be done after solving the update
    void post_update() override{};

    /// wrapper for things that should be done after the output
    void post_output() override{};
  };

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
