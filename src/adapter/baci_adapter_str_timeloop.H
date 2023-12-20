/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_STR_TIMELOOP_H
#define BACI_ADAPTER_STR_TIMELOOP_H

#include "baci_config.H"

#include "baci_adapter_str_wrapper.H"

BACI_NAMESPACE_OPEN

namespace ADAPTER
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
    int Integrate() override;

    /// wrapper for things that should be done before PrepareTimeStep is called
    void PrePredict() override{};

    /// wrapper for things that should be done before solving the nonlinear iterations
    void PreSolve() override{};

    /// wrapper for things that should be done before updating
    void PreUpdate() override{};

    /// wrapper for things that should be done after solving the update
    void PostUpdate() override{};

    /// wrapper for things that should be done after the output
    void PostOutput() override{};
  };

}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
