/*----------------------------------------------------------------------*/
/*!
\file ad_str_timeloop.H

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

\level 0

<pre>
\maintainer Michael Hiermeier
            hiermeier@lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_TIMELOOP_H
#define ADAPTER_STRUCTURE_TIMELOOP_H

#include "ad_str_wrapper.H"

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
    virtual int Integrate();

    /// wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict(){};

    /// wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict(){};

    /// wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve(){};

    /// wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve(){};

    /// wrapper for things that should be done before updating
    virtual void PreUpdate(){};

    /// wrapper for things that should be done after solving the update
    virtual void PostUpdate(){};

    /// wrapper for things that should be done before the output
    virtual void PreOutput(){};

    /// wrapper for things that should be done after the output
    virtual void PostOutput(){};
  };

}  // namespace ADAPTER

#endif
