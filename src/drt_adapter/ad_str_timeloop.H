/*----------------------------------------------------------------------*/
/*!
\file ad_str_timeloop.H

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_TIMELOOP_H
#define ADAPTER_STRUCTURE_TIMELOOP_H

#include "ad_str_wrapper.H"

namespace ADAPTER
{

  /// Just wrap, do nothing new, meant to be derived from
  class StructureTimeLoop : public StructureWrapper
  {
  public:

    /// constructor
    explicit StructureTimeLoop(Teuchos::RCP<Structure> structure) : StructureWrapper(structure) {}

    /// actual time loop
    virtual int Integrate();

    /// empty wrapper for things that should be done before PrepareTimeStep is called
    virtual void PrePredict(){};

    /// empty wrapper for things that should be done after PrepareTimeStep is called
    virtual void PostPredict(){};

    /// empty wrapper for things that should be done before solving the nonlinear iterations
    virtual void PreSolve(){};

    /// empty wrapper for things that should be done after solving the nonlinear iterations
    virtual void PostSolve(){};

    /// empty wrapper for things that should be done before updating
    virtual void PreUpdate(){};

    /// empty wrapper for things that should be done after solving the update
    virtual void PostUpdate(){};

    /// empty wrapper for things that should be done before the output
    virtual void PreOutput(){};

    /// empty wrapper for things that should be done after the output
    virtual void PostOutput(){};

    /// empty wrapper for things that should be done after the actual time loop is finished
    virtual void PostTimeLoop(){};

  };

}

#endif
