/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.H

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_H

#include "../drt_scatra/scatra_timint_elch.H"

#include "scatra_timint_meshtying_strategy_s2i.H"

namespace SCATRA
{

/*!
\brief Scatra-scatra interface coupling strategy for electrochemistry problems

To keep the scalar transport time integrator class and derived classes as plain as possible,
several algorithmic parts have been encapsulated within separate meshtying strategy classes.
These algorithmic parts include initializing the system matrix and other relevant objects,
computing meshtying residual terms and their linearizations, and solving the resulting
linear system of equations. By introducing a hierarchy of strategies for these algorithmic
parts, a bunch of unhandy if-else selections within the time integrator classes themselves
can be circumvented. This class contains the scatra-scatra interface coupling strategy for
electrochemistry problems.

\author fang (fang@lnm.mw.tum.de)
*/

class MeshtyingStrategyS2IElch : public MeshtyingStrategyS2I
{
  public:

    //! constructor
    explicit MeshtyingStrategyS2IElch(
        SCATRA::ScaTraTimIntElch*       elchtimint,   //! elch time integrator
        const Teuchos::ParameterList&   parameters    //! input parameters for scatra-scatra interface coupling
        );

    //! destructor
    virtual ~MeshtyingStrategyS2IElch(){return;};

    //! compute meshtying residual terms and their linearizations
    void EvaluateMeshtying() const;

  private:

    //! copy constructor
    MeshtyingStrategyS2IElch(const MeshtyingStrategyS2IElch& old);

    //! return pointer to elch time integrator after cast
    SCATRA::ScaTraTimIntElch* ElchTimInt() const {return dynamic_cast<SCATRA::ScaTraTimIntElch*>(scatratimint_);};
}; // class MeshtyingStrategyS2IElch
} // namespace SCATRA
#endif // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_S2I_ELCH_H
