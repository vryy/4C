/*----------------------------------------------------------------------*/
/*!
 \file scatra_timint_heterogeneous_reaction_strategy.H

 \brief Solution strategy for heterogeneous reactions

<pre>
   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
</pre>

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_H_
#define SRC_DRT_SCATRA_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_H_

#include "scatra_timint_meshtying_strategy_std.H"


namespace DRT
{
class Discretization;
}

namespace SCATRA
{

/*!
\brief Solution strategy for heterogeneous reactions

To keep the scalar transport time integrator class and derived classes as plain as possible,
several algorithmic parts have been encapsulated within separate meshtying strategy classes.
These algorithmic parts include initializing the system matrix and other relevant objects,
computing meshtying residual terms and their linearizations, and solving the resulting
linear system of equations. By introducing a hierarchy of strategies for these algorithmic
parts, a bunch of unhandy if-else selections within the time integrator classes themselves
can be circumvented. This class contains the solution strategy for heterogeneous reactions.

\author fang (fang@lnm.mw.tum.de)
*/

class HeterogeneousReactionStrategy : public MeshtyingStrategyStd
{
  public:

    //! constructor
    explicit HeterogeneousReactionStrategy(
        SCATRA::ScaTraTimIntImpl* scatratimint
        );

    //! destructor
    virtual ~HeterogeneousReactionStrategy(){return;};

    //! compute residual terms and their linearizations
    void EvaluateMeshtying();

    //! initialize meshtying objects
    void InitMeshtying();

  private:

    //! the discretization for the reaction
    Teuchos::RCP<DRT::Discretization>  discret_;

    //! private copy constructor
    HeterogeneousReactionStrategy(const HeterogeneousReactionStrategy& old);
}; // class MeshtyingStrategyStd
} // namespace SCATRA



#endif /* SRC_DRT_SCATRA_SCATRA_TIMINT_HETEROGENEOUS_REACTION_STRATEGY_H_ */
