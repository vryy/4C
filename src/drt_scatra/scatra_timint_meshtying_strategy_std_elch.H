/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_std_elch.H

\brief Standard solution strategy for electrochemistry problems (without meshtying)

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

\level 2
</pre>

*----------------------------------------------------------------------*/
#ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_STD_ELCH_H
#define SCATRA_TIMINT_MESHTYING_STRATEGY_STD_ELCH_H

#include "../drt_scatra/scatra_timint_elch.H"

#include "scatra_timint_meshtying_strategy_std.H"

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

  \author fang (fang@lnm.mw.tum.de)
  */

  class MeshtyingStrategyStdElch : public MeshtyingStrategyStd
  {
   public:
    //! constructor
    explicit MeshtyingStrategyStdElch(SCATRA::ScaTraTimIntElch* elchtimint);

    // destructor
    virtual ~MeshtyingStrategyStdElch() { return; };

    //! initialize system matrix
    Teuchos::RCP<LINALG::SparseOperator> InitSystemMatrix() const;
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
    void InitConvCheckStrategy();
  };  // class MeshtyingStrategyStdElch
}  // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_MESHTYING_STRATEGY_STD_ELCH_H
