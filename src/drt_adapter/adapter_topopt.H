/*!------------------------------------------------------------------------------------------------*
\file adapter_fluid_topopt.H

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef ADAPTER_TOPOPT_H_
#define ADAPTER_TOPOPT_H_


/// forward declarations
namespace TOPOPT
{
  class Optimizer;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /*!
  \brief general topology optimization field interface

  This adapter is used to be able to connect the optimizer with other fields.

  It shall preserve data of the optimization process for the fields which are
  computed within, f.ex. the primary and adjoint equations or the objective's
  gradient

  \author winklmaier
  \date 12/11
  */

  /// basic scalar transport solver
  class TopOptBaseAlgorithm
  {
  public:

    /// constructor
    explicit TopOptBaseAlgorithm(
        const Teuchos::ParameterList& prbdyn, ///< problem-dependent parameters
        const int disnum                   ///< scatra discretization number (default: 0)
      );

    /// virtual destructor to support polymorph destruction
    virtual ~TopOptBaseAlgorithm();

    /// access to the scalar transport field solver
    Teuchos::RCP<TOPOPT::Optimizer> Optimizer();

  private:

    /// scalar transport field solver
    Teuchos::RCP<TOPOPT::Optimizer> optimizer_;

  }; // class TopOptBaseAlgorithm

} // namespace ADAPTER

#endif /* ADAPTER_FLUID_TOPOPT_H_ */
