/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_base_algorithm.H

\brief scalar transport field base algorithm

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_SCATRA_BASE_ALGORITHM_H
#define ADAPTER_SCATRA_BASE_ALGORITHM_H

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_lib/drt_globalproblem.H"

// forward declarations
namespace DRT
{
  class ResultTest;
  class ScaTraResultTest;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general scalar transport field interface for multiphysics problems
  /*!

  \author gjb
  \date 07/08
  */

  /// basic scalar transport solver
  class ScaTraBaseAlgorithm
  {
  public:

    /// constructor
    explicit ScaTraBaseAlgorithm(
        const Teuchos::ParameterList& prbdyn, ///< problem-dependent parameters
        const bool isale,                     ///< do we need the ALE formulation?
        const int disnum=0,                   ///< scatra discretization number (default: 0)
        const Teuchos::ParameterList& solverparams=DRT::Problem::Instance()->ScalarTransportFluidSolverParams() ///< BACI solver block parameter list from dat file (default: FLUID SCALAR TRANSPORT SOLVER)
      );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraBaseAlgorithm();

    /// access to the scalar transport field solver
    SCATRA::ScaTraTimIntImpl& ScaTraField();

    /// access to the reinitialization field solver
    SCATRA::ScaTraTimIntImpl& ScaTraReinitField();

    /// create result test for scalar transport field
    Teuchos::RCP<DRT::ResultTest> CreateScaTraFieldTest();

  private:

    /// scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

    /// scalar reinitialization solver for level set methods
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> reinit_;


  }; // class ScaTraBaseAlgorithm

} // namespace ADAPTER

#endif
#endif
