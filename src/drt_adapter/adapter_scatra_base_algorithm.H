/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_base_algorithm.H

\brief scalar transport field base algorithm

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_SCATRA_BASE_ALGORITHM_H
#define ADAPTER_SCATRA_BASE_ALGORITHM_H

#include <Teuchos_RCP.hpp>

// forward declarations
namespace DRT
{
  class ResultTest;
}

namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace Teuchos
{
  class ParameterList;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general scalar transport field interface for multiphysics problems
  /*!
  \date 07/08
  */

  /// basic scalar transport solver
  class ScaTraBaseAlgorithm
  {
  public:

    /// constructor
    explicit ScaTraBaseAlgorithm(
        const Teuchos::ParameterList&   prbdyn,               ///< parameter list for global problem
        const Teuchos::ParameterList&   scatradyn,            ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList&   solverparams,         ///< parameter list for scalar transport solver
        const std::string&              disname = "scatra",   ///< name of scalar transport discretization
        const bool                      isale = false         ///< ALE flag
        );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraBaseAlgorithm(){};

    /// access to the scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() { return scatra_; }

    /// create result test for scalar transport field
    Teuchos::RCP<DRT::ResultTest> CreateScaTraFieldTest();

  private:

    /// scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

  }; // class ScaTraBaseAlgorithm

} // namespace ADAPTER

#endif

