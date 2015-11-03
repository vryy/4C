/*--------------------------------------------------------------------------*/
/*!
\file adapter_reynolds.H

\brief Reynolds field base algorithm

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_ADAPTER_ADAPTER_REYNOLDS_H_
#define SRC_DRT_ADAPTER_ADAPTER_REYNOLDS_H_

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

  /// general Reynolds field interface for multiphysics problems
  /*!
  \date 11/15
  */

  /// basic Reynolds solver
  class ReynoldsBaseAlgorithm
  {
  public:

    /// constructor
    ReynoldsBaseAlgorithm(){};

    /// setup
    void Setup(
        const Teuchos::ParameterList&   prbdyn,               ///< parameter list for global problem
        const Teuchos::ParameterList&   reynoldsdyn,            ///< parameter list for Reynolds subproblem
        const Teuchos::ParameterList&   solverparams,         ///< parameter list for Reynolds solver
        const std::string&              disname = "reynolds", ///< name of Reynolds discretization
        const bool                      isale = false         ///< ALE flag
        );

    /// virtual destructor to support polymorph destruction
    virtual ~ReynoldsBaseAlgorithm(){};

    /// access to the Reynolds field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ReynoldsField() { return reynolds_; }

    /// create result test for Reynolds field
    Teuchos::RCP<DRT::ResultTest> CreateReynoldsFieldTest();

  private:

    /// Reynolds field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> reynolds_;

  }; // class ReynoldsBaseAlgorithm

} // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_ADAPTER_REYNOLDS_H_ */
