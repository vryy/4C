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

#ifndef ADAPTER_SCATRA_BASE_ALGORITHM_H
#define ADAPTER_SCATRA_BASE_ALGORITHM_H

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_lib/drt_dserror.H"

// forward declarations
namespace DRT
{
  class ResultTest;
  class ScaTraResultTest;
}
namespace SCATRA
{
  class ScaTraTimIntImpl;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// general scalar transport field interface for multiphysics problems
  /*!

  \author Georg Bauer
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
        const std::string disname,             ///< scatra discretization number
        const Teuchos::ParameterList& solverparams ///< BACI solver block parameter list from dat file (default: fluid scatra solver)
    );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraBaseAlgorithm(){};

    /// access to the scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() { return scatra_; }

    /// create result test for scalar transport field
    Teuchos::RCP<DRT::ResultTest> CreateScaTraFieldTest();

    /// set fluid displacement vector due to biofilm growth
    void SetScFldGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp)
    {
      dserror("Not implemented in the base class!");
    };

    /// set structure displacement vector due to biofilm growth
    void SetScStrGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp)
    {
      dserror("Not implemented in the base class!");
    };

  private:

    /// scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

  }; // class ScaTraBaseAlgorithm

} // namespace ADAPTER

#endif

