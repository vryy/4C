/*----------------------------------------------------------------------*/
/*! \file

\brief scalar transport field base algorithm

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_SCATRA_BASE_ALGORITHM_H
#define BACI_ADAPTER_SCATRA_BASE_ALGORITHM_H

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class ResultTest;
  class Discretization;
}  // namespace DRT

namespace SCATRA
{
  class ScaTraTimIntImpl;
}

namespace CORE::LINALG
{
  class Solver;
}

namespace ADAPTER
{
  /// general scalar transport field interface for multiphysics problems
  /*!
  \date 07/08
  */

  /// basic scalar transport solver
  class ScaTraBaseAlgorithm
  {
   public:
    /// constructor
    ScaTraBaseAlgorithm(
        const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            scatradyn,  ///< parameter list for scalar transport subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
        const std::string& disname = "scatra",       ///< name of scalar transport discretization
        const bool isale = false                     ///< ALE flag
    );

    /// virtual destructor to support polymorph destruction
    virtual ~ScaTraBaseAlgorithm() = default;

    /// initialize this class
    virtual void Init();

    /// setup this class
    virtual void Setup();

    /// access to the scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> ScaTraField() { return scatra_; }

    /// create result test for scalar transport field
    Teuchos::RCP<DRT::ResultTest> CreateScaTraFieldTest();

   private:
    /// scalar transport field solver
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool IsSetup() const { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool IsInit() const { return isinit_; };

    //! check if \ref Setup() was called
    void CheckIsSetup() const;

    //! check if \ref Init() was called
    void CheckIsInit() const;

   private:
    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class ScaTraBaseAlgorithm

}  // namespace ADAPTER

BACI_NAMESPACE_CLOSE

#endif
