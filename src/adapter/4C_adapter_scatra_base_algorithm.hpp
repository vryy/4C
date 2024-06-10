/*----------------------------------------------------------------------*/
/*! \file

\brief scalar transport field base algorithm

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_SCATRA_BASE_ALGORITHM_HPP
#define FOUR_C_ADAPTER_SCATRA_BASE_ALGORITHM_HPP

#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace ScaTra
{
  class ScaTraTimIntImpl;
}

namespace Core::LinAlg
{
  class Solver;
}

namespace Adapter
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
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> ScaTraField() { return scatra_; }

    /// create result test for scalar transport field
    Teuchos::RCP<Core::UTILS::ResultTest> create_sca_tra_field_test();

   private:
    /// scalar transport field solver
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra_;

   private:
    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

   protected:
    //! returns true if Setup() was called and is still valid
    bool is_setup() const { return issetup_; };

    //! returns true if Init(..) was called and is still valid
    bool is_init() const { return isinit_; };

    //! check if \ref Setup() was called
    void check_is_setup() const;

    //! check if \ref Init() was called
    void check_is_init() const;

   private:
    //! set flag true after setup or false if setup became invalid
    void set_is_setup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void set_is_init(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class ScaTraBaseAlgorithm

}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
