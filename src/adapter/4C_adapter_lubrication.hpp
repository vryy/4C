/*--------------------------------------------------------------------------*/
/*! \file

\brief Lubrication field base algorithm

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_LUBRICATION_HPP
#define FOUR_C_ADAPTER_LUBRICATION_HPP

#include "4C_config.hpp"

#include "4C_utils_result_test.hpp"

#include <Teuchos_RCP.hpp>

namespace Teuchos
{
  class ParameterList;
}

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class ResultTest;
}

namespace LUBRICATION
{
  class TimIntImpl;
}

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Adapter
{
  /// general Lubrication field interface for multiphysics problems
  /*!
  \date 11/15
  */

  /// basic Lubrication solver
  class LubricationBaseAlgorithm
  {
   public:
    /// constructor
    LubricationBaseAlgorithm(){};

    /// setup
    void setup(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            lubricationdyn,                          ///< parameter list for Lubrication subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for Lubrication solver
        const std::string& disname = "lubrication",  ///< name of Lubrication discretization
        const bool isale = false                     ///< ALE flag
    );

    /// virtual destructor to support polymorph destruction
    virtual ~LubricationBaseAlgorithm() = default;

    /// access to the Lubrication field solver
    Teuchos::RCP<LUBRICATION::TimIntImpl> LubricationField() { return lubrication_; }

    /// create result test for Lubrication field
    Teuchos::RCP<Core::UTILS::ResultTest> create_lubrication_field_test();

    virtual Teuchos::RCP<Core::IO::DiscretizationWriter> DiscWriter();

   private:
    /// Lubrication field solver
    Teuchos::RCP<LUBRICATION::TimIntImpl> lubrication_;

  };  // class LubricationBaseAlgorithm

}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
