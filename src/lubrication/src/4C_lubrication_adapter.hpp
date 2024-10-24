// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LUBRICATION_ADAPTER_HPP
#define FOUR_C_LUBRICATION_ADAPTER_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"
#include "4C_utils_result_test.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Discret
{
  class ResultTest;
}

namespace Lubrication
{
  class TimIntImpl;
}

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace Lubrication
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
    Teuchos::RCP<Lubrication::TimIntImpl> lubrication_field() { return lubrication_; }

    /// create result test for Lubrication field
    Teuchos::RCP<Core::Utils::ResultTest> create_lubrication_field_test();

    virtual Teuchos::RCP<Core::IO::DiscretizationWriter> disc_writer();

   private:
    /// Lubrication field solver
    Teuchos::RCP<Lubrication::TimIntImpl> lubrication_;

  };  // class LubricationBaseAlgorithm

}  // namespace Lubrication


FOUR_C_NAMESPACE_CLOSE

#endif
