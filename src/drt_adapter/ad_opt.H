/*----------------------------------------------------------------------*/
/*! \file

\brief topology optimization adapter

\level 3

\maintainer Martin Kronbichler
 *------------------------------------------------------------------------------------------------*/


#ifndef AD_OPT_H_
#define AD_OPT_H_


#include <Teuchos_RCP.hpp>


// forward declarations
namespace Teuchos
{
  class ParameterList;
}

namespace TOPOPT
{
  class Optimizer;
}

namespace ADAPTER
{
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
        const Teuchos::ParameterList& prbdyn,  ///< problem-dependent parameters
        const std::string disname              ///< scatra discretization name (default: "opti")
    );

    /// virtual destructor to support polymorph destruction
    virtual ~TopOptBaseAlgorithm(){};

    /// access to the optimizer
    Teuchos::RCP<TOPOPT::Optimizer>& Optimizer() { return optimizer_; };

   private:
    /// scalar transport field solver
    Teuchos::RCP<TOPOPT::Optimizer> optimizer_;

  };  // class TopOptBaseAlgorithm

}  // namespace ADAPTER

#endif /* ADAPTER_FLUID_TOPOPT_H_ */
