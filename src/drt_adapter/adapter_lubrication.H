/*--------------------------------------------------------------------------*/
/*!
\file adapter_lubrication.H

\brief Lubrication field base algorithm

\level 3

\maintainer Alexander Seitz
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_ADAPTER_ADAPTER_LUBRICATION_H_
#define SRC_DRT_ADAPTER_ADAPTER_LUBRICATION_H_

#include <Teuchos_RCP.hpp>

// forward declarations
namespace DRT
{
  class ResultTest;
}

namespace LUBRICATION
{
  class TimIntImpl;
}

namespace IO
{
  class DiscretizationWriter;
}

namespace Teuchos
{
  class ParameterList;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER
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
    void Setup(const Teuchos::ParameterList& prbdyn,  ///< parameter list for global problem
        const Teuchos::ParameterList&
            lubricationdyn,                          ///< parameter list for Lubrication subproblem
        const Teuchos::ParameterList& solverparams,  ///< parameter list for Lubrication solver
        const std::string& disname = "lubrication",  ///< name of Lubrication discretization
        const bool isale = false                     ///< ALE flag
    );

    /// virtual destructor to support polymorph destruction
    virtual ~LubricationBaseAlgorithm(){};

    /// access to the Lubrication field solver
    Teuchos::RCP<LUBRICATION::TimIntImpl> LubricationField() { return lubrication_; }

    /// create result test for Lubrication field
    Teuchos::RCP<DRT::ResultTest> CreateLubricationFieldTest();

    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter();

   private:
    /// Lubrication field solver
    Teuchos::RCP<LUBRICATION::TimIntImpl> lubrication_;

  };  // class LubricationBaseAlgorithm

}  // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_ADAPTER_LUBRICATION_H_ */
