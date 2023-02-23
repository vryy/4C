/*----------------------------------------------------------------------*/
/*! \file
 \brief factory for optimization algorithms

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef INV_ANALYSIS_OPTIMIZER_FACTORY_H
#define INV_ANALYSIS_OPTIMIZER_FACTORY_H

#include "inv_analysis_optimizer_base.H"

#include "Teuchos_RCP.hpp"


// forward declarations
namespace Teuchos
{
  class ParameterList;
}

namespace INVANA
{
  class InvanaBase;

  /*!
  \brief A factory for the inverse analysis optimzation algorithm

  \author kehl \date 08/2014
 */

  class OptimizerFactory
  {
   public:
    //! standard constructor
    OptimizerFactory();

    //! destructor
    virtual ~OptimizerFactory(){};

    //! Return a reference to the optimization problem
    Teuchos::RCP<INVANA::OptimizerBase> Create(const Teuchos::ParameterList& invp);
  };
}  // namespace INVANA

#endif  // INV_ANALYSIS_OPTIMIZER_FACTORY_H
