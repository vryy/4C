/*----------------------------------------------------------------------*/
/*! \file
 \brief factory for optimization algorithms

\level 3

*/
/*----------------------------------------------------------------------*/
#ifndef BACI_INV_ANALYSIS_OPTIMIZER_FACTORY_H
#define BACI_INV_ANALYSIS_OPTIMIZER_FACTORY_H

#include "baci_inv_analysis_optimizer_base.H"

#include <Teuchos_RCP.hpp>


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
