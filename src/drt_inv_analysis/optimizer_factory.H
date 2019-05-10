/*----------------------------------------------------------------------*/
/*!
 \file optimizer_factory.H
 \brief factory for optimization algorithms

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef OPTIMIZER_FACTORY_H
#define OPTIMIZER_FACTORY_H

#include "optimizer_base.H"

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

#endif /*OPTIMIZER_FACTORY_H_*/
