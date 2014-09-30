/*----------------------------------------------------------------------*/
/*!
 \file optimizer_factory.H
 \brief factory for optimization algorithms

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
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
namespace STR
{

namespace INVANA
{
class InvanaBase;

  /*!
  \brief A factory for the inverse analysis optimzation algorithm

  \author kehl
 */

class OptimizerFactory {

public:

  //! standard constructor
  OptimizerFactory();

  //! destructor
  virtual ~OptimizerFactory() {};

  //! Return a reference to the optimization problem
  Teuchos::RCP<STR::INVANA::OptimizerBase> Create(const Teuchos::ParameterList& invp);


};
}  // namespace INVANA
}  // namespace STR

#endif /*OPTIMIZER_FACTORY_H_*/
