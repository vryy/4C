/*----------------------------------------------------------------------*/
/*!
 * \file invana_factory.H
 * \brief factory for the inverse analysis
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_FACTORY_H
#define INVANA_FACTORY_H

#include "Teuchos_RCP.hpp"

// forward declarations
namespace Teuchos
{
class ParameterList;
}
namespace DRT
{
class Discretization;
}

namespace STR
{

namespace INVANA
{

class InvanaBase;

  /*!
  \brief A factory for the inverse analysis optimization problem
  \author kehl
 */

class InvanaFactory {

public:

  //! standard constructor
  InvanaFactory();

  //! destructor
  virtual ~InvanaFactory() {;};

  //! Return a reference to the optimization problem
  Teuchos::RCP<STR::INVANA::InvanaBase> Create(Teuchos::RCP<DRT::Discretization> discret, const Teuchos::ParameterList& invp);


};
}  // namespace INVANA
}  // namespace STR

#endif /*INVANA_FACTORY_H_*/
