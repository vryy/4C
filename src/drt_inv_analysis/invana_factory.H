/*----------------------------------------------------------------------*/
/*!
\file invana_factory.H

\brief Factory for the inverse analysis

<pre>
\level 3
\maintainer Sebastian Brandstaeter
            brandstaeter@lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef INVANA_FACTORY_H
#define INVANA_FACTORY_H

#include <Teuchos_RCP.hpp>

// forward declarations
namespace Teuchos
{
  class ParameterList;
}
namespace DRT
{
  class Discretization;
}

namespace INVANA
{
  class InvanaBase;

  /*!
  \brief A factory for the inverse analysis optimization problem
  \author kehl
 */

  class InvanaFactory
  {
   public:
    //! standard constructor
    InvanaFactory();

    //! destructor
    virtual ~InvanaFactory() { ; };

    //! Return a reference to the optimization problem
    Teuchos::RCP<INVANA::InvanaBase> Create(const Teuchos::ParameterList& invp);
  };
}  // namespace INVANA

#endif /*INVANA_FACTORY_H_*/
