/*----------------------------------------------------------------------*/
/*! \file
\brief Factory for the inverse analysis

\level 3

\maintainer Sebastian Brandstaeter
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
