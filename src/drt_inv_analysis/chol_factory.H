/*----------------------------------------------------------------------*/
/*!
\file chol_factory.H

\brief Factory for the cholesky factorization

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef CHOL_FACTORY_H
#define CHOL_FACTORY_H

#include "Teuchos_RCP.hpp"

class Epetra_CrsMatrix;

namespace Teuchos{
  class ParameterList;
}

namespace INVANA
{

class CholFactorBase;

/*!
 * \brief Factory for the Cholesky factorization
 * \author kehl
 */
class CholFactory {

public:

  //! standard constructor
  CholFactory(){}

  //! destructor
  virtual ~CholFactory() {}

  //! Return the factorization base pointer
  Teuchos::RCP<CholFactorBase> Create(Teuchos::RCP<Epetra_CrsMatrix> A,
      Teuchos::ParameterList& p
      );
};
}  // namespace INVANA

#endif /*CHOL_FACTORY_H_*/
