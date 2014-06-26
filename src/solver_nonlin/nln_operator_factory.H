/*----------------------------------------------------------------------*/
/*!
\file nln_operator_factory.H

\brief Factory for nonlinear operators

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_FACTORYH_
#define NLN_OPERATOR_FACTORYH_

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* forward declarations */

namespace NLNSOL
{
  class NlnOperator;
}

namespace Teuchos
{
  class ParameterList;
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* definition of classes */
//! NLNSOL: a nonlinear solver tool box
namespace NLNSOL
{
  /*! \brief Factory for nonlinear operators
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorFactory
  {
    public:

      //! @name Construction / Destruction
      //@{

      //! Constructor (empty)
      NlnOperatorFactory();

      //! Destructor
      virtual ~NlnOperatorFactory(){};

      //@}

      //! @name Creation on nonlinear operators
      //@{

      /*! \brief Returns the result of a NlnOperatir inverse applied to an Epetra_MultiVector X in Y.
       *
       *  \return Integer error code, set to 0 if successful.
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual Teuchos::RCP<NLNSOL::NlnOperator> Create(
          const Teuchos::ParameterList& params ///< parameter list to configure the nonlinear operator
          );

      //@}

    protected:

    private:

  };
}

#endif /* NLN_OPERATOR_FACTORYH_ */
