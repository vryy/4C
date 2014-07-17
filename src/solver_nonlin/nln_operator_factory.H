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

#ifndef NLN_OPERATOR_FACTORY_H
#define NLN_OPERATOR_FACTORY_H

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* forward declarations */

namespace NLNSOL
{
  class NlnOperator;
} // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
} // namespace Teuchos

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

      /*! \brief Returns the result of a NlnOperator inverse applied to an
       *  Epetra_MultiVector X in Y.
       *
       *  \return Integer error code, set to 0 if successful.
       *
       *  \author mayr.mt \date 06/2014
       */
      virtual Teuchos::RCP<NLNSOL::NlnOperatorBase> Create(
          const Teuchos::ParameterList& params ///< parameter list to configure the nonlinear operator
          );

      //@}

    protected:

    private:

  }; // class NlnOperatorFactory
} // namespace NLNSOL

#endif /* NLN_OPERATOR_FACTORY_H */
