/*----------------------------------------------------------------------*/
/*! \file
\brief Common cholesky factor base class

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#ifndef CHOL_FACTOR_BASE_H_
#define CHOL_FACTOR_BASE_H_

#include "Teuchos_RefCountPtr.hpp"

class Epetra_CrsMatrix;
class Epetra_MultiVector;

namespace Teuchos
{
  class ParameterList;
}

namespace INVANA
{
  /*! \class CholFactor
   * \brief Cholesky factorization interface
   *
   * \author kehl \date 08/2016
   */
  class CholFactorBase
  {
   public:
    //! @name Construction and setup
    //@{

    //! Constructor
    CholFactorBase() {}

    //! Destructor
    virtual ~CholFactorBase() {}

    //! Set parameters using a Teuchos::ParameterList object.
    virtual int SetParameters(Teuchos::ParameterList& parameterlis) = 0;

    //! initialize from given parameters
    virtual int Initialize() = 0;

    //@}


    /*! \brief Compute the factorization
     *
     * Compute upper and lower factor.
     */
    virtual int Compute() = 0;


    //! @name Mathematical functions.
    //@{

    /*! \brief Inverse Application
     *
     * Returns the result of a forward/back solve on a Epetra_MultiVector X in Y.
     */
    virtual int ApplyInverse(const Epetra_MultiVector& X,  ///[in] to apply the inverse to
        Epetra_MultiVector& Y                              /// [out] result
        ) const = 0;

    //@}


    //! @name Access to the matrix and the factorizations thereof
    //@{

    //! The lower factor
    virtual const Epetra_CrsMatrix& H() const = 0;


    //! The matrix to be factored
    virtual Teuchos::RCP<Epetra_CrsMatrix> Matrix() const = 0;

    //@}


    //! Prints basic information on iostream. This function is used by operator<<.
    virtual std::ostream& Print(std::ostream& os) const = 0;
  };

}  // namespace INVANA
/*----------------------------------------------------------------------*/
#endif /* CHOL_FACTOR_BASE_H_ */
