/*----------------------------------------------------------------------*/
/*! \file
\brief Cholesky factorization using Epetra_SerialSpdDenseSolver

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef CHOL_FACTOR_H_
#define CHOL_FACTOR_H_

#include "chol_factor_base.H"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

class Epetra_SerialSpdDenseSolver;
class Epetra_SerialSymDenseMatrix;
class Epetra_Import;
class Epetra_Export;


namespace INVANA
{
  /*! \class CholFactor
   * \brief Cholesky factorization using Epetra_SerialSpdDenseSolver
   *
   * Takes a potentially distributed matrix, exports it to proc 0,
   * perform the factorization and redistributed the lower factor to
   * the parallel layout.
   *
   * Also provides ApplyInverse functionality using the factorization
   * to perform a forward/backword solve.
   *
   * \author kehl \date 08/2016
   */
  class CholFactor : public CholFactorBase
  {
   public:
    //! @name Construction and setup
    //@{

    //! constructor
    CholFactor(Teuchos::RCP<Epetra_CrsMatrix> A);

    //! destructor
    virtual ~CholFactor(){};

    /*! \brief Set parameters
     *
     * nothing to be set for this class
     */
    int SetParameters(Teuchos::ParameterList& List) { return 0; }

    //! initialize from given parameters
    int Initialize();

    //@}


    /*! \brief Compute the factorization
     *
     * Compute upper and lower factor using the given INVANA::DcsMatrix.
     */
    int Compute();


    //! @name Mathematical functions.
    //@{

    /*! \brief Inverse Application
     *
     * Returns the result of a forward/backward solve
     * on a Epetra_MultiVector X in Y.
     */
    int ApplyInverse(const Epetra_MultiVector& X,  ///[in] to apply the inverse to
        Epetra_MultiVector& Y                      /// [out] result
        ) const;

    //@}


    //! @name Access to the matrix and the lower factor
    //@{

    //! The lower factor
    const Epetra_CrsMatrix& H() const { return (*H_); };

    //! The matrix to be factorized.
    Teuchos::RCP<Epetra_CrsMatrix> Matrix() const { return A_; }
    //@}

    //! the associated comunictaor
    const Epetra_Comm& Comm() const { return (Comm_); };

    //! Prints basic information
    virtual std::ostream& Print(std::ostream& os) const;

   private:
    //! @name Matrices and factorization upon input and output
    //@{

    //! Reference to the matrix to be preconditioned, supposed symmetric.
    Teuchos::RCP<Epetra_CrsMatrix> A_;

    //! Contains the lower cholesky factor after Compute() has been called
    Teuchos::RCP<Epetra_CrsMatrix> H_;

    //! has Compute() been called
    bool IsComputed_;
    //@}

    //! @name Machinery for factorization and ApplyInverse
    //@{

    //! SerialDenseSolver
    Teuchos::RCP<Epetra_SerialSpdDenseSolver> solver_;

    /*! \brief The matrix and factor (after factorization)
     *
     * It has to be kept as a member to make sure solver_
     * is able to properly use it throughout the lifetime
     * if this class.
     */
    Teuchos::RCP<Epetra_SerialSymDenseMatrix> C_;

    //! Importer to bring stuff to proc0
    Teuchos::RCP<Epetra_Import> importer_;

    //! Exporter to export from poc0
    Teuchos::RCP<Epetra_Export> exporter_;

    //! Serial map
    Teuchos::RCP<Epetra_Map> serialmap_;

    //! Distributed map
    Teuchos::RCP<Epetra_Map> distributedmap_;
    //@}

    //! a commnunicator
    const Epetra_Comm& Comm_;
  };

}  // namespace INVANA
/*----------------------------------------------------------------------*/
#endif /* CHOL_FACTOR_H_ */
