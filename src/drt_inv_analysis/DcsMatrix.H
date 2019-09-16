/*----------------------------------------------------------------------*/
/*! \file
\brief Dense column storage matrix

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#ifndef DCSMATRIX_H_
#define DCSMATRIX_H_

#include "Teuchos_RCP.hpp"

#include "../drt_timestepping/timintmstep.H"


namespace INVANA
{
  /*! \class DcsMatrix
   * \brief Dense column storage matrix
   *
   * This is an interface to a matrix which is dense and distributed and
   * possibly so big that you don't want to have it stored as a
   * full matrix.
   *
   * E.g. this applies to covariance matrices in an optimization setup. These are
   * then needed to draw samples from the distribution they represent or to evaluate
   * it's value given a specific sample.
   *
   * This class particularly applies to the context of limited memory BFGS methods,
   * where the matrix, or at least its action on a vector, is approximately stored
   * with limited memory. Dense columns of these matrices can then be reproduced by
   * applying it to the unit vectors of the parameter space.
   *
   * It also provides access to itself as (not at all sparse) Epetra_CrsMatrix.
   * This function is only to be used in case the matrix is reasonably small or
   * for testing and visulization in case of small matrix size.
   *
   * todo: If you want get eigenvectors and eigenvalues. All you have to do is
   * to derive from Epetra_Operator and then it can be plugged into to the
   * AnasaziEpetra Adapter. So at least the computation of the first few vectors
   * should be computationally feasible.
   *
   * \author kehl \date 06/2016
   */
  class DcsMatrix
  {
   public:
    //! @name Construction, destruction, setup
    //@{

    //! constructor
    DcsMatrix(Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> sstore,
        Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> ystore,
        bool initscal,       /// [in] initial scaling of the hessian approximation
        bool probscal,       /// [in] scaling of functional?
        double probscalefac  /// [in] scale factor of the functional
    );

    //! destructor
    virtual ~DcsMatrix() {}

    //@}

    //! max num entries in a row
    int MaxNumEntries() const;

    //! the column map
    const Epetra_Map& RowMap() const { return *rowmap_; }

    //! Extract column copy
    int ExtractGlobalColumnCopy(
        int Col, int Length, int& NumEntries, double* Values, int* Indices) const;

    //! Extract column copy as Epetra_Vector
    int ExtractGlobalColumnCopy(int Col, Epetra_Vector& column);

    //! only for testing and presenting small examples
    Teuchos::RCP<Epetra_CrsMatrix> FillMatrix() const;

    //! extract diagonal
    int ExtractDiagonalCopy(Epetra_Vector& diagonal) const;

    //! apply matrix (inverse) to vector
    int Apply(Epetra_Vector& vector) const;

   private:
    //! precompute some quantities
    void Precompute();

    //! apply recursion to vector
    void twolooprecursion(Epetra_Vector& p) const;

    //! number of row/column entries
    int numentries_;

    //! storage of incremental solution delta
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> sstore_;

    //! storage of gradient delta
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector>> ystore_;

    //! precompute sstore'*ystore
    TIMINT::TimIntMStep<double> sty_;

    //! want initial scaling
    bool initscal_;

    //! precompute initial scaling factor
    double init_scale_fac_;

    //! was the objective function scaled
    bool probscal_;

    //! scaling factor for the problem
    double prob_scale_fac_;

    //! a column
    mutable Teuchos::RCP<Epetra_Vector> column_;

    //! latest storage index
    int past_;

    //! most recent storage index
    int future_;

    //! rowmap
    Teuchos::RCP<Epetra_Map> rowmap_;

    //! colmap
    Teuchos::RCP<Epetra_Map> colmap_;
  };
}  // namespace INVANA

#endif /* DCSMATRIX_H_ */
