/*!----------------------------------------------------------------------

\brief A common interface for ifpack, ml and simpler preconditioners.
       This interface allows a modification of the vector returned
       by the ApplyInverse call, which is necessary to do a solution on
       a Krylov space krylovized to certain (for example rigid body
       or zero pressure) modes.

\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

*----------------------------------------------------------------------*/

#ifndef LINALG_KRYLOV_PROJECTOR_H_
#define LINALG_KRYLOV_PROJECTOR_H_

#include "Teuchos_RCP.hpp"

// forward declarations
class Epetra_MultiVector;
class Epetra_Operator;
class Epetra_BlockMap;

namespace LINALG
{
  class SerialDenseMatrix;
  class SerialDenseVector;
  class SparseMatrix;

  /*!
  A class providing a Krylov projectors. Used for projected preconditioner,
  projected operator, and directly in direct solver.

  \author Keijo Nissen
  \date Feb13
  */

  class KrylovProjector
  {
   public:
    /*!
    \brief Standard Constructor, sets mode-ids and weighttype. kernel and weight
           vector as well as their inner product matrix are allocated, but still
           have to be set. Use GetNonConstKernel() and GetNonConstWeights() to
           get pointer and set the vectors and call FillComplete() tobe able to
           use projector.
    */
    KrylovProjector(const std::vector<int>
                        modeids,        //! ids of to-be-projected modes according element nullspace
        const std::string* weighttype,  //! type of weights: integration or pointvalues
        const Epetra_BlockMap* map      //! map for kernel and weight vectors
    );

    //! destructor
    ~KrylovProjector();

    //! give out Teuchos::RCP to c_ for change
    Teuchos::RCP<Epetra_MultiVector> GetNonConstKernel();

    //! give out Teuchos::RCP to w_ for change
    Teuchos::RCP<Epetra_MultiVector> GetNonConstWeights();
    // set c_ and w_ from outside
    void SetCW(Teuchos::RCP<Epetra_MultiVector> c0, Teuchos::RCP<Epetra_MultiVector> w0,
        const Epetra_BlockMap* newmap);
    void SetCW(Teuchos::RCP<Epetra_MultiVector> c0, Teuchos::RCP<Epetra_MultiVector> w0);
    //! compute (w^T c)^(-1) and completes projector for use
    void FillComplete();

    //! give out projector matrix - build it if not yet built (thus not const)
    const LINALG::SparseMatrix GetP();

    //! give out transposed projector matrix - build it if not yet built (thus not const)
    const LINALG::SparseMatrix GetPT();

    //! wrapper for applying projector to vector for iterative solver
    int ApplyP(Epetra_MultiVector& Y) const;

    //! wrapper for applying transpose of projector to vector for iterative solver
    int ApplyPT(Epetra_MultiVector& Y) const;

    //! give out projection P^T A P
    Teuchos::RCP<LINALG::SparseMatrix> Project(const LINALG::SparseMatrix& A) const;

    //! return dimension of nullspace
    int Nsdim() const { return nsdim_; }

    //! return mode-ids corresponding to element nullspace
    std::vector<int> Modes() const { return modeids_; }

    //! return type of projection weights: integration or pointvalues
    const std::string* WeightType() const { return weighttype_; }

   private:
    //! creates actual projector matrix P (or its transpose) for use in direct solver
    void CreateProjector(Teuchos::RCP<LINALG::SparseMatrix>& P,
        const Teuchos::RCP<Epetra_MultiVector>& v1, const Teuchos::RCP<Epetra_MultiVector>& v2,
        const Teuchos::RCP<LINALG::SerialDenseMatrix>& inv_v1Tv2);

    //! applies projector (or its transpose) to vector for iterative solver
    int ApplyProjector(Epetra_MultiVector& Y, const Teuchos::RCP<Epetra_MultiVector>& v1,
        const Teuchos::RCP<Epetra_MultiVector>& v2,
        const Teuchos::RCP<LINALG::SerialDenseMatrix>& inv_v1Tv2) const;

    //! multiplies Epetra_MultiVector times Epetra_SerialDenseMatrix
    Teuchos::RCP<Epetra_MultiVector> MultiplyMultiVecterDenseMatrix(
        const Teuchos::RCP<Epetra_MultiVector>& mv,
        const Teuchos::RCP<LINALG::SerialDenseMatrix>& dm) const;

    //! outer product of two Epetra_MultiVectors
    Teuchos::RCP<LINALG::SparseMatrix> MultiplyMultiVecterMultiVector(
        const Teuchos::RCP<Epetra_MultiVector>& mv1,  //! first MultiVector
        const Teuchos::RCP<Epetra_MultiVector>& mv2,  //! second MultiVector
        const int id = 1,  //! id of MultiVector form which sparsity of output matrix is estimated
        const bool fill = true  //! bool for completing matrix after computation
        ) const;

    /*
      (Modified) ApplyInverse call

      This method calls ApplyInverse on the actual preconditioner and, the
      solution is krylovized against a set of weight vectors provided in a
      multivector.

      This is done using a projector P defined by

                                      T
                                     x * w
                          P  x = x - ------ c
                                      T
                                     w * c

      w is the vector of weights, c a vector of ones (in the dofs under
      consideration) corresponding to the matrix kernel.

      The system we are solving with this procedure is not Au=b for u (since A
      might be singular), but we are solving

                          / T \         T
                         | P A | P u = P b ,
                          \   /

      for the projection of the solution Pu, i.e. in the preconditioned case


                                                            -+
             / T   \     /      -1 \          T              |
            | P * A | * |  P * M    | * xi = P  * b          |
             \     /     \         /                         |
                                                  -1         |
                                         x = P * M  * xi     |
                                                            -+


      Hence, P is always associated with the apply inverse call of the
      preconditioner (the right bracket) and always called after the call
      to ApplyInverse.


      Properties of P are:

      1) c defines the kernel of P, i.e. P projects out the matrix kernel

                            T
                           c * w
                P c = c - ------- c = c - c = 0
                            T
                           w * c

      2) The space spanned by P x is krylov to the weight vector

                         /      T      \              T
       T   /   \     T  |      x * w    |    T       x * w     T       T       T
      w * | P x | = w * | x - ------- c | = w * x - ------- * w * c = w * x - w * x = 0
           \   /        |       T       |             T
                         \     w * c   /             w * c


      This modified Apply call is for singular matrices A when c is
      a vector defining A's nullspace. The preceding projection
      operator ensures
                              |           |
                             -+-         -+-T
                    A u = A u     where u    * c =0,

      even if A*c != 0 (for numerical inaccuracies during the computation
      of A)

      See the following article for further reading:

      @article{1055401,
       author = {Bochev,, Pavel and Lehoucq,, R. B.},
       title = {On the Finite Element Solution of the Pure Neumann Problem},
       journal = {SIAM Rev.},
       volume = {47},
       number = {1},
       year = {2005},
       issn = {0036-1445},
       pages = {50--66},
       doi = {http://dx.doi.org/10.1137/S0036144503426074},
       publisher = {Society for Industrial and Applied Mathematics},
       address = {Philadelphia, PA, USA},
       }

    */

    //! flag whether inverse of (w_^T c_) was computed after w_ or c_ have been
    // given out for change with GetNonConstW() or GetNonConstC(). This is not
    // fool proof since w and c can also be changed after having called
    // FillComplete().
    bool complete_;

    //! dimension of nullspace
    int nsdim_;

    const std::vector<int> modeids_;

    const std::string* weighttype_;

    //! projector matrix - only built if necessary (e.g. for direct solvers)
    Teuchos::RCP<LINALG::SparseMatrix> P_;

    //! transposed projector matrix - only built if necessary (e.g. for direct solvers)
    Teuchos::RCP<LINALG::SparseMatrix> PT_;

    //! a set of vectors defining weighted (basis integral) vector for the projector
    Teuchos::RCP<Epetra_MultiVector> w_;

    //! a set of vectors defining the vectors of ones (in the respective components)
    //! for the matrix kernel
    Teuchos::RCP<Epetra_MultiVector> c_;

    //! inverse of product (c_^T * w_), computed once after setting c_ and w_
    Teuchos::RCP<LINALG::SerialDenseMatrix> invwTc_;
  };

}  // namespace LINALG

#endif  // LINALG_KRYLOV_PROJECTOR_H_
