/*!----------------------------------------------------------------------
\file  linalg_krylov_projector.H

\brief A common interface for ifpack, ml and simpler preconditioners.
       This interface allows a modification of the vector returned
       by the ApplyInverse call, which is necessary to do a solution on
       a Krylov space orthogonalized to certain (for example rigid body
       or zero pressure) modes.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef LINALG_KRYLOV_PROJECTOR_H_
#define LINALG_KRYLOV_PROJECTOR_H_

#include "Teuchos_RCP.hpp"

//!@name forward decarations
class Epetra_MultiVector;
class Epetra_Operator;

namespace LINALG
{
  class SerialDenseVector;
  /*!

  A class providing Krylov projectors. Used for projected preconditioner
  and projected operator.

  \author gammi
  \date 07/09
  */
class KrylovProjector
{
public:

  /*!
  \brief Standard Constructor, sets kernel basis and weight vectors.

  The projectors P, P^T are completely defined by w and c.
  Scalar products c_i o w_i are computed once after setting c_ and w_
  and stored.

  If A is provided, A*c_i=0 is checked for all kernel basis vectors.

  \param project explicit flag if projection is desired
  \param w vector(s) of weights
  \param c kernel basis vector(s)
  \param A operator (optional) for checking Ac=0 (kernel OK)

  */
  KrylovProjector(
      const bool                       project,
      Teuchos::RCP<Epetra_MultiVector> w,
      Teuchos::RCP<Epetra_MultiVector> c,
      Teuchos::RCP<Epetra_Operator>    A = Teuchos::null
    );


  /*!
  \brief Destructor

  */
  ~KrylovProjector();

  //! @name Mathematical functions required to support the Epetra_Operator interface (pass through)
  int ApplyP(Epetra_MultiVector &Y) const;

  //! @name Mathematical functions required to support the Epetra_Operator interface (pass through)
  int ApplyPT (Epetra_MultiVector &Y) const;

  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (modified)
  /*
    (Modified) ApplyInverse call

    This method calls ApplyInverse on the actual preconditioner and, the
    solution is orthogonalized against a set of weight vectors provided in a
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

    2) The space spanned by P x is orthogonal to the weight vector

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

  //! return dimension of nullspace
  int Nsdim(){return nsdim_;}


private:
  //! flag whether to do a projection or just pass through
  bool project_;

  //! dimension of nullspace
  int  nsdim_;

  //! a set of vectors defining weighted (basis integral) vector for the projector
  Teuchos::RCP<Epetra_MultiVector> w_;

  //! a set of vectors defining the vectors of ones (in the respective components)
  //! for the matrix kernel
  Teuchos::RCP<Epetra_MultiVector> c_;

  //! values for scalar products c_i o w_i, computed once after setting c_ and w_
  Teuchos::RCP<LINALG::SerialDenseVector> cTw_;

};

} // namespace LINALG

#endif // LINALG_KRYLOV_PROJECTOR_H_

