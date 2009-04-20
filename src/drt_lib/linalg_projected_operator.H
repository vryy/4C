/*!----------------------------------------------------------------------
\file  linalg_projected_operator.H

\brief 

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef LINALG_PROJECTED_OPERATOR_H_
#define LINALG_PROJECTED_OPERATOR_H_

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "drt_dserror.H"


namespace LINALG
{
  /*! 
  \brief A slightly modified Epetra_Operator class allowing the 
         application of a 'rank 1' projector P before applying the 
         underlying operator

  Only the Apply method is modified. All other calls a 'pass-through'.
  The modification of the apply call is such that for singular matrices
  the matrix kernel can be projected out before the apply call in order
  to become independent of numerical errors, i.e. in the iterative 
  solution process for a singular matrix A this modified apply call
  appears in the following equation:  
  
                   T                    /             \ 
                  P * InversePrecond * | A * P * u -f  | = 0
                                        \             / 
                                         |   |
                                         +---+
                                      included in
                                     the  modified 
                                       Apply call

  \author gammi
  \date 04/09
  */
class LinalgProjectedOperator : public Epetra_Operator
{
public:
  
  /*!
  \brief Standard Constructor

  */
  LinalgProjectedOperator(Teuchos::RCP<Epetra_Operator> A,bool project);
    

  /*!
  \brief Destructor

  */
  ~LinalgProjectedOperator();

  //! @name Special atribute set methods
  /*!

  \brief Set kernel space and weight vectors.
  
  The projector P is completely defined by w and c.

  \param w vector(s) of weights 
  \param c kernel basis vector(s)

  */
  void SetKernelSpace(
    Teuchos::RCP<Epetra_MultiVector> w,
    Teuchos::RCP<Epetra_MultiVector> c
    )
  {
    if(project_)
    {
      if(w==Teuchos::null || c==Teuchos::null)
      {
        dserror("no kernel/weight supplied for projection (but projection flag was set)");
      }
    }

    w_=w;
    c_=c;
    return;
  }
  //! @}

  //! @name Atribute set methods required to support the Epetra_Operator interface
  
  virtual int SetUseTranspose (bool UseTranspose)
  {
    return(A_->SetUseTranspose(UseTranspose));
  }
  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (modified)
  /*!
    \brief (Modified) Apply call

    Given the underlying Epetra_Operator A, we apply AP using 

                                    T   
                                   x * w
                        P x = x - ------- c
                                    T
                                   w * c

    instead of A if projetion_ flag is set. If projetion_ is false,
    this is just the standard Apply call.

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
                       
    See the following article for further reading

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
    
     and the documentation to LinalgPrecondOperator as well as the source
     to LINALG::Solver::Solve_aztec.

  */
  //! @}
  virtual int Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

  virtual double NormInf () const
  {
    return(A_->NormInf());
  }
  //! @}

  //! @name Atribute access functions required to support the Epetra_Operator interface
  virtual const char * Label () const
  {
    return(A_->Label());
  }

  virtual bool UseTranspose () const 
  {
    return(A_->UseTranspose());
  }

  virtual bool HasNormInf () const 
  {
    return(A_->HasNormInf());
  }

  virtual const Epetra_Comm & Comm () const
  {
    return(A_->Comm());
  }

  virtual const Epetra_Map & OperatorDomainMap () const 
  {
    return(A_->OperatorDomainMap());
  }

  virtual const Epetra_Map & OperatorRangeMap () const 
  {
    return(A_->OperatorRangeMap());
  }

  virtual int ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    return(A_->ApplyInverse(X,Y));
  }

  Teuchos::RCP<Epetra_Operator> UnprojectedOperator ()
  {
    return(A_);
  }

private:
  //! flag whether to do a projection or just pass through
  bool                               project_;

  //! the actual unprojected operator
  Teuchos::RCP<Epetra_Operator>      A_;

  //! a set of vectors defining weighted basis mean vector for the projector
  Teuchos::RCP<Epetra_MultiVector>   w_;

  //! a set of vectors defining the vectors of ones (in the respective components) 
  //! for the matrix kernel
  Teuchos::RCP<Epetra_MultiVector>   c_;

};

} // namespace LINALG

#endif // LINALG_PROJECTED_OPERATOR_H_

#endif // CCADISCRET
