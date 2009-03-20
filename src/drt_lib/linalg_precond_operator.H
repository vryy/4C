/*!----------------------------------------------------------------------
\file  linalg_precond_operator.H

\brief A common interface for ifpack, ml and simpler preconditioners. 
       This interface allows a modification of the vector returned
       by the ApplyInverse call, which is necessary to do a solution on
       a Krylov space orthogonalized to certain (for example rigid body)
       modes.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef LINALG_PRECOND_OPERATOR_H_
#define LINALG_PRECOND_OPERATOR_H_


#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "drt_dserror.H"


namespace LINALG
{
  /*!

  A common interface for ifpack, ml and simpler preconditioners. 
  This interface allows a modification of the vector returned
  by the ApplyInverse call, which is necessary to do a solution on
  a Krylov space orthogonalized to certain (for example rigid body)
  modes.

  The linalg preconditioner interface class holds a pointer (rcp)
  to the actual preconditioner. All methods implemented to support
  the Epetra_Operator interface just call the corresponding functions
  of the actual preconditioner.

  Only the ApplyInverse method is modified.

  \author gammi
  \date 03/09
  */
class LinalgPrecondOperator : public Epetra_Operator
{
public:
  
  /*!
  \brief Standard Constructor

  */
  LinalgPrecondOperator(Teuchos::RCP<Epetra_Operator> precond,bool project);
    

  /*!
  \brief Destructor

  */
  ~LinalgPrecondOperator();

  //! @name Special atribute set methods
  void SetKernelSpace (Teuchos::RCP<Epetra_MultiVector> w)
  {
    if(project_)
    {
      if(w==Teuchos::null)
      {
        dserror("no kernel supplied for projection (but projection flag was set)");
      }
    }

    w_=w;
    return;
  }
  
  //! @}

  //! @name Atribute set methods required to support the Epetra_Operator interface
  
  virtual int SetUseTranspose (bool UseTranspose)
  {
    return(precond_->SetUseTranspose(UseTranspose));
  }
  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (pass through)
  virtual int Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const 
  {
    return(precond_->Apply(X,Y));
  }

  virtual double NormInf () const
  {
    return(precond_->NormInf());
  }
  //! @}

  //! @name Atribute access functions required to support the Epetra_Operator interface
  virtual const char * Label () const
  {
    return(precond_->Label());
  }

  virtual bool UseTranspose () const 
  {
    return(precond_->UseTranspose());
  }

  virtual bool HasNormInf () const 
  {
    return(precond_->HasNormInf());
  }

  virtual const Epetra_Comm & Comm () const
  {
    return(precond_->Comm());
  }

  virtual const Epetra_Map & OperatorDomainMap () const 
  {
    return(precond_->OperatorDomainMap());
  }

  virtual const Epetra_Map & OperatorRangeMap () const 
  {
    return(precond_->OperatorRangeMap());
  }

  //! @}

  //! @name Mathematical functions required to support the Epetra_Operator interface (modified)
  /*
    (Modified) ApplyInverse call

    This method calls ApplyInverse on the actual preconditioner. Afterwards, 
    the solution is orthogonalized against a set of vectors provided in a
    multivector, for example rigid body modes etc., using the projector 
    defined by

                                      T
                         T       w * c
                        P  = 1 - ------
                                  T
                                 w * c

    w is the vector of weights, c a vector of ones.

    The system we are solving with this procedure is not Au=f for u (since A 
    might be singular), but we are solving 

                        / T  \       T
                       | P AP | u = P f
                        \    /
                        
    for the projection of the solution Pu.
                       
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
  virtual int ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;


private:
  //! flag whether to do a projection or just pass through
  bool                               project_;

  //! the actual preconditioner
  Teuchos::RCP<Epetra_Operator>      precond_;

  //! a set of vectors defining the projector
  Teuchos::RCP<Epetra_MultiVector>   w_;
};


} // namespace LINALG






#endif // LINALG_PRECOND_OPERATOR_H_

#endif // CCADISCRET
