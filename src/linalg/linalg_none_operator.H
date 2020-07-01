/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef LINALG_NONE_OPERATOR_H
#define LINALG_NONE_OPERATOR_H

#include <Epetra_Operator.h>
#include <Epetra_Map.h>

namespace LINALG
{
  /// identity matrix operator
  /*!
    used for do nothing preconditioning
  */
  class NoneOperator : public Epetra_Operator
  {
   public:
    /// construct from matrix
    explicit NoneOperator(const Epetra_Operator *A)
        : comm_(A->Comm()),
          domainmap_(A->OperatorDomainMap()),
          rangemap_(A->OperatorRangeMap()),
          transpose_(false)
    {
    }

    virtual int SetUseTranspose(bool t)
    {
      transpose_ = t;
      return 0;
    }

    virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
      if (&X != &Y)
      {
        return Y.Update(1., X, 0.);
      }
      return 0;
    }

    virtual int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
    {
      if (&X != &Y)
      {
        return Y.Update(1., X, 0.);
      }
      return 0;
    }

    virtual double NormInf() const { return 1.; }

    virtual const char *Label() const { return "NoneOperator"; }

    virtual bool UseTranspose() const { return transpose_; }

    virtual bool HasNormInf() const { return true; }

    virtual const Epetra_Comm &Comm() const { return comm_; }

    virtual const Epetra_Map &OperatorDomainMap() const { return domainmap_; }

    virtual const Epetra_Map &OperatorRangeMap() const { return rangemap_; }

   private:
    const Epetra_Comm &comm_;
    Epetra_Map domainmap_;
    Epetra_Map rangemap_;
    bool transpose_;
  };
}  // namespace LINALG

#endif
