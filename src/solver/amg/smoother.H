/*
 * smoother.H
 *
 *  Created on: Jun 16, 2010
 *      Author: wiesner
 */

/*!
  \file smoother.H
  \brief basis class for multigrid smoothers
 */


#ifndef SMOOTHER_H_
#define SMOOTHER_H_

// Trilinos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "../../drt_lib/drt_dserror.H"

using Teuchos::ParameterList;

namespace LINALG
{

class MapExtractor;
class BlockSparseMatrixBase;
class SparseMatrix;

/*!
 *  \brief basis class for multigrid smoothers
 *
 *  class provides a virtual interface for concrete smoother implementations.
 */
class Smoother : public virtual Epetra_Operator
{
  public:
    explicit Smoother(string type, const Teuchos::RCP<Epetra_CrsMatrix>& A, ParameterList& params, FILE* outfile = NULL);
    virtual ~Smoother();

    virtual const char* Label() const { return "Smoother"; }
    virtual const Epetra_Comm& Comm() const { return A_->Comm(); }

    virtual const Epetra_Map& OperatorDomainMap() const { return A_->OperatorDomainMap(); } // TODO fix me
    virtual const Epetra_Map& OperatorRangeMap()  const { return A_->OperatorRangeMap();  } // TODO fix me

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("ApplyInverse not implemented"); return -1;
    }

    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("Apply not implemented"); return -1;
    }

    virtual int SetUseTranspose(bool UseTranspose)
    {
      dserror("SetUseTranspose not implemented");
      return -1;
    }

    virtual bool UseTranspose() const
    {
      dserror("UseTranspose not implemented");
      return false;
    }

    virtual bool HasNormInf() const
    {
      dserror("HasNormInf not implemented");
      return false;
    }

    virtual double NormInf() const
    {
      dserror("NormInf not implemented");
      return -1.0;
    }

    virtual ostream& Print(std::ostream& os) const;

  protected:
    const Teuchos::RCP<Epetra_CrsMatrix>  A_;              ///< Input matrix
    ParameterList                         params_;
    string                                type_;
    FILE*                                 outfile_;


};

class Smoother_Ifpack : public virtual Smoother
{
  public:
    explicit Smoother_Ifpack(string type, const Teuchos::RCP<Epetra_CrsMatrix>& A, ParameterList& params, FILE* outfile = NULL);
    virtual ~Smoother_Ifpack();

    virtual const char* Label() const { return "Smoother_Ifpack"; }

    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  private:
    Teuchos::RCP<Ifpack_Preconditioner> prec_;
};

class SmootherFactory
{
  public:
    SmootherFactory() {};
    ~SmootherFactory() {};

    static Teuchos::RCP<Smoother> Create(const string SmootherType, const Teuchos::RCP<SparseMatrix>& A, ParameterList& params, FILE* outfile = NULL);
};

} /* end namespace LINALG */

#endif /* SMOOTHER_H_ */
