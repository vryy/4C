/*
 * smoother.cpp
 *
 *  Created on: Jun 16, 2010
 *      Author: wiesner
 */



// Trilinos headers
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#include <Ifpack.h>

// BACI headers
#include "../../linalg/linalg_sparsematrix.H" // for SmootherFactory::Create

#include "smoother.H"

LINALG::Smoother::Smoother(string type, const Teuchos::RCP<Epetra_CrsMatrix>& A, ParameterList& params, FILE* outfile)
: Epetra_Operator(),
  A_(A),
  params_(params),
  type_(type)
{

}

LINALG::Smoother::~Smoother()
{

}

ostream& LINALG::Smoother::Print(std::ostream& os) const
{
  os << Label() << endl;
  return(os);
}

////////////////////////////////////////////////
// SMOOTHER IFPACK

LINALG::Smoother_Ifpack::Smoother_Ifpack(string type, const Teuchos::RCP<Epetra_CrsMatrix>& A, ParameterList& params, FILE* outfile)
: Smoother(type,A,params,outfile)
{
  Ifpack factory;
  Ifpack_Preconditioner* prec = factory.Create(type,A.get(),params_.get("smoother: ifpack overlap",0));
  prec->SetParameters(params_);
  prec->Initialize();
  prec->Compute();

  // make sure that preconditioner is computed
  if(prec->IsComputed()==false)
    dserror("Smoother_Ifpack: smoother is not computed?");

  prec_ = Teuchos::rcp(prec);
}

LINALG::Smoother_Ifpack::~Smoother_Ifpack()
{

}

int LINALG::Smoother_Ifpack::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // check if type is ILU
  if(type_ == "ILU")
  {
    // check if Y == 0.0 -> incremental ILU
    double normY = 0.0;
    Y.NormInf(&normY);
    if(normY != 0.0)
    {
      Teuchos::RCP<Epetra_MultiVector> rhs_tmp = Teuchos::rcp(new Epetra_MultiVector(X));
      Teuchos::RCP<Epetra_MultiVector> sol_tmp = Teuchos::rcp(new Epetra_MultiVector(Y));

      A_->Apply(Y,*rhs_tmp);
      rhs_tmp->Update(-1.0,X,1.0); // rhs_tmp is difference of new rhs and old rhs

      prec_->ApplyInverse(*rhs_tmp,*sol_tmp);

      Y.Update(-1.0,*sol_tmp,1.0);          // update solution sol ONLY ILU
    }
    else
    {
      // standard ILU
      prec_->ApplyInverse(X,Y);
    }

  }
  else
  {
    // it's just a jacobi or Gauss-Seidel type smoother
    prec_->ApplyInverse(X,Y);
  }

  return 0;
}

////////////////////////////////////////////////
// SMOOTHER IFPACK

Teuchos::RCP<LINALG::Smoother> LINALG::SmootherFactory::Create(const string SmootherType, const Teuchos::RCP<SparseMatrix>& A, ParameterList& params, FILE* outfile)
{
  try
  {
    LINALG::Smoother* sm = NULL;

    if(SmootherType == "ILU" ||
       SmootherType == "point relaxation" ||
       SmootherType == "Amesos")
      sm = new LINALG::Smoother_Ifpack(SmootherType,A->EpetraMatrix(),params,outfile);
    else
        throw string("SmootherType not known");

    return Teuchos::rcp(sm);
  }
  catch(string str)
  {
    cout << "Error: SmootherFactory::Create: " << str << endl;
    dserror("upps");
  }
  return Teuchos::null;
}

