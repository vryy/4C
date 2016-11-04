/*----------------------------------------------------------------------*/
/*!
 * \file chol_factor.cpp
 * \brief Cholesky factorization using Epetra_SerialSpdDenseSolver
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "chol_factor.H"

#include "Epetra_Comm.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialSpdDenseSolver.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
INVANA::CholFactor::CholFactor(Teuchos::RCP<Epetra_CrsMatrix> A) :
  A_(A),
  IsComputed_(false),
  importer_(Teuchos::null),
  exporter_(Teuchos::null),
  serialmap_(Teuchos::null),
  distributedmap_(Teuchos::null),
  Comm_(A->Comm())
{
  // do nothing here
}

/*----------------------------------------------------------------------*/
int INVANA::CholFactor::Initialize()
{
  // map living on proc 0
  serialmap_ = LINALG::AllreduceEMap(A_->RowMap(),0);

  // original matrix layout
  distributedmap_ = Teuchos::rcp(new Epetra_Map(A_->RowMap()));

  // importer to bring stuff to proc 0
  importer_ = Teuchos::rcp(new Epetra_Import(*serialmap_, A_->RowMap()));

  //exporter to distribute stuff from proc 0
  exporter_ = Teuchos::rcp(new Epetra_Export(*serialmap_, A_->RowMap()));

  return 0;
}

/*----------------------------------------------------------------------*/
int INVANA::CholFactor::Compute()
{
  // number of rows
  int numelements = serialmap_->NumMyElements();

  //check for distribution
  Teuchos::RCP<Epetra_CrsMatrix> SerialMatrix = A_;
  if(A_->RowMap().DistributedGlobal())
  {
    SerialMatrix = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *serialmap_, 0));
    SerialMatrix->Import(*A_, *importer_,Insert);
    SerialMatrix->FillComplete();
  }

  // fill EpetraSymSerialDenseMatrix
  C_ = Teuchos::rcp(new Epetra_SerialSymDenseMatrix());
  C_->Shape(numelements); // this also zeroes out
  C_->SetLower();

  int numentries;
  double* vals;
  int* ind;
  for (int i=0; i<numelements; i++)
  {
    // get column
    SerialMatrix->ExtractMyRowView(i,numentries,vals,ind);
    for (int j=0; j<numentries; j++)
    {
      if (ind[j]<=i)
        (*C_)(i,ind[j]) = vals[j];
    }
  }

  // do the factorization
  solver_ = Teuchos::rcp(new Epetra_SerialSpdDenseSolver());
  solver_->SetMatrix(*C_);

  if (numelements != 0)
  {
    int err = solver_->Factor();
    if (err)
      dserror("Covariance matrix factorization failed with err %d", err);
  }
  Comm().Barrier();

  // -> CrsMatrix of the lower factor
  Epetra_SerialSymDenseMatrix* L = solver_->SymFactoredMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> HSerial = Teuchos::rcp(
      new Epetra_CrsMatrix(Copy, *serialmap_,0));
  for (int i=0; i<numelements; i++)
  {
    std::vector<int> indices(i+1);
    std::vector<double> values(i+1);

    for (int j=0; j<=i; j++)
    {
      indices[j] = j;
      values[j] = (*L)(i,j);
    }
    HSerial->InsertGlobalValues(i,values.size(),&values[0],&indices[0]);
  }
  HSerial->FillComplete();

  H_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *distributedmap_,0));
  H_->Export(*HSerial,*exporter_,Insert);
  H_->FillComplete();

  IsComputed_ = true;
  return 0;
}

/*----------------------------------------------------------------------*/
int INVANA::CholFactor::ApplyInverse(const Epetra_MultiVector& X,
                             Epetra_MultiVector& Y) const
{
  // suppose numvectors is just 1
  if (X.NumVectors() != 1 or Y.NumVectors() != 1)
    dserror("Use only 1 Vector per MultiVector.");

  // check maps mutually
  if ( not X.Map().PointSameAs(Y.Map()) )
    dserror("X and Y don't share the same maps which they"
        "should in case of Cholesky factorization.");

  // check against map of the original matrix
  if ( not X.Map().PointSameAs(*distributedmap_) )
    dserror("Maps of the vectors don't match the map of the matrix.");

  // bring them both to proc 0
  Epetra_Vector X0(*serialmap_);
  Epetra_Vector Y0(*serialmap_);
  X0.Import(*X(0),*importer_,Insert);
  Y0.Import(*Y(0),*importer_,Insert);

  int numrows = serialmap_->NumMyElements();

  // reorder to Epetra_SerialDenseMatrix
  Epetra_SerialDenseMatrix x(numrows,1);
  Epetra_SerialDenseMatrix y(numrows,1);
  for (int i=0; i<numrows; i++)
    x(i,0) = X0[i];

  int err = 0;
  // let only proc 0 solve
  if (Comm().MyPID()==0)
  {
    //solve
    solver_->SetVectors(y,x);
    err = solver_->Solve();
  }

  // only proc 0 matters
  Comm().Broadcast(&err,1,0);

  // put result in Y
  for (int i=0; i<numrows; i++)
    Y0[i] = y(i,0);

  // export back to distributed layout
  Y(0)->Export(Y0,*exporter_,Insert);

  return err;
}

/*----------------------------------------------------------------------*/
std::ostream& INVANA::CholFactor::Print(std::ostream& os) const
{
  if (!Comm().MyPID()) {
    os << std::endl;
    os << "================================================================================" << std::endl;
    os << "CholFactor: Local" << std::endl;
    os << "Global number of rows             = " << A_->NumGlobalRows() << std::endl;
    if (IsComputed_) {
      os << "Number of nonzeros of H           = " << H_->NumGlobalNonzeros() << std::endl;
      os << "nonzeros / rows                   = "
         << 1.0 * H_->NumGlobalNonzeros() / H_->NumGlobalRows() << std::endl;
    }
    os << "================================================================================" << std::endl;
    os << std::endl;
  }
  return(os);
}
