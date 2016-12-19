/*----------------------------------------------------------------------*/
/*!
 * \file DcsMatrix.cpp
 * \brief Dense column storage matrix
 *
<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "DcsMatrix.H"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*/
INVANA::DcsMatrix::DcsMatrix(
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > sstore,
    Teuchos::RCP<TIMINT::TimIntMStep<Epetra_Vector> > ystore,
    bool initscal, bool probscal, double probscalefac) :
sty_(0,0,0.0),
initscal_(initscal),
init_scale_fac_(1.0),
probscal_(probscal),
prob_scale_fac_(probscalefac)
{
  // stupid conversion from Block map to non block map
  Teuchos::RCP<Epetra_Vector> avec = (*sstore)(0);
  int *ind;
  avec->Map().MyGlobalElementsPtr(ind);
  rowmap_ = Teuchos::rcp(new Epetra_Map(
      avec->Map().NumGlobalElements(),avec->Map().NumMyElements(),ind,0,avec->Comm()));

  colmap_ = LINALG::AllreduceEMap(*rowmap_);

  numentries_=rowmap_->NumGlobalElements();
  column_ = Teuchos::rcp(new Epetra_Vector(*rowmap_,false));

  // get size of storage
  std::pair<int,int> steps = sstore->GetSteps();
  past_ = steps.first;
  future_ = steps.second;

  sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(past_, future_, rowmap_.get(), true));
  ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(past_, future_, rowmap_.get(), true));
  sty_.Resize(past_,future_,0.0);

  for (int i=past_; i<=future_; i++)
  {
    sstore_->UpdateSteps((*sstore)[i]);
    ystore_->UpdateSteps((*ystore)[i]);
  }

  Precompute();
}

/*----------------------------------------------------------------------*/
int INVANA::DcsMatrix::MaxNumEntries() const
{
  return numentries_;
}

/*----------------------------------------------------------------------*/
int INVANA::DcsMatrix::ExtractGlobalColumnCopy(int Col, int Length, int& NumEntries, double* Values, int* Indices) const
{
  column_->PutScalar(0.0);
  column_->ReplaceGlobalValue(Col,0,1.0);

  twolooprecursion(*column_);

  // extract value and indices
  column_->ExtractCopy(Values);
  rowmap_->MyGlobalElements(Indices);
  NumEntries=numentries_;

  return 0;
}

/*----------------------------------------------------------------------*/
int INVANA::DcsMatrix::ExtractGlobalColumnCopy(int Col, Epetra_Vector& column)
{
  if ( not column_->Map().PointSameAs(column.Map()) )
    dserror("maps for extraction don't agree");

  column.PutScalar(0.0);
  column.ReplaceGlobalValue(Col,0,1.0);

  twolooprecursion(column);

  return 0;
}

/*----------------------------------------------------------------------*/
int INVANA::DcsMatrix::Apply(Epetra_Vector& vector) const
{
  if ( not column_->Map().PointSameAs(vector.Map()) )
    dserror("maps for extraction don't agree");

  twolooprecursion(vector);

  return 0;
}

/*----------------------------------------------------------------------*/
void INVANA::DcsMatrix::Precompute()
{
  // sstore'*ystore
  for (int i=past_; i<=future_; i++)
  {
    double dprod;
    (*sstore_)(i)->Dot((*ystore_)[i],&dprod);
    sty_.UpdateSteps(dprod);
  }

  // initial scaling
  double nomi=0.0;
  double denomi=0.0;
  (*sstore_)(future_)->Dot((*ystore_)[future_],&nomi);
  (*ystore_)(future_)->Dot((*ystore_)[future_],&denomi);
  init_scale_fac_ = nomi/denomi;
}

/*----------------------------------------------------------------------*/
void INVANA::DcsMatrix::twolooprecursion(Epetra_Vector& p) const
{
  std::vector<double> alpha;

  // loop steps
  for (int i=future_; i>=past_; i--)
  {
    double a=sty_[i];
    double b=0.0;
    p.Dot((*sstore_)[i],&b);
    alpha.push_back(1/a*b);

    p.Update(-1.0*alpha.back(), (*ystore_)[i],1.0 );
  }

  // initial scaling: see Nocedal, "Numerical Optimization", 2006, p. 178, formula (7.20)
  if (initscal_)
  {
    p.Scale(init_scale_fac_);
  }

  for (int i=past_; i<=future_; i++)
  {
    double a=sty_[i];
    double b=0.0;
    p.Dot((*ystore_)[i],&b);

    double alphac=alpha.back();
    alpha.pop_back();

    p.Update(alphac-1/a*b, (*sstore_)[i],1.0 );
  }

  // if the original problem was somehow scaled
  if (probscal_)
    p.Scale(prob_scale_fac_);

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix> INVANA::DcsMatrix::FillMatrix() const
{
  // a new matrix (going to be dense)
  Teuchos::RCP<Epetra_CrsMatrix> matrix = Teuchos::rcp(new
      Epetra_CrsMatrix(Copy,*rowmap_,numentries_,false));

  // loop columns
  int mylength = rowmap_->NumMyElements();
  std::vector<double> colvalues(mylength);
  std::vector<int> colindices(mylength);
  for (int col=0; col<numentries_; col++)
  {
    int col_GID = colmap_->GID(col);

    // get colum
    int colnnz;
    this->ExtractGlobalColumnCopy(col_GID,mylength,colnnz,&colvalues[0],&colindices[0]);

    for (int row=0; row<mylength; row++)
    {
      int row_GID = rowmap_->GID(row);
      double v = colvalues[row];
      matrix->InsertGlobalValues(row_GID,1,&v,&col_GID);
    }
  }
  matrix->FillComplete();

  return matrix;
}

/*----------------------------------------------------------------------*/
int INVANA::DcsMatrix::ExtractDiagonalCopy(Epetra_Vector& diagonal) const
{
  if (not rowmap_->PointSameAs(diagonal.Map()))
    dserror("Map of supplied vector doesn't fit.");

  // loop columns
  int mylength = rowmap_->NumMyElements();
  std::vector<double> colvalues(mylength);
  std::vector<int> colindices(mylength);
  for (int col=0; col<numentries_; col++)
  {
    int col_GID = colmap_->GID(col);

    // get colum
    int colnnz;
    this->ExtractGlobalColumnCopy(col_GID,mylength,colnnz,&colvalues[0],&colindices[0]);

    int lrow = rowmap_->LID(col_GID);
    if (lrow!=-1)
    {
      // (only) the proc having row col_GID has to insert
      diagonal.ReplaceGlobalValue(col_GID,0,colvalues[lrow]);
    }
  }

  return 0;
}
