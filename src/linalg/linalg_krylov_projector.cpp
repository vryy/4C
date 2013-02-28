/*!---------------------------------------------------------------------------

\file linalg_orthogonal_projector.cpp

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#include "linalg_krylov_projector.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"
#include "linalg_sparsematrix.H"
#include "linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "../drt_lib/drt_dserror.H"

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::KrylovProjector::KrylovProjector(
  const std::vector<int> modeids,
  const std::string* weighttype,
  const Epetra_BlockMap* map
  ) :
  modeids_(modeids),
  weighttype_(weighttype)
{
  nsdim_ = modeids_.size();
  c_ = Teuchos::rcp(new Epetra_MultiVector(*map,nsdim_,false));
  if (*weighttype_=="integration")
    w_ = Teuchos::rcp(new Epetra_MultiVector(*map,nsdim_,false));
  else if (*weighttype_=="pointvalues")
    w_ = c_;
  else
    dserror("No permissible weight type.");

  invwTc_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(nsdim_,nsdim_));
}// LINALG::KrylovProjector::KrylovProjector


/* --------------------------------------------------------------------
                  Give out RCP to c_ for change
   -------------------------------------------------------------------- */
Teuchos::RCP<Epetra_MultiVector> LINALG::KrylovProjector::GetNonConstKernel()
{
  // since c_ will be changed, need to call FillComplete() to recompute invwTc_
  complete_ = false;

  return c_;
}

/* --------------------------------------------------------------------
                  Give out RCP to w_ for change
   -------------------------------------------------------------------- */
Teuchos::RCP<Epetra_MultiVector> LINALG::KrylovProjector::GetNonConstWeights()
{
  if ((*weighttype_) == "pointvalues")
    dserror("For weight type 'pointvalues' weight vector equals kernel vector and can thus only be changed implicitely by changing the kernel.");

  // since w_ will be changed, need to call FillComplete() to recompute invwTc_
  complete_ = false;

  return w_;
}

/* --------------------------------------------------------------------
            Compute (w_^T c_)^(-1) and set complete flag
   -------------------------------------------------------------------- */
void LINALG::KrylovProjector::FillComplete()
{
  if(c_==Teuchos::null)
  {
    dserror("No kernel vector supplied for projection");
  }

  if(w_==Teuchos::null)
  {
    dserror("No weight vector supplied for projection");
  }

  if(c_->NumVectors() != nsdim_)
  {
    dserror("Number of kernel vectors has been changed.");
  }

  if(w_->NumVectors() != nsdim_)
  {
    dserror("Number of weight vectors has been changed.");
  }

  // loop all kernel basis vectors
  for(int mm=0;mm<nsdim_;++mm)
  {
    // loop all weight vectors
    for(int rr=0;rr<nsdim_;++rr)
    {
      /*
        Compute dot product of all different combinations of c_ and w_ and
        put result in dense matrix. In case that all <w_i,c_j>=0 for all
        i!=j, wTc_ is diagonal.

               T
              w * c
       */
      double wTc;
      ((*w_)(mm))->Dot(*((*c_)(rr)),&wTc);

      // make sure c_i and w_i must not be krylov.
      if ((rr==mm) and (abs(wTc)<1e-14))
      {
        // not sure whether c_i and w_i must not be krylov.
        // delete dserror in case you are sure what you are doing!
        dserror("weight vector w_%i must not be orthogonal to c_%i",rr,mm);
      }
      // fill matrix (w_^T * c_) - not yet inverted!
      (*invwTc_)(mm,rr) = wTc;
    }
  }

  // invert wTc-matrix (also done if it's only a scalar - check with Micheal
  // Gee before changing this)
  Epetra_SerialDenseSolver densesolver;
  densesolver.SetMatrix(*invwTc_);
  int err = densesolver.Invert();
  if (err)
    dserror("Error inverting dot-product matrix of kernels and weights for orthogonal (\"krylov\") projection.");

  complete_ = true;

  return;
}
// LINALG::KrylovProjector::FillComplete

/* --------------------------------------------------------------------
                          Destructor
   -------------------------------------------------------------------- */
LINALG::KrylovProjector::~KrylovProjector()
{
  return;
} // LINALG::KrylovProjector::~KrylovProjector

/* --------------------------------------------------------------------
                    Create projector P (for direct solvers)
   -------------------------------------------------------------------- */
Teuchos::RCP<LINALG::SparseMatrix> LINALG::KrylovProjector::CreateP() const
{
  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete().");

  /*
   *               / T   \ -1   T
   * P = I - c_ * | w_ c_ |  * w_
   *               \     /
   *
   *      `--------v--------´
   *              temp1
   */

  // initialization of P with map of w_
  Teuchos::RCP<LINALG::SparseMatrix> P = Teuchos::rcp(new LINALG::SparseMatrix(w_->Map(),81));

  //------------------------------------
  // compute temp1 (including sign "-"):
  //------------------------------------
  // create empty multivector temp1
  Epetra_MultiVector temp1(c_->Map(),nsdim_);
  // loop over all vectors of temp1
  for(int rr=0;rr<nsdim_;++rr)
  {
    // extract i-th (rr-th) vector of temp1
    Epetra_Vector temp1i(View,temp1,rr);
    // loop over all vectors of c_
    for(int mm=0;mm<nsdim_;++mm)
    {
      // scale j-th (mm-th) vector of c_ with corresponding entry of invwTc_
      // and add to i-th (rr-th) vector of temp1
      temp1i.Update(-(*invwTc_)(mm,rr),*((*c_)(mm)),1.0);
    }
  }

  //-------------------------------
  // make w_ redundant on al procs:
  //-------------------------------
  // auxiliary variables
  const int nummyrows = w_->MyLength();
  const int numvals = w_->GlobalLength();
  const double one = 1;

  // this is brutal, yet we need a Epetra_Map and not a Epetra_BlockMap
  // and we (hopefully) never ever ever use Epetra_BlockMap in BACI
  const Epetra_Map& w_map = static_cast<const Epetra_Map&>(w_->Map() );
  // fully redundant/overlapping map
  Teuchos::RCP<Epetra_Map> redundant_map =  LINALG::AllreduceEMap(w_map);
  // initialize global w without setting to 0
  Epetra_MultiVector wglob(*redundant_map,nsdim_);
  // create importer with redundant target map and distributed source map
  Epetra_Import importer(*redundant_map,w_->Map());
  // import values to global w
  wglob.Import(*w_,importer,Insert);

  //--------------------------------------------------------
  // compute P by multiplying upright temp1 with lying w_^T:
  //--------------------------------------------------------
  // loop over all proc-rows
  for(int rr=0; rr<nummyrows; ++rr)
  {
    // get global row id of current local row id
    const int grid = P->EpetraMatrix()->GRID(rr);

    // vector of all row values - prevented from growing in following loops
    std::vector<double> rowvals;
    rowvals.reserve(numvals);

    // vector of indices cooresponding to vector of rowvalues
    std::vector<int>    indices;
    indices.reserve(numvals);

    // loop over all entries of global w
    for(int mm=0; mm<numvals; ++mm)
    {
      double sum = 0;
      // loop over all kernel/weight vector
      for(int vv=0; vv<nsdim_; ++vv)
      {
        sum += (*(temp1(vv)))[rr] * (*(wglob(vv)))[mm];
      }

      // add value to vector only if non-zero
      if (sum != 0)
      {
        rowvals.push_back(sum);
        indices.push_back(w_map.GID(mm));
      }
    }
    // insert values in P
    int err = P->EpetraMatrix()->InsertGlobalValues(grid,indices.size(),rowvals.data(),indices.data());
    if (err < 0)
    {
      dserror("insertion error when trying to computekrylov projection matrix.");
    }

    // add identity matrix by adding 1 on diagonal entries
    err = P->EpetraMatrix()->InsertGlobalValues(grid,1,&one,&grid);
    if (err < 0)
    {
      err = P->EpetraMatrix()->SumIntoGlobalValues(grid,1,&one,&grid);
      if (err < 0)
      {
        dserror("insertion error when trying to computekrylov projection matrix.");
      }
    }
  }

  // call fill complete
  P->Complete();

  return P;
}


/* --------------------------------------------------------------------
                  Apply projector P(T) (for iterative solvers)
   -------------------------------------------------------------------- */
int LINALG::KrylovProjector::ApplyP(Epetra_MultiVector& Y) const
{
  /*
   *                  / T   \ -1   T
   * P(x) = x - c_ * | w_ c_ |  * w_ * x
   *                  \     /
   *                `----v-----´
   *                   invwTc_
   */

  invwTc_->SetUseTranspose(false);
  return ApplyProjector(Y,w_,c_,invwTc_);
}

int LINALG::KrylovProjector::ApplyPT(Epetra_MultiVector& Y) const
{
  /*
   *  T                / T   \ -1   T
   * P (x) = x - w_ * | c_ w_ |  * c_ * x
   *                   \     /
   *                 `----v-----´
   *                  (invwTc_)^T
   */
  invwTc_->SetUseTranspose(true);
  return ApplyProjector(Y,c_,w_,invwTc_);
}

int LINALG::KrylovProjector::ApplyProjector(
  Epetra_MultiVector& Y,
  const Teuchos::RCP<Epetra_MultiVector>& v1,
  const Teuchos::RCP<Epetra_MultiVector>& v2,
  const Teuchos::RCP<LINALG::SerialDenseMatrix>& inv_v1Tv2
  ) const
{
  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete().");

  int ierr=0;

  // if necessary, project out matrix kernel to maintain well-posedness
  // of problem
  // there is only one solution vector --- so solution
  // vector index is zero
  if(Y.NumVectors()!=1)
  {
    dserror("expecting only one solution vector during AZTEC Apply call\n");
  }

  /*
   *  (T)                /  T  \ -1    T
   * P   (x) = x - v2 * | v1 v2 |  * v1 * x
   *                     \     /
   *                                `---v---´
   *                                 =:temp1
   *                   `----------v----------´
   *                           =:temp2
   */

  // compute dot product of solution vector with all projection vectors
  // temp1(rr) = v1(rr)^T * Y
  LINALG::SerialDenseVector temp1(nsdim_);
  for(int rr=0;rr<nsdim_;++rr)
  {
    (Y(0))->Dot(*((*v1)(rr)),&(temp1(rr)));
  }

  // compute temp2 from matrix-vector-product:
  // temp2 = (v1^T v2)^(-1) * temp1
  LINALG::SerialDenseVector temp2(nsdim_);
  inv_v1Tv2->Apply(temp1,temp2);

  // loop
  for(int rr=0;rr<nsdim_;++rr)
  {
    (Y(0))->Update(-temp2(rr),*((*v2)(rr)),1.0);
  }

  return(ierr);
} // LINALG::KrylovProjector::ApplyProjector
