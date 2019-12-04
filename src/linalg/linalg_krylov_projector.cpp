/*----------------------------------------------------------------------*/
/*! \file


\brief Krylov projector

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de

*---------------------------------------------------------------------------*/

#include "linalg_krylov_projector.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"
#include "linalg_sparsematrix.H"
#include "linalg_utils_sparse_algebra_create.H"
#include "linalg_utils_densematrix_communication.H"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "Epetra_Map.h"
#include "../drt_lib/drt_dserror.H"

/* ====================================================================
    public
   ==================================================================== */

/* --------------------------------------------------------------------
                          Constructor
   -------------------------------------------------------------------- */
LINALG::KrylovProjector::KrylovProjector(
    const std::vector<int> modeids, const std::string* weighttype, const Epetra_BlockMap* map)
    : complete_(false),
      modeids_(modeids),
      weighttype_(weighttype),
      P_(Teuchos::null),
      PT_(Teuchos::null)
{
  nsdim_ = modeids_.size();
  c_ = Teuchos::rcp(new Epetra_MultiVector(*map, nsdim_, false));
  if (*weighttype_ == "integration")
    w_ = Teuchos::rcp(new Epetra_MultiVector(*map, nsdim_, false));
  else if (*weighttype_ == "pointvalues")
    w_ = c_;
  else
    dserror("No permissible weight type.");

  invwTc_ = Teuchos::rcp(new LINALG::SerialDenseMatrix(nsdim_, nsdim_));
}  // LINALG::KrylovProjector::KrylovProjector

/* --------------------------------------------------------------------
                          Destructor
   -------------------------------------------------------------------- */
LINALG::KrylovProjector::~KrylovProjector()
{
  return;
}  // LINALG::KrylovProjector::~KrylovProjector

/* --------------------------------------------------------------------
                  Give out Teuchos::RCP to c_ for change
   -------------------------------------------------------------------- */
Teuchos::RCP<Epetra_MultiVector> LINALG::KrylovProjector::GetNonConstKernel()
{
  // since c_ will be changed, need to call FillComplete() to recompute invwTc_
  complete_ = false;

  // projector matrices will change
  P_ = Teuchos::null;
  PT_ = Teuchos::null;

  return c_;
}

/* --------------------------------------------------------------------
                  Give out Teuchos::RCP to w_ for change
   -------------------------------------------------------------------- */
Teuchos::RCP<Epetra_MultiVector> LINALG::KrylovProjector::GetNonConstWeights()
{
  if ((*weighttype_) == "pointvalues")
    dserror(
        "For weight type 'pointvalues' weight vector equals kernel vector and can thus only be "
        "changed implicitely by changing the kernel.");

  // since w_ will be changed, need to call FillComplete() to recompute invwTc_
  complete_ = false;

  // projector matrices will change
  P_ = Teuchos::null;
  PT_ = Teuchos::null;

  return w_;
}

void LINALG::KrylovProjector::SetCW(Teuchos::RCP<Epetra_MultiVector> c0,
    Teuchos::RCP<Epetra_MultiVector> w0, const Epetra_BlockMap* newmap)
{
  c_ = Teuchos::null;
  w_ = Teuchos::null;

  c_ = Teuchos::rcp(new Epetra_MultiVector(*newmap, nsdim_, false));
  w_ = Teuchos::rcp(new Epetra_MultiVector(*newmap, nsdim_, false));
  *c_ = *c0;
  *w_ = *w0;
  return;
}

void LINALG::KrylovProjector::SetCW(
    Teuchos::RCP<Epetra_MultiVector> c0, Teuchos::RCP<Epetra_MultiVector> w0)
{
  *c_ = *c0;
  *w_ = *w0;
  return;
}

/* --------------------------------------------------------------------
            Compute (w_^T c_)^(-1) and set complete flag
   -------------------------------------------------------------------- */
void LINALG::KrylovProjector::FillComplete()
{
  if (c_ == Teuchos::null)
  {
    dserror("No kernel vector supplied for projection");
  }

  if (w_ == Teuchos::null)
  {
    dserror("No weight vector supplied for projection");
  }

  if (c_->NumVectors() != nsdim_)
  {
    dserror("Number of kernel vectors has been changed.");
  }

  if (w_->NumVectors() != nsdim_)
  {
    dserror("Number of weight vectors has been changed.");
  }

  // projector matrices will change
  P_ = Teuchos::null;
  PT_ = Teuchos::null;

  // loop all kernel basis vectors
  for (int mm = 0; mm < nsdim_; ++mm)
  {
    // loop all weight vectors
    for (int rr = 0; rr < nsdim_; ++rr)
    {
      /*
        Compute dot product of all different combinations of c_ and w_ and
        put result in dense matrix. In case that all <w_i,c_j>=0 for all
        i!=j, wTc_ is diagonal.

               T
              w * c
       */
      double wTc;
      ((*w_)(mm))->Dot(*((*c_)(rr)), &wTc);

      // make sure c_i and w_i must not be krylov.
      if ((rr == mm) and (abs(wTc) < 1e-14))
      {
        // not sure whether c_i and w_i must not be krylov.
        // delete dserror in case you are sure what you are doing!
        dserror("weight vector w_%i must not be orthogonal to c_%i", rr, mm);
      }
      // fill matrix (w_^T * c_) - not yet inverted!
      (*invwTc_)(mm, rr) = wTc;
    }
  }

  // invert wTc-matrix (also done if it's only a scalar - check with Micheal
  // Gee before changing this)
  Epetra_SerialDenseSolver densesolver;
  densesolver.SetMatrix(*invwTc_);
  int err = densesolver.Invert();
  if (err)
    dserror(
        "Error inverting dot-product matrix of kernels and weights for orthogonal (\"krylov\") "
        "projection.");

  complete_ = true;

  return;
}
// LINALG::KrylovProjector::FillComplete

/* --------------------------------------------------------------------
                    Create projector P(^T) (for direct solvers)
   -------------------------------------------------------------------- */
const LINALG::SparseMatrix LINALG::KrylovProjector::GetP()
{
  /*
   *               / T   \ -1   T
   * P = I - c_ * | w_ c_ |  * w_
   *               \     /
   *             `----v-----´
   *                invwTc_
   */
  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete() after changing c_ or w_.");

  if (P_ == Teuchos::null)
  {
    invwTc_->SetUseTranspose(false);
    CreateProjector(P_, w_, c_, invwTc_);
  }

  return *P_;
}

const LINALG::SparseMatrix LINALG::KrylovProjector::GetPT()
{
  /*
   *  T             / T   \ -1   T
   * P  = x - w_ * | c_ w_ |  * c_
   *
   *                \     /
   *              `----v-----´
   *               (invwTc_)^T
   */
  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete() after changing c_ or w_.");

  if (PT_ == Teuchos::null)
  {
    if ((*weighttype_) == "pointvalues")
    {
      if (P_ == Teuchos::null)
        dserror("When using type pointvalues, first get P_ than PT_. Don't ask - do it!");
      PT_ = P_;
    }
    else
    {
      invwTc_->SetUseTranspose(true);
      CreateProjector(PT_, c_, w_, invwTc_);
    }
  }

  return *PT_;
}

/* --------------------------------------------------------------------
                  Apply projector P(^T) (for iterative solvers)
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

  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete() after changing c_ or w_.");

  invwTc_->SetUseTranspose(false);
  return ApplyProjector(Y, w_, c_, invwTc_);
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

  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete() after changing c_ or w_.");

  invwTc_->SetUseTranspose(true);
  return ApplyProjector(Y, c_, w_, invwTc_);
}

/* --------------------------------------------------------------------
                  give out projection P^T A P
   -------------------------------------------------------------------- */
Teuchos::RCP<LINALG::SparseMatrix> LINALG::KrylovProjector::Project(
    const LINALG::SparseMatrix& A) const
{
  /*
   * P^T A P = A - { A c (w^T c)^-1 w^T + w (c^T w)^-1 c^T A } + w (c^T w)^-1 (c^T A c) (w^T c)^-1
   * w^T
   *                `--------v--------´   `--------v--------´ `-----------------v------------------´
   *                        mat1                mat2                            mat3
   *
   *
   *
   */

  if (!complete_)
    dserror("Krylov space projector is not complete. Call FillComplete() after changing c_ or w_.");

  // auxiliary preliminary products

  Teuchos::RCP<Epetra_MultiVector> w_invwTc = MultiplyMultiVecterDenseMatrix(w_, invwTc_);

  // here: matvec = A c_;
  Teuchos::RCP<Epetra_MultiVector> matvec =
      Teuchos::rcp(new Epetra_MultiVector(c_->Map(), nsdim_, false));
  A.EpetraMatrix()->Multiply(false, *c_, *matvec);

  // compute serial dense matrix c^T A c
  Teuchos::RCP<LINALG::SerialDenseMatrix> cTAc =
      Teuchos::rcp(new LINALG::SerialDenseMatrix(nsdim_, nsdim_, false));
  for (int i = 0; i < nsdim_; ++i)
    for (int j = 0; j < nsdim_; ++j) (*c_)(i)->Dot(*((*matvec)(j)), &((*cTAc)(i, j)));
  // std::cout << *matvec << std::endl;
  // std::cout << A << std::endl;

  // std::cout << *w_invwTc << std::endl;
  // compute and add matrices
  Teuchos::RCP<LINALG::SparseMatrix> mat1 =
      MultiplyMultiVecterMultiVector(matvec, w_invwTc, 1, false);
  {
    // put in brackets to delete mat2 immediately after being added to mat1
    // here: matvec = A^T c_;
    A.EpetraMatrix()->Multiply(true, *c_, *matvec);
    Teuchos::RCP<LINALG::SparseMatrix> mat2 =
        MultiplyMultiVecterMultiVector(w_invwTc, matvec, 2, true);
    mat1->Add(*mat2, false, 1.0, 1.0);
    mat1->Complete();
  }

  // here: matvec = w (c^T w)^-1 (c^T A c);
  matvec = MultiplyMultiVecterDenseMatrix(w_invwTc, cTAc);
  Teuchos::RCP<LINALG::SparseMatrix> mat3 =
      MultiplyMultiVecterMultiVector(matvec, w_invwTc, 1, false);
  mat3->Add(*mat1, false, -1.0, 1.0);
  mat3->Add(A, false, 1.0, 1.0);

  mat3->Complete();
  return mat3;
}

/* ====================================================================
    private methods
   ==================================================================== */

/* --------------------------------------------------------------------
                    Create projector (for direct solvers)
   -------------------------------------------------------------------- */
void LINALG::KrylovProjector::CreateProjector(Teuchos::RCP<LINALG::SparseMatrix>& P,
    const Teuchos::RCP<Epetra_MultiVector>& v1, const Teuchos::RCP<Epetra_MultiVector>& v2,
    const Teuchos::RCP<LINALG::SerialDenseMatrix>& inv_v1Tv2)
{
  /*
   *               /  T  \ -1    T
   * P = I - v2 * | v1 v2 |  * v1
   *               \     /
   *      `--------v--------´
   *              temp1
   */

  // compute temp1
  Teuchos::RCP<Epetra_MultiVector> temp1 = MultiplyMultiVecterDenseMatrix(v2, inv_v1Tv2);
  temp1->Scale(-1.0);


  // compute P by multiplying upright temp1 with lying v1^T:
  P = LINALG::KrylovProjector::MultiplyMultiVecterMultiVector(temp1, v1, 1, false);

  //--------------------------------------------------------
  // Add identity matrix
  //--------------------------------------------------------
  const int nummyrows = v1->MyLength();
  const double one = 1.0;
  // loop over all proc-rows
  for (int rr = 0; rr < nummyrows; ++rr)
  {
    // get global row id of current local row id
    const int grid = P->EpetraMatrix()->GRID(rr);

    // add identity matrix by adding 1 on diagonal entries
    int err = P->EpetraMatrix()->InsertGlobalValues(grid, 1, &one, &grid);
    if (err < 0)
    {
      err = P->EpetraMatrix()->SumIntoGlobalValues(grid, 1, &one, &grid);
      if (err < 0)
      {
        dserror("insertion error when trying to computekrylov projection matrix.");
      }
    }
  }

  // call fill complete
  P->Complete();

  return;
}


/* --------------------------------------------------------------------
                  Apply projector P(T) (for iterative solvers)
   -------------------------------------------------------------------- */
int LINALG::KrylovProjector::ApplyProjector(Epetra_MultiVector& Y,
    const Teuchos::RCP<Epetra_MultiVector>& v1, const Teuchos::RCP<Epetra_MultiVector>& v2,
    const Teuchos::RCP<LINALG::SerialDenseMatrix>& inv_v1Tv2) const
{
  if (!complete_) dserror("Krylov space projector is not complete. Call FillComplete().");

  int ierr = 0;

  // if necessary, project out matrix kernel to maintain well-posedness
  // of problem
  // there is only one solution vector --- so solution
  // vector index is zero
  if (Y.NumVectors() != 1)
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
  for (int rr = 0; rr < nsdim_; ++rr)
  {
    (Y(0))->Dot(*((*v1)(rr)), &(temp1(rr)));
  }

  // compute temp2 from matrix-vector-product:
  // temp2 = (v1^T v2)^(-1) * temp1
  LINALG::SerialDenseVector temp2(nsdim_);
  inv_v1Tv2->Apply(temp1, temp2);

  // loop
  for (int rr = 0; rr < nsdim_; ++rr)
  {
    (Y(0))->Update(-temp2(rr), *((*v2)(rr)), 1.0);
  }

  return (ierr);
}  // LINALG::KrylovProjector::ApplyProjector



/* --------------------------------------------------------------------
         multiplies Epetra_MultiVector times Epetra_SerialDenseMatrix
   -------------------------------------------------------------------- */
Teuchos::RCP<Epetra_MultiVector> LINALG::KrylovProjector::MultiplyMultiVecterDenseMatrix(
    const Teuchos::RCP<Epetra_MultiVector>& mv,
    const Teuchos::RCP<LINALG::SerialDenseMatrix>& dm) const
{
  if (mv == Teuchos::null or dm == Teuchos::null)
    dserror("Either the multivector or the densematrix point to Teuchos::null");

  // create empty multivector mvout
  Teuchos::RCP<Epetra_MultiVector> mvout = Teuchos::rcp(new Epetra_MultiVector(mv->Map(), nsdim_));

  // use depending on whether dm is set Transpose or not:
  if (dm->UseTranspose())
    // loop over all vectors of mvout
    for (int rr = 0; rr < nsdim_; ++rr)
    {
      // extract i-th (rr-th) vector of mvout
      Epetra_Vector mvouti(::View, *mvout, rr);
      // loop over all vectors of mv
      for (int mm = 0; mm < nsdim_; ++mm)
      {
        // scale j-th (mm-th) vector of mv with corresponding entry of dm
        // and add to i-th (rr-th) vector of mvout
        mvouti.Update((*dm)(rr, mm), *((*mv)(mm)), 1.0);
      }
    }
  else
  {
    // loop over all vectors of mvout
    for (int rr = 0; rr < nsdim_; ++rr)
    {
      // extract i-th (rr-th) vector of mvout
      Epetra_Vector mvouti(::View, *mvout, rr);
      // loop over all vectors of mv
      for (int mm = 0; mm < nsdim_; ++mm)
      {
        // scale j-th (mm-th) vector of mv with corresponding entry of dm
        // and add to i-th (rr-th) vector of mvout
        mvouti.Update((*dm)(mm, rr), *((*mv)(mm)), 1.0);
      }
    }
  }

  return mvout;
}

/* --------------------------------------------------------------------
                  outer product of two Epetra_MultiVectors
   -------------------------------------------------------------------- */
Teuchos::RCP<LINALG::SparseMatrix> LINALG::KrylovProjector::MultiplyMultiVecterMultiVector(
    const Teuchos::RCP<Epetra_MultiVector>& mv1, const Teuchos::RCP<Epetra_MultiVector>& mv2,
    const int id, const bool fill) const
{
  if (mv1 == Teuchos::null or mv2 == Teuchos::null)
    dserror("At least one multi-vector points to Teuchos::null.");

  // compute information about density of P^T A P
  Teuchos::RCP<Epetra_MultiVector> temp = Teuchos::null;
  if (id == 1)
    temp = mv1;
  else if (id == 2)
    temp = mv2;
  else
    dserror("id must be 1 or 2");

  Epetra_Vector prod(*((*temp)(0)));
  for (int i = 1; i < nsdim_; ++i) prod.Multiply(1.0, *((*temp)(i)), prod, 1.0);
  int numnonzero = 0;
  for (int i = 0; i < prod.MyLength(); ++i)
    if (prod[i] != 0.0) numnonzero++;

  int glob_numnonzero = 0;
  prod.Comm().SumAll(&numnonzero, &glob_numnonzero, 1);

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> mv1map = Teuchos::rcp(new Epetra_Map(mv1->Map().NumGlobalElements(),
      mv1->Map().NumMyElements(), mv1->Map().MyGlobalElements(), 0, mv1->Map().Comm()));
  // initialization of mat with map of mv1
  Teuchos::RCP<LINALG::SparseMatrix> mat =
      Teuchos::rcp(new LINALG::SparseMatrix(*mv1map, glob_numnonzero, false));

  //-------------------------------
  // make mv2 redundant on all procs:
  //-------------------------------
  // auxiliary variables
  const int nummyrows = mv1->MyLength();
  const int numvals = mv2->GlobalLength();

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> mv2map = Teuchos::rcp(new Epetra_Map(mv2->Map().NumGlobalElements(),
      mv2->Map().NumMyElements(), mv2->Map().MyGlobalElements(), 0, mv2->Map().Comm()));

  // fully redundant/overlapping map
  Teuchos::RCP<Epetra_Map> redundant_map = LINALG::AllreduceEMap(*mv2map);
  // initialize global mv2 without setting to 0
  Epetra_MultiVector mv2glob(*redundant_map, nsdim_);
  // create importer with redundant target map and distributed source map
  Epetra_Import importer(*redundant_map, mv2->Map());
  // import values to global mv2
  mv2glob.Import(*mv2, importer, Insert);

  //--------------------------------------------------------
  // compute mat by multiplying upright mv1 with lying mv2^T:
  //--------------------------------------------------------
  // loop over all proc-rows
  for (int rr = 0; rr < nummyrows; ++rr)
  {
    // get global row id of current local row id
    const int grid = mat->EpetraMatrix()->GRID(rr);

    // vector of all row values - prevented from growing in following loops
    std::vector<double> rowvals;
    rowvals.reserve(numvals);

    // vector of indices corresponding to vector of row values
    std::vector<int> indices;
    indices.reserve(numvals);

    // loop over all entries of global w
    for (int mm = 0; mm < numvals; ++mm)
    {
      double sum = 0.0;
      // loop over all kernel/weight vector
      for (int vv = 0; vv < nsdim_; ++vv)
      {
        sum += (*((*mv1)(vv)))[rr] * (*(mv2glob(vv)))[mm];
      }

      // add value to vector only if non-zero
      if (sum != 0.0)
      {
        rowvals.push_back(sum);
        indices.push_back(redundant_map->GID(mm));
      }
    }

    // insert values in mat
    int err = mat->EpetraMatrix()->InsertGlobalValues(
        grid, indices.size(), rowvals.data(), indices.data());
    if (err < 0)
    {
      dserror(
          "insertion error when trying to compute krylov projection matrix (error code: %i).", err);
    }
  }

  // call fill complete
  if (fill) mat->Complete();

  return mat;
}
