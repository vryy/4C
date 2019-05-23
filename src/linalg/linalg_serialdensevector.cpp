/*----------------------------------------------------------------------*/
/*!

\brief A class that wraps Epetra_SerialDenseVector with minor modifications
       in the constructor

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/

#include "linalg_serialdensevector.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector() : Epetra_SerialDenseVector(), allocatedSize_(0)
{
  SetLabel("LINALG::SerialDenseVector");
}


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(int Length, bool init)
    : Epetra_SerialDenseVector(), allocatedSize_(0)
{
  SetLabel("LINALG::SerialDenseVector");
  if (Length < 0) throw ReportError("Length = " + toString(Length) + ". Should be >= 0", -1);

  int errorcode = 0;
  if (init == true)
    errorcode = Size(Length);
  else
    errorcode = LightSize(Length);
  if (errorcode != 0) throw ReportError("LightSize returned non-zero value", errorcode);
}



/*----------------------------------------------------------------------*
 | ctor (public)                                              nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(Epetra_DataAccess CV, double* Values, int Length)
    : Epetra_SerialDenseVector(CV, Values, Length), allocatedSize_(0)
{
  SetLabel("LINALG::SerialDenseVector");
}

/*----------------------------------------------------------------------*
 | ctor (public)                                              nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(Epetra_SerialDenseVector& Source, Epetra_DataAccess CV)
    : Epetra_SerialDenseVector(CV, Source.Values(), Source.Length()), allocatedSize_(0)
{
  SetLabel("LINALG::SerialDenseVector");
}

/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(const SerialDenseVector& Source)
    : Epetra_SerialDenseVector(Source)
{
  if (CV_ == Copy) allocatedSize_ = M_;
}

/*----------------------------------------------------------------------*
 | copy-ctor (public)                                         nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(const Epetra_SerialDenseVector& Source)
    : Epetra_SerialDenseVector(Source)
{
  if (CV_ == Copy) allocatedSize_ = M_;
  SetLabel("LINALG::SerialDenseVector");
}

/*----------------------------------------------------------------------*
 | dtor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::~SerialDenseVector() {}

// << operator
// ostream& operator << (ostream& os, const LINALG::SerialDenseVector& vector)
//{
//  vector.Print(os);
//  return os;
//}



/*----------------------------------------------------------------------*
 |  size the matrix but do not init to zero  (public)        mwgee 05/07|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseVector::LightSize(int Length)
{
  if (Length < 0) return (-1);

  if (Length > allocatedSize_)
  {
    CleanupData();  // Get rid of anything that might be already allocated
    A_ = new double[Length];
    A_Copied_ = true;
    allocatedSize_ = Length;
  }
  else if (Length == 0)
    CleanupData();

  M_ = Length;
  N_ = 1;  // this is a vector, therefore ONE column
  LDA_ = M_;

  return (0);
}



/*----------------------------------------------------------------------*
 |  size the matrix but and init to zero  (public)     kronbichler 08/14|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseVector::Size(int Length)
{
  int err = LightSize(Length);
  Zero();
  return err;
}



/*----------------------------------------------------------------------*
 |  resize the matrix but do not init excess space to zero  mwgee 05/07|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseVector::LightResize(int Length)
{
  if (Length < 0) return (-1);

  if (Length == 0)
  {
    CleanupData();
    allocatedSize_ = 0;
    return 0;
  }

  double* A_tmp = 0;
  const int newsize = Length;  // ONE column

  const bool morememory = newsize > allocatedSize_;
  if (morememory)
  {
    // Allocate space for new matrix
    A_tmp = new double[newsize];
    allocatedSize_ = newsize;

    int M_tmp = EPETRA_MIN(M_, Length);
    int N_tmp = EPETRA_MIN(N_, 1);
    if (A_ != 0) CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, Length);

    CleanupData();
    A_ = A_tmp;
    A_Copied_ = true;
  }

  M_ = Length;
  N_ = 1;
  LDA_ = M_;

  return (0);
}



/*----------------------------------------------------------------------*
 |  resize the matrix and init excess space to zero    kronbichler 08/14|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseVector::Resize(int Length)
{
  const int oldsize = M_;
  int err = LightResize(Length);
  if (oldsize < Length) memset(A_ + oldsize, 0, (Length - oldsize) * sizeof(double));

  return err;
}


/*----------------------------------------------------------------------*
 |   Update vector components with scaled values of B,                  |
 |   this = ScalarThis * this + ScalarB * B         (public) bborn 08/08|
 *----------------------------------------------------------------------*/
void LINALG::SerialDenseVector::Update(const double& ScalarB, /*!< scale for input vector */
    const Epetra_SerialDenseVector& B,                        /*!< input vector */
    const double& ScalarThis                                  /*!< scale for this vector */
)
{
  Scale(ScalarThis);
  AXPY(M() * N(), ScalarB, B.A(), A());
}


/*----------------------------------------------------------------------*
 |   Set vector components to zero                                      |
 |   this = 0.0                                     (public) a.ger 11/08|
 *----------------------------------------------------------------------*/
void LINALG::SerialDenseVector::Zero() { memset(A(), 0, M() * N() * sizeof(double)); }
