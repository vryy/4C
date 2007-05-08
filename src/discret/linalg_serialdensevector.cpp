/*!----------------------------------------------------------------------
\file linalg_serialdensevector.cpp
\brief A class that wraps Epetra_SerialDenseVector with minor modifications
       in the constructor

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "linalg_serialdensevector.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector() :
Epetra_SerialDenseVector()
{
  SetLabel("LINALG::SerialDenseVector");
}


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(int Length, bool init) :
Epetra_SerialDenseVector()                                     
{
  SetLabel("LINALG::SerialDenseVector");
  if(Length < 0)
  throw ReportError("Length = " + toString(Length) + ". Should be >= 0", -1);

  int errorcode = 0;
  if (init==true)
    errorcode = Shape(Length, 1);
  else
    errorcode = LightSize(Length);
  if(errorcode != 0)
    throw ReportError("LightSize returned non-zero value", errorcode);
}



/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(Epetra_DataAccess CV,
                                                   double* Values,
                                                   int Length) :
Epetra_SerialDenseVector(CV,Values,Length)
{
  SetLabel("LINALG::SerialDenseVector");
}



/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(const SerialDenseVector& Source) :
Epetra_SerialDenseVector(Source)
{
}


/*----------------------------------------------------------------------*
 | dtor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::~SerialDenseVector()
{
}

// << operator
//ostream& operator << (ostream& os, const LINALG::SerialDenseVector& vector)
//{
//  vector.Print(os);
//  return os;
//}

/*----------------------------------------------------------------------*
 |  size the matrix but do not init to zero  (public)       mwgee 05/07|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseVector::LightSize(int Length) 
{
  if(Length < 0) return(-1);

  CleanupData(); // Get rid of anything that might be already allocated
  M_ = Length;
  N_ = 1;        // this is a vector, therefore ONE column
  LDA_ = M_;
	const int newsize = LDA_ * N_;
	if(newsize > 0) {
		A_ = new double[newsize];
		A_Copied_ = true;
	}

  return(0);
}


/*----------------------------------------------------------------------*
 |  resize the matrix but do not init excess space to zero  mwgee 05/07|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseVector::LightResize(int Length) 
{
	if(Length < 0)
		return(-1);

	double* A_tmp = 0;
	const int newsize = Length;    // ONE column

	if(newsize > 0) {
		// Allocate space for new matrix
		A_tmp = new double[newsize];
		int M_tmp = EPETRA_MIN(M_, Length);
		int N_tmp = EPETRA_MIN(N_, 1);
		if (A_ != 0) 
			CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, Length); // Copy principal submatrix of A to new A
  }
  CleanupData(); // Get rid of anything that might be already allocated  
  M_ = Length;
  N_ = 1;
  LDA_ = M_;
	if(newsize > 0) {
		A_ = A_tmp; // Set pointer to new A
		A_Copied_ = true;
	}

  return(0);
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
