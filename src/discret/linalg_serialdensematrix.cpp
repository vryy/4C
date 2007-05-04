/*!----------------------------------------------------------------------
\file drt_node.cpp
\brief A virtual class for a node

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "linalg_serialdensematrix.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(bool set_object_label) :
Epetra_SerialDenseMatrix(set_object_label)
{
}


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(int NumRows, int NumCols,
                                             bool init, bool set_object_label) :
Epetra_SerialDenseMatrix(set_object_label)                                     
{
  if(NumRows < 0)
	throw ReportError("NumRows = " + toString(NumRows) + ". Should be >= 0", -1);
  if(NumCols < 0)
	throw ReportError("NumCols = " + toString(NumCols) + ". Should be >= 0", -1);

  int errorcode = 0;
  if (init==true)
    errorcode = Shape(NumRows, NumCols);
  else
    errorcode = LightShape(NumRows, NumCols);
  if(errorcode != 0)
    throw ReportError("Shape returned non-zero value", errorcode);
}



/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(Epetra_DataAccess CV, double* A, int LDA,
                                                   int NumRows, int NumCols,
                                                   bool set_object_label) :
Epetra_SerialDenseMatrix(CV,A,LDA,NumRows,NumCols,set_object_label)
{
}



/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(const SerialDenseMatrix& Source) :
Epetra_SerialDenseMatrix(Source)
{
}


/*----------------------------------------------------------------------*
 | dtor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::~SerialDenseMatrix()
{
}

/*----------------------------------------------------------------------*
 |  shape the matrix but do not init to zero  (public)       mwgee 05/07|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseMatrix::LightShape(int NumRows, int NumCols) 
{
  if(NumRows < 0 || NumCols < 0) return(-1);

  CleanupData(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
	const int newsize = LDA_ * N_;
	if(newsize > 0) {
		A_ = new double[newsize];
		A_Copied_ = true;
	}

  return(0);
}


/*----------------------------------------------------------------------*
 |  reshape the matrix but do not init excess space to zero  mwgee 05/07|
 *----------------------------------------------------------------------*/
int LINALG::SerialDenseMatrix::LightReshape(int NumRows, int NumCols) 
{
	if(NumRows < 0 || NumCols < 0)
		return(-1);

	double* A_tmp = 0;
	const int newsize = NumRows * NumCols;

	if(newsize > 0) {
		// Allocate space for new matrix
		A_tmp = new double[newsize];
		int M_tmp = EPETRA_MIN(M_, NumRows);
		int N_tmp = EPETRA_MIN(N_, NumCols);
		if (A_ != 0) 
			CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, NumRows); // Copy principal submatrix of A to new A
  }
  CleanupData(); // Get rid of anything that might be already allocated  
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
	if(newsize > 0) {
		A_ = A_tmp; // Set pointer to new A
		A_Copied_ = true;
	}

  return(0);
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
