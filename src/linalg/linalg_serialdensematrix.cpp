/*!----------------------------------------------------------------------
\file linalg_serialdensematrix.cpp
\brief A class that wraps Epetra_SerialDenseMatrix with minor modifications
       in the constructor

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "linalg_serialdensematrix.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(bool set_object_label) :
Epetra_SerialDenseMatrix(false)
{
  if (set_object_label) SetLabel("LINALG::SerialDenseMatrix");
}


/*----------------------------------------------------------------------*
 | ctor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(int NumRows, int NumCols,
                                             bool init, bool set_object_label) :
Epetra_SerialDenseMatrix(false)
{
  if (set_object_label) SetLabel("LINALG::SerialDenseMatrix");
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
Epetra_SerialDenseMatrix(CV,A,LDA,NumRows,NumCols,false)
{
  if (set_object_label) SetLabel("LINALG::SerialDenseMatrix");
}



/*----------------------------------------------------------------------*
 | ctor (public)                                              nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(Epetra_SerialDenseMatrix& Source,
                                             Epetra_DataAccess CV,
                                             bool set_object_label) :
Epetra_SerialDenseMatrix(CV, Source.A(), Source.LDA(), Source.M(), Source.N(), false)
{
  if (set_object_label) SetLabel("LINALG::SerialDenseMatrix");
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(const SerialDenseMatrix& Source) :
Epetra_SerialDenseMatrix(Source)
{
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                         nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::SerialDenseMatrix(const Epetra_SerialDenseMatrix& Source) :
Epetra_SerialDenseMatrix(Source)
{
  if (Source.Label()) {
    SetLabel("LINALG::SerialDenseMatrix");
  }
}


/*----------------------------------------------------------------------*
 | dtor (public)                                             mwgee 05/07|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix::~SerialDenseMatrix()
{
}

/*----------------------------------------------------------------------*
 |                                                         vlf 06/07    |
 | recursive computation of determinant of a  matrix using Sarrus rule  |
 | (uses long double to boost accuracy                                  |
 *----------------------------------------------------------------------*/
long double LINALG::SerialDenseMatrix::Det_long(void )
{
 if (N()==1)
 {
   return (*this)(0,0);
 }
 else if (N()==2)
 {
   long double out_det;
   out_det = ((long double)((*this)(0,0)) * (long double)((*this)(1,1))) - \
     ((long double)((*this)(0,1)) * (long double)((*this)(1,0)));
   return out_det;
 }
 else if (N()>2)
 {
   long double out_det=0;
   int sign=1;
   for (int i_col=0;i_col < N();i_col++)
   {

     SerialDenseMatrix temp_matrix(N()-1,N()-1);
     for (int c_col=0;c_col < i_col;c_col++)
     {
       for(int row=1;row<N();row++)
         temp_matrix(row-1,c_col)=(*this)(row,c_col);
     }
     for (int c_col=i_col+1;c_col < N();c_col++)
     {
       for(int row=1;row<N();row++)
         temp_matrix(row-1,c_col-1)=(*this)(row,c_col);
     }
     out_det = out_det + ((long double )( sign)* (long double) ((*this)(0,i_col)) \
                          * (long double) (temp_matrix.Det_long()));
     sign*=-1;
   }
   return out_det;
 }
 else return 0;
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

/*----------------------------------------------------------------------*
 |   Update matrix components with scaled values of B,                  |
 |   this = ScalarThis * this + ScalarB * B         (public) bborn 08/08|
 *----------------------------------------------------------------------*/
void LINALG::SerialDenseMatrix::Update(
  const double& ScalarB,  /*!< scale for input matrix */
  const Epetra_SerialDenseMatrix& B,  /*!< input matrix */
  const double& ScalarThis  /*!< scale for this matrix */
)
{
  Scale(ScalarThis);
  AXPY(M()*N(), ScalarB, B.A(), A());
}

/*----------------------------------------------------------------------*
 |   Set matrix components to zero                                      |
 |   this = 0.0                                     (public) bborn 08/08|
 *----------------------------------------------------------------------*/
void LINALG::SerialDenseMatrix::Zero()
{
  memset(A(), 0, M()*N()*sizeof(double));
}

