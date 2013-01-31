/*----------------------------------------------------------------------*/
/*!
\file linalg_serialdensevector.cpp
\brief A class that wraps Epetra_SerialDenseVector with minor modifications
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
*/
/*----------------------------------------------------------------------*/

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
 | ctor (public)                                              nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(Epetra_DataAccess CV,
                                                   double* Values,
                                                   int Length) :
Epetra_SerialDenseVector(CV,Values,Length)
{
  SetLabel("LINALG::SerialDenseVector");
}

/*----------------------------------------------------------------------*
 | ctor (public)                                              nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(Epetra_SerialDenseVector& Source,
                                             Epetra_DataAccess CV) :
Epetra_SerialDenseVector(CV,Source.Values(),Source.Length())
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
 | copy-ctor (public)                                         nis Jan13 |
 *----------------------------------------------------------------------*/
LINALG::SerialDenseVector::SerialDenseVector(const Epetra_SerialDenseVector& Source) :
Epetra_SerialDenseVector(Source)
{
  SetLabel("LINALG::SerialDenseVector");
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


/*----------------------------------------------------------------------*
 |   Update vector components with scaled values of B,                  |
 |   this = ScalarThis * this + ScalarB * B         (public) bborn 08/08|
 *----------------------------------------------------------------------*/
void LINALG::SerialDenseVector::Update(
  const double& ScalarB,  /*!< scale for input vector */
  const Epetra_SerialDenseVector& B,  /*!< input vector */
  const double& ScalarThis  /*!< scale for this vector */
)
{
  Scale(ScalarThis);
  AXPY(M()*N(), ScalarB, B.A(), A());
}


/*----------------------------------------------------------------------*
 |   Set vector components to zero                                      |
 |   this = 0.0                                     (public) a.ger 11/08|
 *----------------------------------------------------------------------*/
void LINALG::SerialDenseVector::Zero()
{
  memset(A(), 0, M()*N()*sizeof(double));
}

