/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of eigenvalue methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_densematrix_eigen.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  compute all eigenvalues of a real symmetric matrix A        lw 04/08|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigenValues(
    Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& L, const bool postproc)
{
  LINALG::SymmetricEigen(A, L, 'N', postproc);
}

/*----------------------------------------------------------------------*
 |  compute all eigenvalues and eigenvectors of a real symmetric        |
 |  matrix A (eigenvectors are stored in A, i.e. original matrix        |
 |  is destroyed!!!)                                            lw 04/08|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigenProblem(
    Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& L, const bool postproc)
{
  LINALG::SymmetricEigen(A, L, 'V', postproc);
}

/*----------------------------------------------------------------------*
 |  compute all eigenvalues and, optionally,                            |
 |  eigenvectors of a real symmetric matrix A  (public)        maf 06/07|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigen(
    Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& L, const char jobz, const bool postproc)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != L.Length()) dserror("Dimension of eigenvalues does not match");

  double* a = A.A();
  double* w = L.A();
  const char uplo = {'U'};
  const int lda = A.LDA();
  const int dim = A.M();

  int liwork = 0;
  if (dim == 1)
    liwork = 1;
  else
  {
    if (jobz == 'N')
      liwork = 1;
    else if (jobz == 'V')
      liwork = 3 + 5 * dim;
  }
  std::vector<int> iwork(liwork);

  int lwork = 0;
  if (dim == 1)
    lwork = 1;
  else
  {
    if (jobz == 'N')
      lwork = 2 * dim + 1;
    else if (jobz == 'V')
      lwork = 2 * dim * dim + 6 * dim + 1;
  }
  std::vector<double> work(lwork);
  int info = 0;

  Epetra_LAPACK lapack;

  lapack.SYEVD(jobz, uplo, dim, a, lda, w, &(work[0]), lwork, &(iwork[0]), liwork, &info);

  if (!postproc)
  {
    if (info > 0) dserror("Lapack algorithm syevd failed");
    if (info < 0) dserror("Illegal value in Lapack syevd call");
  }
  // if we only calculate eigenvalues/eigenvectors for postprocessing,
  // a warning might be sufficient
  else
  {
    if (info > 0)
      std::cout
          << "Lapack algorithm syevd failed: " << info
          << " off-diagonal elements of intermediate tridiagonal form did not converge to zero"
          << std::endl;
    if (info < 0) std::cout << "Illegal value in Lapack syevd call" << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  compute all eigenvalues for the generalized eigenvalue problem
 |  Ax =  lambda Bx via QZ-algorithm (B is singular) and returns the
 |  maximum eigenvalue                              shahmiri  05/13
 *----------------------------------------------------------------------*/
double LINALG::GeneralizedEigen(Epetra_SerialDenseMatrix& A, Epetra_SerialDenseMatrix& B)
{
  Epetra_SerialDenseMatrix tmpA(A);
  Epetra_SerialDenseMatrix tmpB(B);

  //--------------------------------------------------------------------
  // STEP 1:
  // Transform the matrix B to upper triangular matrix with the help of
  // QR-factorization
  //--------------------------------------------------------------------

  int N = tmpA.M();
  double* a = tmpA.A();
  double* b = tmpB.A();
  int LDA = tmpA.LDA();
  int LDB = tmpB.LDA();

  // the order of permutation matrix
  int jpvt[N];
  // factor uses for calculating orthogonal matrix Q
  double tau[N];
  for (int i = 0; i < N; ++i)
  {
    jpvt[i] = 0;
    tau[i] = 0.;
  }

  int lwork1 = 0;
  if (N == 1)
    lwork1 = 1;
  else
  {
    lwork1 = N * 3 + 1;
  }
  double work1[lwork1];
  int info;
  dgeqp3(&N, &N, b, &LDB, jpvt, tau, work1, &lwork1, &info);

  if (info < 0)
    std::cout << "Lapack algorithm dgeqp3: The " << info << "-th argument had an illegal value"
              << std::endl;

  // calculate the matrix Q from multiplying householder transformations
  // Q = H(1)*H(2) ... H(k)
  // H(i) = I - tau-v*v**T
  // v is a vector with v(1:i-1) = 0 and v(i+1:m) is stored on exit in B(i+1:m,i)

  // Q is initialized as an unit matrix
  Epetra_SerialDenseMatrix Q_new(true);
  Q_new.Shape(N, N);
  for (int i = 0; i < N; ++i) Q_new(i, i) = 1.0;

  for (int i = 0; i < N; ++i)
  {
    Epetra_SerialDenseVector v;
    v.Shape(N, 1);
    v(i, 0) = 1.;
    for (int j = i + 1; j < N; ++j) v(j, 0) = tmpB(j, i);

    Epetra_SerialDenseMatrix H;
    H.Shape(N, N);

    H.Multiply('N', 'T', tau[i], v, v, 0.);
    H.Scale(-1.);
    for (int k = 0; k < N; ++k) H(k, k) = 1. + H(k, k);

    Epetra_SerialDenseMatrix Q_help;
    Q_help.Shape(N, N);
    Q_new.Apply(H, Q_help);
    Q_new = Q_help;
  }

  // permutation matrix
  Epetra_SerialDenseMatrix P(true);
  P.Shape(N, N);
  for (int i = 0; i < N; ++i)
  {
    int w = jpvt[i];
    P(w - 1, i) = 1.;
  }

  // annul the under-diagonal elements of B
  // loop of columns
  for (int i = 0; i < N; ++i)
  {
    // loop of rows
    for (int j = i + 1; j < N; ++j) tmpB(j, i) = 0.;
  }

  // the new A:= Q**T A P
  Epetra_SerialDenseMatrix A_tmp;
  A_tmp.Shape(N, N);
  // A_tt.Multiply('T','N',1.,Q_qr_tt,A,0.);
  A_tmp.Multiply('T', 'N', 1., Q_new, tmpA, 0.);

  Epetra_SerialDenseMatrix A_new;
  A_new.Shape(N, N);
  A_new.Multiply('N', 'N', 1., A_tmp, P, 0.);

  a = A_new.A();

  //--------------------------------------------------------
  // STEP 2
  // transform A to a upper hessenberg matrix and keep B as
  // an upper diagonal matrix
  //--------------------------------------------------------

  // balance problem to obtain better accuracy
  char job = 'P';
  int ILO;
  int IHI;
  double lscale[N];
  double rscale[N];
  double work0[6 * N];
  dggbal(&job, &N, a, &N, b, &N, &ILO, &IHI, lscale, rscale, work0, &info);
  if (info != 0) dserror("error dggbal");

  job = 'E';
  char COMPQ = 'I';
  char COMPZ = 'I';

  int lwork = 0;
  if (N == 1)
    lwork = 1;
  else
  {
    lwork = N;
  }
  double work[lwork];

  Epetra_SerialDenseMatrix A1(true);
  Epetra_SerialDenseMatrix A2(true);
  A1.Shape(N, N);
  A2.Shape(N, N);
  double* Q = A1.A();
  int LDQ = A1.LDA();
  double* Z = A2.A();
  int LDZ = A2.LDA();

  dgghrd(&COMPQ, &COMPZ, &N, &ILO, &IHI, a, &LDA, b, &LDB, Q, &LDQ, Z, &LDZ, &info);

  if (info < 0)
    std::cout << "Lapack algorithm dgghrd: The " << info << "-th argument had an illegal value"
              << std::endl;

  //--------------------------------------------------------
  // STEP 3
  // transform A which is an upper hessenberg matrix to an upper
  // diagonal matrix and keep B an upper diagonal matrix via a
  // QZ-transformation
  //--------------------------------------------------------
  // vectors which contain the eigenvalues of the problem
  Epetra_SerialDenseVector L1(true);
  Epetra_SerialDenseVector L2(true);
  Epetra_SerialDenseVector L3(true);
  L1.Shape(N, 1);
  L2.Shape(N, 1);
  L3.Shape(N, 1);
  double* ALPHAR = L1.A();
  double* ALPHAI = L2.A();
  double* BETA = L3.A();

  int LDH = A_new.LDA();
  int LDT = tmpB.LDA();

  char COMPQ2 = 'V';
  char COMPZ2 = 'V';

  dhgeqz(&job, &COMPQ2, &COMPZ2, &N, &ILO, &IHI, a, &LDH, b, &LDT, ALPHAR, ALPHAI, BETA, Q, &LDQ, Z,
      &LDZ, work, &lwork, &info);

  if (info < 0)
    std::cout << "Lapack algorithm dhgeqz: The " << info << "-th argument haa an illegal value!"
              << std::endl;
  else if (info > N)
    std::cout << "Lapack algorithm dhgeqz: The QZ iteration did not converge. (H,T) is not in "
                 "Schur Form, but the Eigenvalues should be correct!"
              << std::endl;

  /*cout << "--------Final----------" << std::endl;
   std::cout << std::setprecision(16) << "A 2" << A_new << std::endl;
   std::cout << std::setprecision(16) <<  "B 2" << tmpB << std::endl;
   std::cout << std::setprecision(16) << "Q 2 " << Q_2 << std::endl;
   std::cout << std::setprecision(16) << "Z 2 " << Z_2 << std::endl;*/

  double maxlambda = 0.;
  for (int i = 0; i < N; ++i)
  {
    if (BETA[i] > 1e-13)
    {
      // Eigenvalues:
      // std::cout << "lambda " << i << ":  " <<  ALPHAR[i]/BETA[i] << std::endl;
      maxlambda = std::max(ALPHAR[i] / BETA[i], maxlambda);
    }
    if (ALPHAI[i] > 1e-12)
    {
      std::cout << " Warning: you have an imaginary EW " << ALPHAI[i] << std::endl;
    }
  }
  return maxlambda;
}
