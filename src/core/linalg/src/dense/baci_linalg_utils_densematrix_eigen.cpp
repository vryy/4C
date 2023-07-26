/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of eigenvalue methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_linalg_fortran_definitions.h"
#include "baci_linalg_utils_densematrix_eigen.H"
#include "baci_utils_exceptions.H"

/*----------------------------------------------------------------------*
 |  compute all eigenvalues of a real symmetric matrix A        lw 04/08|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigenValues(CORE::LINALG::SerialDenseMatrix::Base& A,
    CORE::LINALG::SerialDenseVector::Base& L, const bool postproc)
{
  CORE::LINALG::SymmetricEigen(A, L, 'N', postproc);
}

/*----------------------------------------------------------------------*
 |  compute all eigenvalues and eigenvectors of a real symmetric        |
 |  matrix A (eigenvectors are stored in A, i.e. original matrix        |
 |  is destroyed!!!)                                            lw 04/08|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigenProblem(CORE::LINALG::SerialDenseMatrix::Base& A,
    CORE::LINALG::SerialDenseVector::Base& L, const bool postproc)
{
  CORE::LINALG::SymmetricEigen(A, L, 'V', postproc);
}

/*----------------------------------------------------------------------*
 |  compute all eigenvalues and, optionally,                            |
 |  eigenvectors of a real symmetric matrix A                  maf 06/07|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricEigen(CORE::LINALG::SerialDenseMatrix::Base& A,
    CORE::LINALG::SerialDenseVector::Base& L, const char jobz, const bool postproc)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != L.Length()) dserror("Dimension of eigenvalues does not match");

  double* a = A.A();
  double* w = L.A();
  const char uplo = {'U'};
  const int lda = A.LDA();
  const int dim = A.M();

  int lwork = 0;
  if (dim == 1)
    lwork = 1;
  else
  {
    if (jobz == 'N')
      lwork = 3 * dim + 1;
    else if (jobz == 'V')
      lwork = 2 * dim * dim + 6 * dim + 1;
  }
  std::vector<double> work(lwork);
  int info = 0;

  Teuchos::LAPACK<int, double> lapack;
  lapack.SYEV(jobz, uplo, dim, a, lda, w, work.data(), lwork, &info);

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
double LINALG::GeneralizedEigen(
    CORE::LINALG::SerialDenseMatrix::Base& A, CORE::LINALG::SerialDenseMatrix::Base& B)
{
  CORE::LINALG::SerialDenseMatrix tmpA(A);
  CORE::LINALG::SerialDenseMatrix tmpB(B);

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
  std::vector<int> jpvt(N);
  // factor uses for calculating orthogonal matrix Q
  std::vector<double> tau(N);
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
  std::vector<double> work1(lwork1);
  int info;
  dgeqp3(&N, &N, b, &LDB, jpvt.data(), tau.data(), work1.data(), &lwork1, &info);

  if (info < 0)
    std::cout << "Lapack algorithm dgeqp3: The " << info << "-th argument had an illegal value"
              << std::endl;

  // calculate the matrix Q from multiplying householder transformations
  // Q = H(1)*H(2) ... H(k)
  // H(i) = I - tau-v*v**T
  // v is a vector with v(1:i-1) = 0 and v(i+1:m) is stored on exit in B(i+1:m,i)

  // Q is initialized as an unit matrix
  CORE::LINALG::SerialDenseMatrix Q_new(true);
  Q_new.Shape(N, N);
  for (int i = 0; i < N; ++i) Q_new(i, i) = 1.0;

  for (int i = 0; i < N; ++i)
  {
    CORE::LINALG::SerialDenseVector v;
    v.Shape(N, 1);
    v(i, 0) = 1.;
    for (int j = i + 1; j < N; ++j) v(j, 0) = tmpB(j, i);

    CORE::LINALG::SerialDenseMatrix H;
    H.Shape(N, N);

    H.Multiply('N', 'T', tau[i], v, v, 0.);
    H.Scale(-1.);
    for (int k = 0; k < N; ++k) H(k, k) = 1. + H(k, k);

    CORE::LINALG::SerialDenseMatrix Q_help;
    Q_help.Shape(N, N);
    Q_new.Apply(H, Q_help);
    Q_new = Q_help;
  }

  // permutation matrix
  CORE::LINALG::SerialDenseMatrix P(true);
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
  CORE::LINALG::SerialDenseMatrix A_tmp;
  A_tmp.Shape(N, N);
  // A_tt.Multiply('T','N',1.,Q_qr_tt,A,0.);
  A_tmp.Multiply('T', 'N', 1., Q_new, tmpA, 0.);

  CORE::LINALG::SerialDenseMatrix A_new;
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
  std::vector<double> lscale(N);
  std::vector<double> rscale(N);
  std::vector<double> work0(6 * N);
  dggbal(&job, &N, a, &N, b, &N, &ILO, &IHI, lscale.data(), rscale.data(), work0.data(), &info);
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
  std::vector<double> work(lwork);

  CORE::LINALG::SerialDenseMatrix A1(true);
  CORE::LINALG::SerialDenseMatrix A2(true);
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
  CORE::LINALG::SerialDenseVector L1(true);
  CORE::LINALG::SerialDenseVector L2(true);
  CORE::LINALG::SerialDenseVector L3(true);
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
      &LDZ, work.data(), &lwork, &info);

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
