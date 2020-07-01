/*----------------------------------------------------------------------*/
/*! \file

\brief collection of math tools for the interface determination of trv1o meshes

    ML      math library for the interface computation

    !WARNING: Except Tolerances not used at the moment
    (remove this comment and change level as soon as this functionality is tested again!)


\level 3

*----------------------------------------------------------------------*/


#include "../drt_geometry/intersection_math.H"
#include "../linalg/linalg_serialdensematrix.H"


void GEO::test_svdcmp(LINALG::Matrix<3, 3>& A, LINALG::Matrix<3, 3>& U, LINALG::Matrix<3, 1>& W,
    LINALG::Matrix<3, 3>& V, int dim)
{
  LINALG::Matrix<3, 3> H1;
  LINALG::Matrix<3, 3> H2;

  printf("W U\n");
  for (int i = 0; i < dim; i++)
  {
    printf("W = %f\t", W(i));
    for (int j = 0; j < dim; j++)
    {
      printf("U = %f\t", U(i, j));
    }
    printf("\n");
  }
  printf("\n");

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      H1(i, j) = U(i, j) * W(j);
    }
  }

  printf("H1\n");
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("H1 = %f\t", H1(i, j));
    }
    printf("\n");
  }
  printf("\n");

  printf("V\n");
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("V = %f\t", V(i, j));
    }
    printf("\n");
  }
  printf("\n");

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      H2(i, j) = 0.0;
      for (int k = 0; k < dim; k++)
      {
        H2(i, j) += H1(i, k) * V(j, k);
      }
    }
  }

  printf("system matrix\n");
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("A = %f\t", A(i, j));
    }
    printf("\n");
  }
  printf("\n");

  printf("system matrix SVD\n");
  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < dim; j++)
    {
      printf("H2 = %f\t", H2(i, j));
    }
    printf("\n");
  }
  printf("\n");
}

void GEO::svdcmpSerialDense(
    Epetra_SerialDenseMatrix& A, Epetra_SerialDenseMatrix& W, Epetra_SerialDenseMatrix& V)
{
  // Dimensionen der Matrix A herausfinden
  const int n = A.N();
  const int m = A.M();

  // Prüfen ob die W und V entsprechend richtige Dimension haben
  if (!((W.M() == n) && (V.M() == n && V.N() == n)))
    dserror("Dimensionen der Matrizen nicht korrekt");

  LINALG::SerialDenseMatrix rv1(n, 1);

  // Householder reduction to bidiagonal form.
  double g = 0.0;
  double scale = 0.0;
  double anorm = 0.0;
  for (int i = 0; i < n; ++i)
  {
    const int l = i + 1;
    rv1(i, 0) = scale * g;
    g = scale = 0.0;
    if (i < m)
    {
      for (int k = i; k < m; k++)
      {
        scale += fabs(A(k, i));
      }
      if (scale != 0.0)
      {
        double s = 0.0;
        for (int k = i; k < m; ++k)
        {
          A(k, i) /= scale;
          s += std::pow(A(k, i), 2);
        }
        const double f = A(i, i);
        g = -XSIGN(sqrt(s), f);
        const double h = f * g - s;
        A(i, i) = f - g;
        for (int j = l; j < n; ++j)
        {
          double s = 0.0;
          for (int k = i; k < m; ++k) s += A(k, i) * A(k, j);
          const double f = s / h;
          for (int k = i; k < m; ++k) A(k, j) += f * A(k, i);
        }
        for (int k = i; k < m; ++k) A(k, i) *= scale;
      }
    }
    W(i, 0) = scale * g;
    g = scale = 0.0;
    if (i < m && i != (n - 1))
    {
      for (int k = l; k < n; k++)
      {
        scale += fabs(A(i, k));
      }
      if (scale)
      {
        double s = 0.0;
        for (int k = l; k < n; k++)
        {
          A(i, k) /= scale;
          s += std::pow(A(i, k), 2);
        }
        const double f = A(i, l);
        g = -XSIGN(sqrt(s), f);
        const double h = f * g - s;
        const double h_inv = 1.0 / h;
        A(i, l) = f - g;
        for (int k = l; k < n; k++) rv1(k, 0) = A(i, k) * h_inv;
        for (int j = l; j < m; ++j)
        {
          double s = 0.0;
          for (int k = l; k < n; ++k) s += A(j, k) * A(i, k);
          for (int k = l; k < n; ++k) A(j, k) += s * rv1(k, 0);
        }
        for (int k = l; k < n; k++) A(i, k) *= scale;
      }
    }
    anorm = std::max(anorm, (fabs(W(i, 0)) + fabs(rv1(i, 0))));
  }
  // Accumulation of right-hand transformations.
  for (int i = (n - 1); i >= 0; i--)
  {
    const int l = i + 1;
    if (i < n)
    {
      if (g)
      {
        // Double division to avoid possible underflow.
        const double g_inv = 1.0 / g;
        for (int j = l; j < n; j++) V(j, i) = (A(i, j) / A(i, l)) * g_inv;
        for (int j = l; j < n; j++)
        {
          double s = 0.0;
          for (int k = l; k < n; k++) s += A(i, k) * V(k, j);
          for (int k = l; k < n; k++) V(k, j) += s * V(k, i);
        }
      }
      for (int j = l; j < n; j++) V(i, j) = V(j, i) = 0.0;
    }
    V(i, i) = 1.0;
    g = rv1(i, 0);
  }
  // Accumulation of left-hand transformations.
  for (int i = std::min((m - 1), (n - 1)); i >= 0; i--)
  {
    const int l = i + 1;
    const double g = W(i, 0);
    for (int j = l; j < n; ++j) A(i, j) = 0.0;
    if (g)
    {
      const double g_inv = 1.0 / g;
      for (int j = l; j < n; ++j)
      {
        double s = 0.0;
        for (int k = l; k < m; ++k) s += A(k, i) * A(k, j);
        const double f = (s / A(i, i)) * g_inv;
        for (int k = i; k < m; ++k) A(k, j) += f * A(k, i);
      }
      for (int j = i; j < m; ++j) A(j, i) *= g_inv;
    }
    else
      for (int j = i; j < m; ++j) A(j, i) = 0.0;
    ++A(i, i);
  }
  // Diagonalization of the bidiagonal form: Loop over
  for (int k = (n - 1); k >= 0; k--)
  {
    // singular values, and over allowed iterations.
    for (int its = 1; its <= 30; its++)
    {
      // Test for splitting.
      bool flag = true;
      int nm;
      int l = k;
      for (l = k; l >= 0; l--)
      {
        // Note that rv1(1) is always zero.
        nm = l - 1;
        if ((double)(fabs(rv1(l, 0)) + anorm) == anorm)
        {
          flag = false;
          break;
        }
        if ((double)(fabs(W(nm, 0)) + anorm) == anorm) break;
      }
      if (flag)
      {
        // Cancellation of rv1(l), if l > 1.
        double c = 0.0;
        double s = 1.0;
        for (int i = l; i <= k; i++)
        {
          const double f = s * rv1(i, 0);
          rv1(i, 0) *= c;
          if ((double)(fabs(f) + anorm) == anorm) break;
          const double g = W(i, 0);
          const double h = GEO::pythagoras(f, g);
          W(i, 0) = h;
          const double h_inv = 1.0 / h;
          c = g * h_inv;
          s = -f * h_inv;
          for (int j = 0; j < m; j++)
          {
            const double y = A(j, nm);
            const double z = A(j, i);
            A(j, nm) = y * c + z * s;
            A(j, i) = z * c - y * s;
          }
        }
      }
      const double z = W(k, 0);
      // Convergence.
      if (l == k)
      {
        // Singular value is made nonnegative.
        if (z < 0.0)
        {
          W(k, 0) = -z;
          for (int j = 0; j < n; j++) V(j, k) = -V(j, k);
        }
        break;
      }
      if (its == 30) dserror("no convergence in 30 svdcmp iterations");
      // Shift from bottom 2-by-2 minor.
      double x = W(l, 0);
      const int km = k - 1;
      const double y = W(km, 0);
      double g = rv1(km, 0);
      double h = rv1(k, 0);
      double f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = GEO::pythagoras(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + XSIGN(g, f))) - h)) / x;
      // Next QR transformation:
      double c = 1.0;
      double s = 1.0;
      for (int j = l; j <= km; j++)
      {
        const int i = j + 1;
        double g = rv1(i, 0);
        double y = W(i, 0);
        double h = s * g;
        g = c * g;
        double z = GEO::pythagoras(f, h);
        rv1(j, 0) = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (int jj = 0; jj < n; jj++)
        {
          const double x = V(jj, j);
          const double z = V(jj, i);
          V(jj, j) = x * c + z * s;
          V(jj, i) = z * c - x * s;
        }
        z = GEO::pythagoras(f, h);
        // Rotation can be arbitrary if z = 0.
        W(j, 0) = z;
        if (z)
        {
          const double z_inv = 1.0 / z;
          c = f * z_inv;
          s = h * z_inv;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (int jj = 0; jj < m; jj++)
        {
          const double y = A(jj, j);
          const double z = A(jj, i);
          A(jj, j) = y * c + z * s;
          A(jj, i) = z * c - y * s;
        }
      }
      rv1(l, 0) = 0.0;
      rv1(k, 0) = f;
      W(k, 0) = x;
    }
  }
}
