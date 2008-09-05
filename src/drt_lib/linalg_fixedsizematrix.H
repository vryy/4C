
#ifndef FIXED_SIZE_SERIAL_DENSE_MATRIX_H
#define FIXED_SIZE_SERIAL_DENSE_MATRIX_H

#include <Epetra_LAPACK.h>
#include <Epetra_BLAS.h>
#include <Epetra_SerialDenseMatrix.h>

#include <cstring>
#include <vector>
#include <cmath>
#include <ostream>
#include <iostream>

#include "../drt_lib/drt_dserror.H"


/// Fixed size serial dense matrix definitions
namespace LINALG
{

  /// compile time error check for fixed size serial dense classes
  namespace DENSEERROR
  {

    /// compile time error definitions
    template <bool expr>
    struct Checker
    {
      static inline void Matrix_dimensions_cannot_be_zero() { };
      static inline void Cannot_call_1D_access_function_on_2D_matrix() { };
      static inline void Cannot_compute_determinant_of_nonsquare_matrix() { };
      static inline void Cannot_compute_inverse_of_nonsquare_matrix() { };
      static inline void Transpose_argument_must_be_N_or_T() { };
      static inline void Matrix_size_in_solver_must_be_square() { };
    };

    /// compile time error definitions: missing functions raise errors
    template<>
    struct Checker<false>
    {
    };

  }

  /// Internal templated serial dense functions
  namespace DENSEFUNCTIONS
  {
    // functions taking double*
    template<unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(double* out, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double outfac, double* out, const double infac, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(double* out, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(const double outfac, double* out, const double infac, const double* left, const double* right);

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2, int fac1, int fac2>
    inline void multiply(const double outfac, double* out, const double infac, const double* left, const double* right);

    template <unsigned int i, unsigned int j>
    inline double invert(double* out, const double* in);

    template <unsigned int i, unsigned int j>
    inline double invert(double* mat);

    template <unsigned int i, unsigned int j>
    inline double determinant(const double* mat);

    template <unsigned int i, unsigned int j>
    inline void addMatrix(double* out, const double* left,
                          const double* right);

    template <unsigned int i, unsigned int j>
    inline void addMatrix(double* out, const double leftfac, const double* left, const double rightfac, const double* right);

    template <unsigned int i, unsigned int j>
    inline void addMatrix(const double outfac, double* out, const double leftfac, const double* left,
                          const double rightfac, const double* right);

    template <unsigned int i, unsigned int j, int fac1, int fac2>
    inline void addMatrix(const double outfac, double* out, const double leftfac, const double* left,
                          const double rightfac, const double* right);

    template <unsigned int i, unsigned int j>
    inline void addMatrixTo(double* out, const double* in);

    template <unsigned int i, unsigned int j>
    inline void addMatrixTo(const double outfac, double* out, const double infac, const double* in);

    template <unsigned int i, unsigned int j, int fac1, int fac2>
    inline void addMatrixTo(const double outfac, double* out, const double infac, const double* in);

    template <unsigned int i, unsigned int j>
    inline void scaleMatrix(const double scalar, double* matrix);

    template <unsigned int i, unsigned int j, int fac>
    inline void scaleMatrix(const double scalar, double* matrix);

    template <unsigned int i, unsigned int j>
    inline void copyMatrix(double* out, const double* in);

    template <unsigned int i, unsigned int j>
    inline double dot(const double* left, const double* right);

    template<unsigned int i, unsigned int j>
    inline void clearMatrix(double* mat);

    template <unsigned int i, unsigned int j>
    inline void putScalar(const double scalar, double* data);

    template <unsigned int i, unsigned int j>
    inline void abs(double* dest, const double* src);

    template <unsigned int i, unsigned int j>
    inline void reciprocal(double* dest, const double* src);

    template <unsigned int i,  unsigned int j>
    inline double norm1(const double* data);

    template <unsigned int i,  unsigned int j>
    inline double norm2(const double* data);

    template <unsigned int i,  unsigned int j>
    inline double normInf(const double* data);

    template <unsigned int i,  unsigned int j>
    inline double minValue(const double* data);

    template <unsigned int i,  unsigned int j>
    inline double maxValue(const double* data);

    template <unsigned int i,  unsigned int j>
    inline double meanValue(const double* data);

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(double* dest, const double* left, const double* right);

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(const double destfac, double* dest, const double prodfac, const double* left, const double* right);

    // functions taking Epetra_SerialDenseMatrix. These are only wrapper
    template<unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
        dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      multiply<i,j,k>(out.A(), left.A(), right.A());
    }

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double outfac, Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
        dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      multiply<i,j,k>(outfac, out.A(), infac, left.A(), right.A());
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      switch (trans1) {
      case 'N':
        switch (trans2) {
        case 'N':
          if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.M() or out.N() != right.M() or left.N() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      case 'T':
        switch (trans2) {
        case 'N':
          if (out.M() != left.N() or out.N() != right.N() or left.M() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.N() or out.N() != right.M() or left.M() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      }
#endif
      multiply<i,j,k,trans1,trans2>(out.A(), left.A(), right.A());
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline void multiply(const double outfac, Epetra_SerialDenseMatrix& out, const double infac,
                         const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      switch (trans1) {
      case 'N':
        switch (trans2) {
        case 'N':
          if (out.M() != left.M() or out.N() != right.N() or left.N() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.M() or out.N() != right.M() or left.N() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i) * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      case 'T':
        switch (trans2) {
        case 'N':
          if (out.M() != left.N() or out.N() != right.N() or left.M() != right.M())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        case 'T':
          if (out.M() != left.N() or out.N() != right.M() or left.M() != right.N())
            dserror("Invalid matrix sizes for multiplication, (%i,%i) = (%i,%i)^T * (%i,%i)^T",
                    out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
          break;
        } break;
      }
#endif
      multiply<i,j,k,trans1,trans2>(outfac, out.A(), infac, left.A(), right.A());
    }

    template <unsigned int i, unsigned int j>
    inline double invert(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != out.N() or in.M() != in.N() or out.M() != in.M())
        dserror("Invalid matrix sizes for inversion, (%i,%i) = inv( (%i,%i) )",out.M(),out.N(),in.M(),in.N());
#endif
      return invert<i,j>(out.A(), in.A());
    }

    template <unsigned int i, unsigned int j>
    inline double invert(Epetra_SerialDenseMatrix& mat)
    {
#ifdef DEBUG
      if (mat.M() != mat.N())
        dserror("Invalid matrix sizes for inversion, inv( (%i,%i) )",mat.M(),mat.N());
#endif
      return invert<i,j>(mat.A());
    }

    template <unsigned int i, unsigned int j>
    inline double determinant(const Epetra_SerialDenseMatrix& mat)
    {
#ifdef DEBUG
      if (mat.M() != mat.N())
        dserror("Invalid matrix sizes for determinant, inv( (%i,%i) )",mat.M(),mat.N());
#endif
      return determinant<i,j>(mat.A());
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrix(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& left,
                          const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or left.M() != right.M() or out.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      addMatrix<i,j>(out.A(), left.A(), right.A());
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrix(Epetra_SerialDenseMatrix& out, const double leftfac, const Epetra_SerialDenseMatrix& left,
                          const double rightfac, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or left.M() != right.M() or out.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      addMatrix<i,j>(out.A(), leftfac, left.A(), rightfac, right.A());
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrix(const double outfac, Epetra_SerialDenseMatrix& out, const double leftfac,
                          const Epetra_SerialDenseMatrix& left, const double rightfac, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (out.M() != left.M() or left.M() != right.M() or out.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) = (%i,%i) + (%i,%i)",
                out.M(),out.N(), left.M(), left.N(), right.M(), right.N());
#endif
      addMatrix<i,j>(outfac, out.A(), leftfac, left.A(), rightfac, right.A());
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrixTo(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != in.M() or out.N() != in.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)",
                out.M(),out.N(), in.M(), in.N());
#endif
      addMatrixTo<i,j>(out.A(), in.A());
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrixTo(const double outfac, Epetra_SerialDenseMatrix& out, const double infac, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != in.M() or out.N() != in.N())
        dserror("Invalid matrix sizes for addition, (%i,%i) += (%i,%i)",
                out.M(),out.N(), in.M(), in.N());
#endif
      addMatrixTo<i,j>(outfac, out.A(), infac, in.A());
    }

    template <unsigned int i, unsigned int j>
    inline void scaleMatrix(const double scalar, Epetra_SerialDenseMatrix& matrix)
    { scaleMatrix<i,j>(scalar, matrix.A()); }

    template <unsigned int i, unsigned int j>
    inline void copyMatrix(Epetra_SerialDenseMatrix& out, const Epetra_SerialDenseMatrix& in)
    {
#ifdef DEBUG
      if (out.M() != in.M() or out.N() != in.N())
        dserror("Invalid matrix sizes for copy, (%i,%i) = (%i,%i)",
                out.M(),out.N(), in.M(), in.N());
#endif
      copyMatrix<i,j>(out.A(), in.A());
    }

    template <unsigned int i, unsigned int j>
    inline double dot(const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    { return dot<i,j>(left.A(), right.A()); }

    template<unsigned int i, unsigned int j>
    inline void clearMatrix(Epetra_SerialDenseMatrix& mat)
    { clearMatrix<i,j>(mat.A()); }

    template <unsigned int i, unsigned int j>
    inline void putScalar(const double scalar, Epetra_SerialDenseMatrix& data)
    { putScalar<i,j>(scalar, data.A()); }

    template <unsigned int i, unsigned int j>
    inline void abs(Epetra_SerialDenseMatrix& dest, const Epetra_SerialDenseMatrix& src)
    {
#ifdef DEBUG
      if (dest.M() != dest.M() or src.N() != src.N())
        dserror("Invalid matrix sizes for abs, (%i,%i) = abs( (%i,%i) )",
                dest.M(),dest.N(), src.M(), src.N());
#endif
      abs<i,j>(dest.A(), src.A());
    }

    template <unsigned int i, unsigned int j>
    inline void reciprocal(Epetra_SerialDenseMatrix& dest, const Epetra_SerialDenseMatrix& src)
    {
#ifdef DEBUG
      if (dest.M() != dest.M() or src.N() != src.N())
        dserror("Invalid matrix sizes for reciprocal, (%i,%i) = reciprocal( (%i,%i) )",
                dest.M(),dest.N(), src.M(), src.N());
#endif
      reciprocal<i,j>(dest.A(), src.A());
    }

    template <unsigned int i,  unsigned int j>
    inline double norm1(const Epetra_SerialDenseMatrix& data)
    { return norm1<i,j>(data.A()); }

    template <unsigned int i,  unsigned int j>
    inline double norm2(const Epetra_SerialDenseMatrix& data)
    { return norm2<i,j>(data.A()); }

    template <unsigned int i,  unsigned int j>
    inline double normInf(const Epetra_SerialDenseMatrix& data)
    { return normInf<i,j>(data.A()); }

    template <unsigned int i,  unsigned int j>
    inline double minValue(const Epetra_SerialDenseMatrix& data)
    { return minValue<i,j>(data.A()); }

    template <unsigned int i,  unsigned int j>
    inline double maxValue(const Epetra_SerialDenseMatrix& data)
    { return maxValue<i,j>(data.A()); }

    template <unsigned int i,  unsigned int j>
    inline double meanValue(const Epetra_SerialDenseMatrix& data)
    { return meanValue<i,j>(data.A()); }

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(Epetra_SerialDenseMatrix& dest, const Epetra_SerialDenseMatrix& left,
                                   const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (dest.M() != left.M() or left.M() != right.M() or dest.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for reciprocal multiplication, (%i,%i) = reciprocalMultiply( (%i,%i), (%i,%i) )",
                dest.M(),dest.N(), left.M(), left.N(), right.M(), right.N());
#endif
      reciprocalMultiply<i,j>(dest.A(), left.A(), right.A());
    }

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(const double destfac, Epetra_SerialDenseMatrix& dest, const double prodfac,
                                   const Epetra_SerialDenseMatrix& left, const Epetra_SerialDenseMatrix& right)
    {
#ifdef DEBUG
      if (dest.M() != left.M() or left.M() != right.M() or dest.N() != left.N() or left.N() != right.N())
        dserror("Invalid matrix sizes for reciprocal multiplication, (%i,%i) = reciprocalMultiply( (%i,%i), (%i,%i) )",
                dest.M(),dest.N(), left.M(), left.N(), right.M(), right.N());
#endif
      reciprocalMultiply<i,j>(destfac, dest.A(), prodfac, left.A(), right.A());
    }


/*
 * multiplication functions
 */

    template<unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(double* dest, const double* left, const double* right)
    {
      multiply<i,j,k,'N','N',0,1>(0.0,dest,1.0,left,right);
    }

    template <unsigned int i, unsigned int j, unsigned int k>
    inline void multiply(const double destfac, double* dest, double infac, const double* left, const double* right)
    {
      multiply<i,j,k,'N','N'>(destfac,dest,infac,left,right);
    }


    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline
    void multiply(double* dest, const double* left, const double* right)
    {
      DENSEERROR::Checker<(trans1 == 'N' or trans1 == 'T') and (trans2 == 'N' or trans2 == 'T')>::Transpose_argument_must_be_N_or_T();
      multiply<i,j,k,trans1,trans2,0,1>(0.0, dest, 1.0, left, right);
    }

    template <unsigned int i, unsigned int j, unsigned int k, char trans1, char trans2>
    inline
    void multiply(const double destfac, double* dest, const double infac, const double* left, const double* right)
    {
      DENSEERROR::Checker<(trans1 == 'N' or trans1 == 'T') and (trans2 == 'N' or trans2 == 'T')>::Transpose_argument_must_be_N_or_T();
      if (infac == 0.0) {
        scaleMatrix<i,k>(destfac,dest);
      }
      else if (destfac == -1.0) {
        if (infac == -1.0)
          multiply<i,j,k,trans1,trans2,-1,-1>(destfac, dest, infac, left, right);
        else if (infac == 1.0)
          multiply<i,j,k,trans1,trans2,-1,1>(destfac, dest, infac, left, right);
        else
          multiply<i,j,k,trans1,trans2,-1,2>(destfac, dest, infac, left, right);
      }
      else if (destfac == 0.0) {
        if (infac == -1.0)
          multiply<i,j,k,trans1,trans2,0,-1>(destfac, dest, infac, left, right);
        else if (infac == 1.0)
          multiply<i,j,k,trans1,trans2,0,1>(destfac, dest, infac, left, right);
        else
          multiply<i,j,k,trans1,trans2,0,2>(destfac, dest, infac, left, right);
      }
      else if (destfac == 1.0) {
        if (infac == -1.0)
          multiply<i,j,k,trans1,trans2,1,-1>(destfac, dest, infac, left, right);
        else if (infac == 1.0)
          multiply<i,j,k,trans1,trans2,1,1>(destfac, dest, infac, left, right);
        else
          multiply<i,j,k,trans1,trans2,1,2>(destfac, dest, infac, left, right);
      }
      else {
        if (infac == -1.0)
          multiply<i,j,k,trans1,trans2,2,-1>(destfac, dest, infac, left, right);
        else if (infac == 1.0)
          multiply<i,j,k,trans1,trans2,2,1>(destfac, dest, infac, left, right);
        else
          multiply<i,j,k,trans1,trans2,2,2>(destfac, dest, infac, left, right);
      }
    }


/*
 * The function that does the real work. It isn't exactly easy to
 * read, but it handles all the cases.
 */

    template <unsigned i, unsigned j, unsigned k, char trans1, char trans2, int fac1, int fac2>
    inline
    void multiply(const double destfac, double* dest, const double infac, const double* left, const double* right)
    {
      // only used when fac1 != 0 and fac2 == 2, I hope that the compiler
      // throws it out otherwise...
      double tmp;
      for (unsigned int d1 = 0; d1 < k; ++d1) {
        for (unsigned int d2 = 0; d2 < i; ++d2 ) {
          switch(fac1) {
          case 0:
            switch(fac2) {
            case -1:
              *dest = -(*left)*(*right);
              break;
            default:
              *dest = (*left)*(*right);
              break;
            }
            break;
          case -1:
            switch(fac2) {
            case -1:
              *dest = -(*left)*(*right) - *dest;
              break;
            case 1:
              *dest = (*left)*(*right) - *dest;
              break;
            default: // 2
              tmp = (*left)*(*right);
              *dest = -(*dest);
              break;
            }
            break;
          case 1:
            switch(fac2) {
            case -1:
              *dest -= (*left)*(*right);
              break;
            case 1:
              *dest += (*left)*(*right);
              break;
            default: // 2
              tmp = (*left)*(*right);
              break;
            }
            break;
          default: // 2
            *dest *= destfac;
            switch(fac2) {
            case -1:
              *dest -= (*left)*(*right);
              break;
            case 1:
              *dest += (*left)*(*right);
              break;
            default: // 2
              tmp = (*left)*(*right);
              break;
            }
            break;
          }
          for (unsigned int d3 = 1; d3 < j; ++d3 ) {
            switch(trans1) {
            case 'N':
              left += i;
              break;
            case 'T':
              ++left;
              break;
            }
            switch(trans2) {
            case 'N':
              ++right;
              break;
            case 'T':
              right += k;
              break;
            }
            switch(fac1) {
            case 0:
              switch(fac2) {
              case -1:
                *dest -= (*left)*(*right);
                break;
              default:
                *dest += (*left)*(*right);
                break;
              }
              break;
            default: // 1 or 2
              switch(fac2) {
              case -1:
                *dest -= (*left)*(*right);
                break;
              case 1:
                *dest += (*left)*(*right);
                break;
              default: // 2
                tmp += (*left)*(*right);
              }
            }
          }
          switch(fac2) {
          case -1:
          case 1:
            break;
          default: // 2
            switch(fac1) {
            case 0:
              *dest *= infac;
              break;
            default: // 1 or 2
              *dest += tmp * infac;
            }
          }
          ++dest;
          switch(trans1) {
          case 'N':
            left -= i*(j-1) - 1;
            break;
          case 'T':
            ++left;
            break;
          }
          switch(trans2) {
          case 'N':
            right -= j-1;
            break;
          case 'T':
            right -= (j-1)*k;
            break;
          }
        }
        switch(trans1) {
        case 'N':
          left -= i;
          break;
        case 'T':
          left -= i*j;
          break;
        }
        switch(trans2) {
        case 'N':
          right += j;
          break;
        case 'T':
          ++right;
          break;
        }
      }
    }


/*
 * end of multiplication functions
 */

    template <unsigned int i, unsigned int j>
    inline double invert(double* out, const double* in)
    {
      DENSEERROR::Checker<i == j>::Cannot_compute_inverse_of_nonsquare_matrix();
      // do something smart here...
      return 0.0;
    }

    template<>
    inline double invert<1,1>(double* out, const double* in)
    {
      const double det = in[0];
      if (det == 0.0) dserror("Determinant of 1x1 matrix is zero");
      out[0] = 1.0/in[0];
      return det;
    }

    template<>
    inline double invert<2,2>(double* out, const double* in)
    {
      const double det = in[0]*in[1+1*2] - in[1]*in[1*2];
      if (det == 0.0) dserror("Determinant of 2x2 matrix is zero");
      const double invdet = 1.0/det;
      out[0    ] =  invdet*in[1+1*2];
      out[1    ] = -invdet*in[1    ];
      out[  1*2] = -invdet*in[  1*2];
      out[1+1*2] =  invdet*in[0    ];
      return det;
    }


    template<>
    inline double invert<3,3>(double* out, const double* in)
    {
      out[0] = in[1+1*3]*in[2+2*3] - in[2+1*3]*in[1+2*3];
      out[1] = in[2    ]*in[1+2*3] - in[1    ]*in[2+2*3];
      out[2] = in[1    ]*in[2+1*3] - in[2    ]*in[1+1*3];
      const double det = in[0] * out[0] + in[1*3] * out[1] + in[2*3] * out[2];
      //const double det = in[0]*in[1+3*1]*in[2+3*2] +
      //                   in[0+3*1]*in[1+3*2]*in[2+3*0] +
      //                   in[0+3*2]*in[1+3*0]*in[2+3*1] -
      //                   in[0+3*2]*in[1+3*1]*in[2+3*0] -
      //                   in[0+3*0]*in[1+3*2]*in[2+3*1] -
      //                   in[0+3*1]*in[1+3*0]*in[2+3*2];
      if (det == 0.0) dserror("Determinant of 3x3 matrix is zero");
      const double invdet = 1.0/det;
      out[0    ] *= invdet;
      out[1    ] *= invdet;
      out[2    ] *= invdet;
      out[  1*3] = invdet * (in[2+1*3]*in[  2*3] - in[  1*3]*in[2+2*3]);
      out[1+1*3] = invdet * (in[0    ]*in[2+2*3] - in[2    ]*in[  2*3]);
      out[2+1*3] = invdet * (in[2    ]*in[  1*3] - in[0    ]*in[2+1*3]);
      out[  2*3] = invdet * (in[  1*3]*in[1+2*3] - in[1+1*3]*in[  2*3]);
      out[1+2*3] = invdet * (in[1    ]*in[  2*3] - in[0    ]*in[1+2*3]);
      out[2+2*3] = invdet * (in[0    ]*in[1+1*3] - in[1    ]*in[  1*3]);
      return det;
    }

    template <unsigned int i, unsigned int j>
    inline double invert(double* mat)
    {
      DENSEERROR::Checker<i == j>::Cannot_compute_inverse_of_nonsquare_matrix();
      // ...
      return 0.0;
    }

    template <>
    inline double invert<1,1>(double* mat)
    {
      const double det = mat[0];
      if (det == 0.0) dserror("Determinant of 1x1 matrix is zero");
      mat[0] = 1.0/mat[0];
      return det;
    }

    template<>
    inline double invert<2,2>(double* mat)
    {
      double tmp;
      const double det = mat[0]*mat[1+1*2] - mat[1]*mat[1*2];
      if (det == 0.0) dserror("Determinant of 2x2 matrix is zero");
      const double invdet = 1.0/det;
      tmp = mat[0    ]; mat[0    ] = invdet*mat[1+1*2]; mat[1+1*2] = invdet*tmp;
      mat[1    ] *= -invdet;
      mat[  1*2] *= -invdet;
      return det;
    }


    template<>
    inline double invert<3,3>(double* mat)
    {
      const double tmp00 = mat[1+1*3]*mat[2+2*3] - mat[2+1*3]*mat[1+2*3];
      const double tmp10 = mat[2    ]*mat[1+2*3] - mat[1    ]*mat[2+2*3];
      const double tmp20 = mat[1    ]*mat[2+1*3] - mat[2    ]*mat[1+1*3];
      const double det = mat[0] * tmp00 + mat[1*3] * tmp10 + mat[2*3] * tmp20;
      //const double det = mat[0+3*0]*mat[1+3*1]*mat[2+3*2] +
      //                    mat[0+3*1]*mat[1+3*2]*mat[2+3*0] +
      //                    mat[0+3*2]*mat[1+3*0]*mat[2+3*1] -
      //                    mat[0+3*2]*mat[1+3*1]*mat[2+3*0] -
      //                    mat[0+3*0]*mat[1+3*2]*mat[2+3*1] -
      //                    mat[0+3*1]*mat[1+3*0]*mat[2+3*2];
      if (det == 0.0) dserror("Determinant of 3x3 matrix is zero");
      const double invdet = 1.0/det;
      const double tmp01 = mat[  1*3];
      const double tmp11 = mat[1+1*3];
      const double tmp12 = mat[1+2*3];
      mat[  1*3] = invdet * (mat[2+1*3] * mat[  2*3] -    tmp01   * mat[2+2*3]);
      mat[1+1*3] = invdet * (mat[0    ] * mat[2+2*3] - mat[2    ] * mat[  2*3]);
      mat[1+2*3] = invdet * (mat[1    ] * mat[  2*3] - mat[0    ] *    tmp12  );
      mat[2+1*3] = invdet * (mat[2    ] *    tmp01   - mat[0    ] * mat[2+1*3]);
      mat[  2*3] = invdet * (   tmp01   *    tmp12   -    tmp11   * mat[  2*3]);
      mat[2+2*3] = invdet * (mat[0    ] *    tmp11   - mat[1    ] *    tmp01  );
      mat[0    ] = invdet * tmp00;
      mat[1    ] = invdet * tmp10;
      mat[2    ] = invdet * tmp20;
      return det;
    }


    template <unsigned int i, unsigned int j>
    inline double determinant(const double* mat)
    {
      DENSEERROR::Checker<i == j>::Cannot_compute_determinant_of_nonsquare_matrix();

      // taken from src/drt_lib/linalg_utils.cpp: LINALG::DeterminantLU,
      // only with minor changes.
      double* tmp = new double[i*j];
      copyMatrix<i,j>(tmp,mat);
      Epetra_LAPACK lapack;
      std::vector<int> ipiv(j);
      int info;
      lapack.GETRF(i,j,tmp,i,&ipiv[0],&info);
      if (info<0) dserror("Lapack's dgetrf returned %d",info);
      else if (info>0) return 0.0;
      double d = *tmp;
      for (unsigned int c=1; c<j; ++c) d *= tmp[c+i*c];
      // swapping rows of A changes the sign of the determinant, so we have to
      // undo lapack's permutation w.r.t. the determinant
      // note the fortran indexing convention in ipiv
      for (unsigned int c=0; c<j; ++c)
        if (static_cast<unsigned>(ipiv[c])!=c+1) d *= -1.0;
      return d;
    }

    template <>
    inline double determinant<1,1>(const double* mat)
    {
      return *mat;
    }

    template <>
    inline double determinant<2,2>(const double* mat)
    {
      return mat[0    ]*mat[1+1*2] - mat[1    ]*mat[  1*2];
    }

    template <>
    inline double determinant<3,3>(const double* mat)
    {
      return mat[0    ] * (mat[1+1*3]*mat[2+2*3] - mat[2+1*3]*mat[1+2*3]) +
        mat[  1*3] * (mat[2    ]*mat[1+2*3] - mat[1    ]*mat[2+2*3]) +
        mat[  2*3] * (mat[1    ]*mat[2+1*3] - mat[2    ]*mat[1+1*3]);
      //return mat[0]*mat[1+3*1]*mat[2+3*2] +
      //  mat[0+3*1]*mat[1+3*2]*mat[2+3*0] +
      //  mat[0+3*2]*mat[1+3*0]*mat[2+3*1] -
      //  mat[0+3*2]*mat[1+3*1]*mat[2+3*0] -
      //  mat[0+3*0]*mat[1+3*2]*mat[2+3*1] -
      //  mat[0+3*1]*mat[1+3*0]*mat[2+3*2];
    }

    /* add matrices */

    template <unsigned int i, unsigned int j>
    inline
    void addMatrix(double* out, const double* left, const double* right)
    {
      addMatrix<i,j,1,1>(out,1.0,left,1.0,right);
    }

    template <unsigned int i, unsigned int j>
    inline
    void addMatrix(double* out, const double leftfac, const double* left, const double rightfac, const double* right)
    {
      if (leftfac == 1.0) {
        if (rightfac == 1.0)
          addMatrix<i,j,1,1>(0.0, out,1.0,left,1.0,right);
        else
          addMatrix<i,j,1,2>(0.0, out,1.0,left,rightfac,right);
      } else {
        if (rightfac == 1.0)
          addMatrix<i,j,2,1>(0.0, out,leftfac,left,1.0,right);
        else
          addMatrix<i,j,2,2>(0.0, out,leftfac,left,rightfac,right);
      }
    }

    template <unsigned int i, unsigned int j>
    inline
    void addMatrix(const double outfac, double* out, const double leftfac, const double* left,
                   const double rightfac, const double* right)
    {
      if (leftfac == 1.0) {
        if (rightfac == 1.0)
          addMatrix<i,j,1,1>(outfac, out,1.0,left,1.0,right);
        else
          addMatrix<i,j,1,2>(outfac, out,1.0,left,rightfac,right);
      } else {
        if (rightfac == 1.0)
          addMatrix<i,j,2,1>(outfac, out,leftfac,left,1.0,right);
        else
          addMatrix<i,j,2,2>(outfac, out,leftfac,left,rightfac,right);
      }
    }

    template <unsigned int i, unsigned int j, int fac1, int fac2>
    inline
    void addMatrix(const double outfac, double* out, const double leftfac, const double* left,
                   const double rightfac, const double* right)
    {
      switch(fac1) {
      case 1:
        switch(fac2) {
        case 1:
          *out = (*out)*outfac + *left + *right;
          break;
        default:
          *out = (*out)*outfac + *left + rightfac*(*right);
        }
        break;
      default:
        switch(fac2) {
        case 1:
          *out = (*out)*outfac + leftfac*(*left) + *right;
          break;
        default:
          *out = (*out)*outfac + leftfac*(*left) + rightfac*(*right);
        }
      }
      for (unsigned int c = 1; c < i*j; ++c) {
        switch(fac1) {
        case 1:
          switch(fac2) {
          case 1:
            *(++out) = *(++left) + *(++right);
            break;
          default:
            *(++out) = *(++left) + rightfac*(*(++right));
          }
          break;
        default:
          switch(fac2) {
          case 1:
            *(++out) = leftfac*(*(++left)) + *(++right);
            break;
          default:
            *(++out) = leftfac*(*(++left)) + rightfac*(*(++right));
          }
        }
      }
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrixTo(double* out, const double* in)
    {
      addMatrixTo<i,j,1,1>(1.0,out,1.0,in);
    }

    template <unsigned int i, unsigned int j>
    inline void addMatrixTo(const double outfac, double* out, const double infac, const double* in)
    {
      if (infac == 1.0) {
        if (outfac == 1.0)
          addMatrixTo<i,j,1,1>(1.0,out,1.0,in);
        else
          addMatrixTo<i,j,2,1>(outfac,out,1.0,in);
      } else {
        if (outfac == 1.0)
          addMatrixTo<i,j,1,2>(1.0,out,infac,in);
        else
          addMatrixTo<i,j,2,2>(outfac,out,infac,in);
      }
    }

    template <unsigned int i, unsigned int j, int fac1, int fac2>
    inline void addMatrixTo(const double outfac, double* out, const double infac, const double* in)
    {
      switch(fac1) {
      case 1:
        break;
      default:
        *out *= outfac;
      }
      switch(fac2) {
      case 1:
        *out += *in;
        break;
      default:
        *out += infac*(*in);
      }
      for (unsigned int c = 1; c < i*j; ++c) {
        switch(fac1) {
        case 1:
          switch(fac2) {
          case 1:
            *(++out) += *(++in);
            break;
          default:
            *(++out) += infac*(*(++in));
          }
          break;
        default:
          *(++out) *= outfac;
          switch(fac2) {
          case 1:
            *out += *(++in);
            break;
          default:
            *out += infac*(*(++in));
          }
        }
      }
    }


    /* scale a matrix */
    // should the if be avoided here? Probably it only slows the whole
    // thing down...
    template <unsigned int i, unsigned int j>
    inline
    void scaleMatrix(const double scalar, double* mat)
    {
      if (scalar == -1.0)
        scaleMatrix<i,j,-1>(-1.0,mat);
      else if (scalar == 1.0)
        return;
      else
        scaleMatrix<i,j,2>(scalar,mat);
    }

    template <unsigned int i, unsigned int j, int fac>
    inline
    void scaleMatrix(const double scalar, double* mat)
    {
      switch(fac) {
      case -1:
        *mat = -(*mat);
        break;
      case 1:
        return;
      default:
        *mat *= scalar;
      }
      for (unsigned int c = 1; c < i*j; ++c) {
        ++mat;
        switch(fac) {
        case -1:
          *mat = -(*mat);
          break;
        default:
          *mat *= scalar;
        }
      }
    }

    /* copy a matrix */

    template <unsigned int i, unsigned int j>
    inline
    void copyMatrix(double* out, const double* in)
    {
      std::memcpy(out,in,i*j*sizeof(double));
    }

    template <unsigned int i, unsigned int j>
    inline double dot(const double* left, const double* right)
    {
      double res = (*left)*(*right);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++left; ++right;
        res += (*left)*(*right);
      }
      return res;
    }

    // I depend on the double representation of 0.0 to be 0x0. Probably
    // one shouldn't do that, but it allows me to call the fast function memset
    template <unsigned int i, unsigned int j>
    inline void clearMatrix(double* mat)
    {
      std::memset(mat,0,i*j*sizeof(double));
    }


    template <unsigned int i, unsigned int j>
    inline void putScalar(const double scalar, double* data)
    {
      *data = scalar;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++data;
        *data = scalar;
      }
    }

    template <unsigned int i, unsigned int j>
    inline void abs(double* dest, const double* src)
    {
      *dest = *src >= 0 ? *src : -*src;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++dest; ++src;
        *dest = *src >= 0 ? *src : -*src;
      }
    }

    template <unsigned int i, unsigned int j>
    inline void reciprocal(double* dest, const double* src)
    {
      *dest = 1.0/(*src);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++dest; ++src;
        *dest = 1.0/(*src);
      }
    }

    template <unsigned int i,  unsigned int j>
    inline double norm1(const double* data)
    {
      double result = *data;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++data;
        result += *data;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double norm2(const double* data)
    {
      double result = (*data)*(*data);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++data;
        result += (*data)*(*data);
      }
      return std::sqrt(result);
    }

    template <unsigned int i,  unsigned int j>
    inline double normInf(const double* data)
    {
      double result = *data > 0 ? *data : -*data;
      double tmp;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++data;
        tmp = *data >= 0 ? *data : -*data;
        if (tmp > result) result = tmp;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double minValue(const double* data)
    {
      double result = *data;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++data;
        if (*data < result) result = *data;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double maxValue(const double* data)
    {
      double result = *data;
      for (unsigned int c = 1; c < i*j; ++c) {
        ++data;
        if (*data > result) result = *data;
      }
      return result;
    }

    template <unsigned int i,  unsigned int j>
    inline double meanValue(const double* data)
    {
      double result = *data;
      for (unsigned int c = 1; c < i*j; ++c)
        result += *data;
      return result/(i*j);
    }

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(double* dest, const double* left, const double* right)
    {
      *dest = (*right)/(*left);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++dest; ++right; ++left;
        *dest = (*right)/(*left);
      }
    }

    template<unsigned int i, unsigned int j>
    inline void reciprocalMultiply(const double destfac, double* dest, const double prodfac, const double* left, const double* right)
    {
      *dest = *dest*destfac + prodfac*(*right)/(*left);
      for (unsigned int c = 1; c < i*j; ++c) {
        ++dest; ++right; ++left;
        *dest = *dest*destfac + prodfac*(*right)/(*left);
      }
    }

  }


/// Serial dense matrix with templated dimensions for maximum speed
template <unsigned int rows_, unsigned int cols_ /* =1  ?*/ >
class FixedSizeSerialDenseMatrix
{
private:
  double* data_;
  bool isview_;

public:
  FixedSizeSerialDenseMatrix(bool setzero = false);
  explicit FixedSizeSerialDenseMatrix(double* d, bool view=false);
  FixedSizeSerialDenseMatrix(FixedSizeSerialDenseMatrix<rows_,cols_>& source, bool view);
  FixedSizeSerialDenseMatrix(const FixedSizeSerialDenseMatrix<rows_,cols_>& source);
  ~FixedSizeSerialDenseMatrix();

  // methods
  inline const double* A() const { return data_; }
  inline double* A() { return data_; }
  inline unsigned int Rows() const { return rows_; }
  inline unsigned int Columns() const { return cols_; }
  inline unsigned int Cols() const { return cols_; }
  inline unsigned int M() const { return rows_; }
  inline unsigned int N() const { return cols_; }
  inline bool IsInitialized() const { return data_!=NULL; }

  void SetView(double* data);
  void SetView(FixedSizeSerialDenseMatrix<rows_, cols_>& source);
  void SetCopy(double* data);
  void SetCopy(FixedSizeSerialDenseMatrix<rows_, cols_>& source);
  void Copy(double* data);
  void Copy(FixedSizeSerialDenseMatrix<rows_, cols_>& source);


  inline double Determinant() const;
  inline double Invert();
  inline double Invert(const FixedSizeSerialDenseMatrix<rows_, cols_>& other);

  inline void Clear();

  // Epetra-style Functions
  inline void PutScalar(const double scalar);
  inline double Dot(const FixedSizeSerialDenseMatrix<rows_, cols_>& other);
  inline void Abs(const FixedSizeSerialDenseMatrix<rows_, cols_>& other);
  inline void Reciprocal(const FixedSizeSerialDenseMatrix<rows_, cols_>& other);
  inline void Scale(const double scalar);
  inline void Scale(const double scalar, const FixedSizeSerialDenseMatrix<rows_, cols_>& other);
  inline void Update(const FixedSizeSerialDenseMatrix<rows_, cols_>& other);
  inline void Update(const double scalarOther, const FixedSizeSerialDenseMatrix<rows_, cols_>& other,
                     const double scalarThis);
  inline void Update(const double scalarLeft, const FixedSizeSerialDenseMatrix<rows_, cols_>& left,
                     const double scalarRight, const FixedSizeSerialDenseMatrix<rows_, cols_>& right,
                     const double scalarThis);
  inline double Norm1() const;
  inline double Norm2() const;
  inline double NormInf() const;
  inline double MinValue() const;
  inline double MaxValue() const;
  inline double MeanValue() const;
  template <unsigned int inner>
  inline void Multiply(const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right);
  template <unsigned int inner>
  inline void MultiplyNN(const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right);
  template <unsigned int inner>
  inline void MultiplyNT(const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols_, inner>& right);
  template <unsigned int inner>
  inline void MultiplyTN(const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right);
  template <unsigned int inner>
  inline void MultiplyTT(const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                       const FixedSizeSerialDenseMatrix<cols_, inner>& right);

  template <unsigned int inner>
  inline void Multiply(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right, const double scalarThis);
  template <unsigned int inner>
  inline void MultiplyNN(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right, const double scalarThis);
  template <unsigned int inner>
  inline void MultiplyNT(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols_, inner>& right, const double scalarThis);
  template <unsigned int inner>
  inline void MultiplyTN(const double scalarOthers, const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right, const double scalarThis);
  template <unsigned int inner>
  inline void MultiplyTT(const double scalarOthers, const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                       const FixedSizeSerialDenseMatrix<cols_, inner>& right, const double scalarThis);
  inline void ReciprocalMultiply(const double scalarOthers, const FixedSizeSerialDenseMatrix<rows_, cols_>& left,
                                 const FixedSizeSerialDenseMatrix<rows_, cols_>& right, const double scalarThis);
  inline void ReciprocalMultiply(const FixedSizeSerialDenseMatrix<rows_, cols_>& left,
                                 const FixedSizeSerialDenseMatrix<rows_, cols_>& right);
  void Print(std::ostream& out) const;


  // operators
  /*
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator-() const;

  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator+(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const;
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator-(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const;
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator*(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const;
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator*(const double other) const;
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator/(const double other) const;

  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator+=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other);
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator-=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other);
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator*=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other);
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator*=(const double other);
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator/=(const double other);
  */

  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other);
  inline FixedSizeSerialDenseMatrix<rows_,cols_> operator=(const double other);

  inline bool operator==(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const;
  inline bool operator!=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const;

  inline double& operator()(unsigned int r, unsigned int c);
  inline const double& operator()(unsigned int r, unsigned int c) const;
  inline double& operator()(unsigned int r);  // for vectors, with check at compile-time
  inline const double& operator()(unsigned int r) const;
/*
    inline double* operator[](unsigned int c);
    inline const double* operator[](unsigned int c) const;
*/
};

template <unsigned int cols_, unsigned int rows_>
std::ostream& operator<<(std::ostream& out, const FixedSizeSerialDenseMatrix<rows_,cols_>& matrix);

// Constructors

template<unsigned int rows_, unsigned int cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::FixedSizeSerialDenseMatrix(bool setzero) : data_(0), isview_(false)
{
  DENSEERROR::Checker<(rows_ != 0) and (cols_ != 0)>::Matrix_dimensions_cannot_be_zero();
  data_ = new double[rows_*cols_];
  if (setzero)
    DENSEFUNCTIONS::clearMatrix<rows_,cols_>(data_);
}

template<unsigned int rows_, unsigned int cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::FixedSizeSerialDenseMatrix(double* d, bool view) :
    data_(0), isview_(view)
{
  DENSEERROR::Checker<(rows_ != 0) and cols_ != 0>::Matrix_dimensions_cannot_be_zero();
  if (isview_) {
    data_ = d;
  }
  else {
    data_ = new double[rows_*cols_];
    std::memcpy(data_, d, rows_*cols_*sizeof(double));
  }
}

template<unsigned int rows_, unsigned int cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::FixedSizeSerialDenseMatrix(FixedSizeSerialDenseMatrix<rows_,cols_>& source, bool view) :
    data_(0), isview_(view)
{
  DENSEERROR::Checker<(rows_ != 0) and cols_ != 0>::Matrix_dimensions_cannot_be_zero();
  if (isview_) {
    data_ = source.data_;
  }
  else {
    data_ = new double[rows_*cols_];
    std::memcpy(data_, source.data_, rows_*cols_*sizeof(double));
  }
}

template<unsigned int rows_, unsigned int cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::FixedSizeSerialDenseMatrix(const FixedSizeSerialDenseMatrix<rows_,cols_>& source) :
  data_(0), isview_(false)
{
  data_ = new double[rows_*cols_];
  std::memcpy(data_, source.data_, rows_*cols_*sizeof(double));
}

// Destructor
template<unsigned int rows_, unsigned int cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::~FixedSizeSerialDenseMatrix()
{
  if (not isview_)
    delete[] data_;
}

template<unsigned int rows_, unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_, cols_>::SetView(double* data)
{
  if (not isview_) {
    delete[] data_;
    isview_ = true;
  }
  data_ = data;
}

template<unsigned int rows_, unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_, cols_>::SetView(FixedSizeSerialDenseMatrix<rows_, cols_>& source)
{
  if (not isview_) {
    delete[] data_;
    isview_ = true;
  }
  data_ = source.data_;
}

template<unsigned int rows_, unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_, cols_>::SetCopy(double* data)
{
    if (isview_) {
        data_ = new double[rows_*cols_];
        isview_ = false;
    }
    std::memcpy(data_, data, rows_*cols_*sizeof(double));
}

template<unsigned int rows_, unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_, cols_>::SetCopy(FixedSizeSerialDenseMatrix<rows_, cols_>& source)
{
  if (isview_) {
    data_ = new double[rows_*cols_];
    isview_ = false;
  }
  std::memcpy(data_, source.data_, rows_*cols_*sizeof(double));
}

template<unsigned int rows_, unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_, cols_>::Copy(double* data)
{
    std::memcpy(data_, data, rows_*cols_*sizeof(double));
}

template<unsigned int rows_, unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_, cols_>::Copy(FixedSizeSerialDenseMatrix<rows_, cols_>& source)
{
  std::memcpy(data_, source.data_, rows_*cols_*sizeof(double));
}


// Determinant
template<unsigned int rows_, unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::Determinant() const
{
  DENSEERROR::Checker<rows_ == cols_>::Cannot_compute_determinant_of_nonsquare_matrix();
  return DENSEFUNCTIONS::determinant<rows_, cols_>(data_);
}

// Invert
template<unsigned int rows_, unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::Invert()
{
  DENSEERROR::Checker<rows_ == cols_>::Cannot_compute_inverse_of_nonsquare_matrix();
  return DENSEFUNCTIONS::invert<rows_, cols_>(data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::Invert(const FixedSizeSerialDenseMatrix<rows_,cols_>& other)
{
  DENSEERROR::Checker<rows_ == cols_>::Cannot_compute_inverse_of_nonsquare_matrix();
  return DENSEFUNCTIONS::invert<rows_, cols_>(data_, other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Clear()
{
  DENSEFUNCTIONS::clearMatrix<rows_,cols_>(data_);
}


template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::PutScalar(const double scalar)
{
  DENSEFUNCTIONS::putScalar<rows_, cols_>(scalar, data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::Dot(const FixedSizeSerialDenseMatrix<rows_, cols_>& other)
{
    return DENSEFUNCTIONS::dot<rows_, cols_>(data_, other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Abs(const FixedSizeSerialDenseMatrix<rows_, cols_>& other)
{
  return DENSEFUNCTIONS::abs<rows_, cols_>(data_, other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Reciprocal(const FixedSizeSerialDenseMatrix<rows_, cols_>& other)
{
  return DENSEFUNCTIONS::reciprocal<rows_, cols_>(data_, other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Scale(const double scalar)
{
    DENSEFUNCTIONS::scaleMatrix<rows_, cols_>(scalar, data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Scale(const double scalar,
                                                           const FixedSizeSerialDenseMatrix<rows_, cols_>& other)
{
    DENSEFUNCTIONS::scaleMatrix<rows_, cols_>(data_, scalar, other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Update(const FixedSizeSerialDenseMatrix<rows_, cols_>& other)
{
  DENSEFUNCTIONS::addMatrixTo<rows_,cols_>(data_,other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Update(const double scalarOther,
                                                            const FixedSizeSerialDenseMatrix<rows_, cols_>& other,
                                                            const double scalarThis = 1.0)
{
    DENSEFUNCTIONS::addMatrixTo<rows_, cols_>(scalarThis, data_, scalarOther, other.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Update(const double scalarLeft,
                                                            const FixedSizeSerialDenseMatrix<rows_, cols_>& left,
                                                            const double scalarRight,
                                                            const FixedSizeSerialDenseMatrix<rows_, cols_>& right,
                                                            const double scalarThis = 1.0)
{
    DENSEFUNCTIONS::addMatrix<rows_, cols_>(scalarThis, data_, scalarLeft, left.data_, scalarRight, right.data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::Norm1() const
{
  return DENSEFUNCTIONS::norm1<rows_, cols_>(data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::Norm2() const
{
  return DENSEFUNCTIONS::norm2<rows_, cols_>(data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::NormInf() const
{
  return DENSEFUNCTIONS::normInf<rows_, cols_>(data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::MinValue() const
{
  return DENSEFUNCTIONS::minValue<rows_, cols_>(data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::MaxValue() const
{
  return DENSEFUNCTIONS::maxValue<rows_, cols_>(data_);
}

template <unsigned int rows_,  unsigned int cols_>
inline double FixedSizeSerialDenseMatrix<rows_,cols_>::MeanValue() const
{
  return DENSEFUNCTIONS::meanValue<rows_, cols_>(data_);
}


template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Multiply(const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'N','N'>(data_,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyNN(const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'N','N'>(data_,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyNT(const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                       const FixedSizeSerialDenseMatrix<cols_, inner>& right)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'N','T'>(data_,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyTN(const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                       const FixedSizeSerialDenseMatrix<inner, cols_>& right)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'T','N'>(data_,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyTT(const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                       const FixedSizeSerialDenseMatrix<cols_, inner>& right)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'T','T'>(data_,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
  template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::Multiply(const double scalarOthers,
                                                              const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                                                              const FixedSizeSerialDenseMatrix<inner, cols_>& right,
                                                              const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'N','N'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyNN(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                                                                const FixedSizeSerialDenseMatrix<inner, cols_>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'N','N'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyNT(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<rows_, inner>& left,
                                                                const FixedSizeSerialDenseMatrix<cols_, inner>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'N','T'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyTN(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                                                                const FixedSizeSerialDenseMatrix<inner, cols_>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'T','N'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}

template <unsigned int rows_, unsigned int cols_>
template <unsigned int inner>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::MultiplyTT(const double scalarOthers,
                                                                const FixedSizeSerialDenseMatrix<inner, rows_>& left,
                                                                const FixedSizeSerialDenseMatrix<cols_, inner>& right,
                                                                const double scalarThis)
{
  DENSEFUNCTIONS::multiply<rows_,inner,cols_,'T','T'>(scalarThis,data_,scalarOthers,left.A(),right.A());
}


template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::ReciprocalMultiply(const double scalarOthers,
                                                                        const FixedSizeSerialDenseMatrix<rows_, cols_>& left,
                                                                        const FixedSizeSerialDenseMatrix<rows_, cols_>& right,
                                                                        const double scalarThis)
{
    DENSEFUNCTIONS::reciprocalMultiply<rows_, cols_>(scalarThis, data_, scalarOthers, left.A(), right.A());
}

template <unsigned int rows_,  unsigned int cols_>
inline void FixedSizeSerialDenseMatrix<rows_,cols_>::ReciprocalMultiply(const FixedSizeSerialDenseMatrix<rows_, cols_>& left,
                                                                        const FixedSizeSerialDenseMatrix<rows_, cols_>& right)
{
    DENSEFUNCTIONS::reciprocalMultiply<rows_, cols_>(data_, left.data_, right.data_);
}

//inline void FixedSizeSerialDenseMatrix<rows_,cols_>::ResetView(double* data);

template <unsigned int rows_,  unsigned int cols_>
void FixedSizeSerialDenseMatrix<rows_,cols_>::Print(std::ostream& out) const
{
  out << "FixedSizeSerialDenseMatrix<" << rows_ << ',' << cols_ << '>';
  if (isview_) out << " (view to memory only)";
  if (data_ == NULL) {
    out << " data is NULL!\n";
    return;
  }
  out << "\n[";
  for (unsigned int i = 0; i < rows_; ++i) {
    if (i!=0) out << ' ';
    for (unsigned int j = 0; j < cols_; ++j) {
      out << data_[i+rows_*j];
      if (j < cols_-1) out << ", ";
    }
    if (i < rows_-1) out << ",\n";
    else out << "]\n";
  }
}


template <unsigned int cols_, unsigned int rows_>
std::ostream& operator<<(std::ostream& out, const FixedSizeSerialDenseMatrix<rows_,cols_>& matrix)
{
    matrix.Print(out);
    return out;
}


template<unsigned int rows_, unsigned int cols_>
inline FixedSizeSerialDenseMatrix<rows_,cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::operator=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other)
{
  if (this == &other) return *this;
  if (isview_) {
    isview_ = false;
    data_ = new double[rows_*cols_];
  }
  DENSEFUNCTIONS::copyMatrix<rows_,cols_>(data_, other.data_);
  return *this;
}

template<unsigned int rows_, unsigned int cols_>
inline FixedSizeSerialDenseMatrix<rows_,cols_>
FixedSizeSerialDenseMatrix<rows_, cols_>::operator=(const double other)
{
  DENSEFUNCTIONS::putScalar<rows_,cols_>(other, data_);
}

template<unsigned int rows_, unsigned int cols_>
inline bool
FixedSizeSerialDenseMatrix<rows_, cols_>::operator==(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const
{
  if (data_ == other.data_) return true;
  return (std::memcmp(data_, other.data_, rows_*cols_*sizeof(double)) == 0);
}

template<unsigned int rows_, unsigned int cols_>
inline bool
FixedSizeSerialDenseMatrix<rows_, cols_>::operator!=(const FixedSizeSerialDenseMatrix<rows_,cols_>& other) const
{
  return not (*this == other);
}

// Access operator

template<unsigned int rows_, unsigned int cols_>
inline double& FixedSizeSerialDenseMatrix<rows_, cols_>::operator()(unsigned int r, unsigned int c)
{
  return data_[r+c*rows_];
}

template<unsigned int rows_, unsigned int cols_>
inline const double& FixedSizeSerialDenseMatrix<rows_, cols_>::operator()(unsigned int r, unsigned int c) const
{
  return data_[r+c*rows_];
}

template<unsigned int rows_, unsigned int cols_>
inline double& FixedSizeSerialDenseMatrix<rows_, cols_>::operator()(unsigned int r)
{
  DENSEERROR::Checker<(cols_==1) or (rows_==1)>::Cannot_call_1D_access_function_on_2D_matrix();
  return data_[r];
}

template<unsigned int rows_, unsigned int cols_>
inline const double& FixedSizeSerialDenseMatrix<rows_, cols_>::operator()(unsigned int r) const
{
    DENSEERROR::Checker<(cols_==1) or (rows_==1)>::Cannot_call_1D_access_function_on_2D_matrix();
    return data_[r];
}

/*
template<unsigned int rows_, unsigned int cols_>
inline double* FixedSizeSerialDenseMatrix<rows_, cols_>::operator[](unsigned int c)

{
  return data_+c*rows_;
}

template<unsigned int rows_, unsigned int cols_>
inline const double* FixedSizeSerialDenseMatrix<rows_, cols_>::operator[](unsigned int c) const
{
  return data_+c*rows_;
}
*/


/// Fixed size solver
template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
class FixedSizeSerialDenseSolver
{
private:

  FixedSizeSerialDenseSolver(const FixedSizeSerialDenseSolver<rows_, cols_,dim_rhs_>&);
  FixedSizeSerialDenseSolver& operator=(const FixedSizeSerialDenseSolver<rows_, cols_,dim_rhs_>&);

  void DeleteMatrices();

  Epetra_LAPACK lapack_;
  Epetra_BLAS blas_;

  FixedSizeSerialDenseMatrix<rows_,cols_>* matrix_;
  FixedSizeSerialDenseMatrix<rows_,cols_>* factor_;
  FixedSizeSerialDenseMatrix<cols_,dim_rhs_>* vec_X_;
  FixedSizeSerialDenseMatrix<rows_,dim_rhs_>* vec_B_;

  int* pivot_vec_;
  double* r_;
  double* c_;

  char trans_char_;
  bool equilibrate_;
  bool matEquilibrated_;
  bool rhsEquilibrated_;
  bool inverted_;
  bool solved_;
  bool factored_;
  bool transpose_;

public:
  FixedSizeSerialDenseSolver();
  ~FixedSizeSerialDenseSolver();

  void SetMatrix(FixedSizeSerialDenseMatrix<rows_,cols_>& mat);
  void SetVectors(FixedSizeSerialDenseMatrix<cols_,dim_rhs_>& X, FixedSizeSerialDenseMatrix<rows_,dim_rhs_>& B);
  void FactorWithEquilibration(bool b) { equilibrate_ = b; }
  void SolveWithTranspose(bool b) { transpose_ = b; trans_char_ = b ? 'T' : 'N';}

  int Factor();
  int Solve();
  int Invert();
  int ComputeEquilibrateScaling();
  int EquilibrateMatrix();
  int EquilibrateRHS();
  int UnequilibrateLHS();
};

template<unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::FixedSizeSerialDenseSolver() :
  lapack_(),
  blas_(),
  matrix_(NULL),
  factor_(NULL),
  vec_X_(NULL),
  vec_B_(NULL),
  pivot_vec_(NULL),
  r_(NULL),
  c_(NULL),
  trans_char_('N'),
  equilibrate_(false),
  matEquilibrated_(false),
  rhsEquilibrated_(false),
  inverted_(false),
  solved_(false),
  factored_(false),
  transpose_(false)
{
  DENSEERROR::Checker<rows_==cols_>::Matrix_size_in_solver_must_be_square();
}

template<unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::~FixedSizeSerialDenseSolver()
{
  if (factor_ and not (matrix_ == factor_)) delete factor_;
  if (pivot_vec_) delete[] pivot_vec_;
  if (r_) delete[] r_;
  if (c_) delete[] c_;
}

template<unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
void FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::DeleteMatrices()
{
  if (factor_ and not (matrix_ == factor_)) delete factor_;
  if (pivot_vec_) delete[] pivot_vec_;
  if (r_) delete[] r_;
  if (c_) delete[] c_;
}

template<unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
void FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::SetMatrix(FixedSizeSerialDenseMatrix<rows_,cols_>& mat)
{
  DeleteMatrices();
  c_ = NULL;
  r_ = NULL;
  pivot_vec_ = NULL;
  factor_ = NULL;
  matEquilibrated_ = rhsEquilibrated_ = inverted_ = factored_ = solved_ = false;
  vec_B_ = vec_X_ = NULL;
  matrix_ = &mat;
  factor_ = &mat;
}

template<unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
void FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::SetVectors(FixedSizeSerialDenseMatrix<cols_,dim_rhs_>& X,
                                                         FixedSizeSerialDenseMatrix<rows_,dim_rhs_>& B)
{
  rhsEquilibrated_ = false;
  vec_X_ = &X;
  vec_B_ = &B;
}

template<unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::Factor()
{
  int errnum = 0;
  if (equilibrate_) errnum = EquilibrateMatrix();
  if (errnum != 0) return errnum;
  if (pivot_vec_ == NULL) pivot_vec_ = new int[rows_ < cols_ ? rows_ : cols_];
  lapack_.GETRF(rows_, cols_, factor_->A(), rows_, pivot_vec_, &errnum);
  if (errnum != 0) return errnum;

  factored_ = true;
  return 0;
}

template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::Solve()
{
  int errnum = 0;
  if (equilibrate_) {
    errnum = EquilibrateRHS();
  }
  if (errnum != 0) return errnum;
  if (matEquilibrated_ != rhsEquilibrated_) return -1;
  // if not vec_B_ or not vec_X_

  if (inverted_) {
    if (vec_B_ == vec_X_) return -100;

    blas_.GEMM(trans_char_, 'N', cols_, dim_rhs_, cols_, 1.0, factor_->A(), rows_, vec_B_->A(), rows_, 0.0,
               vec_X_->A(), cols_);
    solved_ = true;
  }
  else {
    if (!factored_) {
      errnum = Factor();
      if (errnum != 0) return errnum;
    }

    if (vec_B_ != vec_X_) *vec_X_ = *vec_B_;
    lapack_.GETRS(trans_char_, cols_, dim_rhs_, factor_->A(), rows_, pivot_vec_, vec_X_->A(), cols_, &errnum);
    if (errnum != 0) return errnum;
    solved_ = true;
  }
  if (equilibrate_) errnum = UnequilibrateLHS();
  if (errnum != 0) return errnum;
  return 0;
}

template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::ComputeEquilibrateScaling()
{
  if (r_ != 0) return 0;
  int errnum;
  double rowcnd, colcnd, amax;
  r_ = new double [rows_];
  c_ = new double [cols_];
  lapack_.GEEQU(rows_, cols_, factor_->A(), rows_, r_, c_, &rowcnd, &colcnd, &amax, &errnum);
  if (errnum != 0) return errnum;

  //if (colcnd_<0.1 || rowcnd_<0.1 || amax_ < Epetra_Underflow || amax_ > Epetra_Overflow) shouldEquilibrate_ = true;
  return 0;
}

template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::EquilibrateMatrix()
{
  if (matEquilibrated_) return 0;
  int errnum = 0;
  if (r_ == 0) errnum = ComputeEquilibrateScaling();
  if (errnum != 0) return errnum;
  if (matrix_ == factor_) {
    double* ptr = matrix_->A();
    double s1;
    for (unsigned j = 0; j < cols_; ++j) {
      s1 = c_[j];
      for (unsigned i = 0; i < rows_; ++i) {
        *ptr *= s1*r_[i];
        ++ptr;
      }
    }
  }
  else {
    double* ptr1 = matrix_->A();
    double* ptr2 = factor_->A();
    double s1;
    for (unsigned j = 0; j < cols_; ++j) {
      s1 = c_[j];
      for (unsigned i = 0; i < rows_; ++i) {
        *ptr1 *= s1*r_[i];
        *ptr2 *= s1*r_[i];
        ++ptr1; ++ptr2;
      }
    }
  }

  matEquilibrated_ = true;

  return 0;
}

template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::EquilibrateRHS()
{
  if (rhsEquilibrated_) return 0;
  int errnum = 0;
  if (r_ == 0) errnum = ComputeEquilibrateScaling();
  if (errnum != 0) return errnum;
  double* r = transpose_ ? c_ : r_;
  double* ptr = vec_B_->A();
  for (unsigned j = 0; j < dim_rhs_; ++j) {
    for (unsigned i = 0; i < cols_; ++i) {
      *ptr *= r[i];
      ++ptr;
    }
  }

  rhsEquilibrated_ = true;
  return 0;
}

template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::UnequilibrateLHS()
{
  if (not rhsEquilibrated_) return 0;
  double* c = transpose_ ? r_ : c_;
  double* ptr = vec_X_->A();
  for (unsigned j = 0; j < dim_rhs_; ++j) {
    for (unsigned i = 0; i < rows_; ++i) {
      *ptr *= c[i];
      ++ptr;
    }
  }

  return 0;
}

template <unsigned int rows_, unsigned int cols_, unsigned int dim_rhs_>
int FixedSizeSerialDenseSolver<rows_,cols_,dim_rhs_>::Invert()
{
  int errnum;

  int lwork = 4*cols_;
  double* work = new double[lwork];
  lapack_.GETRI(cols_, factor_->A(), rows_, pivot_vec_, work, &lwork, &errnum);
  if (errnum != 0) return errnum;
  inverted_ = true;
  factored_ = false;

  return 0;
}

}

#endif // FIXED_SIZE_SERIAL_DENSE_MATRIX_H
