/*----------------------------------------------------------------------*/
/*!
\file FAD_utils.H

\brief A set of utility functions for Forward Automatic Differentiation

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/


#ifndef FAD_UTILS_H
#define FAD_UTILS_H

#include "../linalg/linalg_fixedsizematrix.H"
#include <Sacado.hpp>

typedef Sacado::Fad::DFad<double> FAD;

namespace FADUTILS
{
  /*!
  \brief Cast of FAD to double
  */
  inline double CastToDouble(FAD a)
  {
    return a.val();
  }

  /*!
  \brief Cast of double to double
  */
  inline double CastToDouble(double a)
  {
    return a;
  }

  /*!
  \brief Cast of a FAD matrix to a double matrix
  */
  template<typename type,int dim1, int dim2>
  LINALG::TMatrix<double,dim1,dim2> CastToDouble(LINALG::TMatrix<type,dim1,dim2> a)
  {
    LINALG::TMatrix<double,dim1,dim2> b(true);

    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<dim2;j++)
      {
        b(i,j)=CastToDouble(a(i,j));
      }
    }
    return b;
  }

  /*!
  \brief Calculate signum function of FAD or double quantity
  */
  template<typename type>
  type Signum(type a)
  {
    if (CastToDouble(a) >= 0.0)
      return 1.0;
    else
      return -1.0;
  }

  /*!
  \brief Calculate Norm of a scalar double quantity
  */
  inline double Norm(double a)
  {
    return std::sqrt(a*a);
  }

  /*!
  \brief Calculate Norm of a scalar FAD quantity
  */
  inline FAD Norm(FAD a)
  {
    return std::pow(a*a,0.5);
  }

  /*!
  \brief Calculate Norm of a double vector
  */
  template<int length>
  double VectorNorm(LINALG::TMatrix<double,length,1> v)
  {
    return v.Norm2();
  }

  /*!
  \brief Calculate Norm of a FAD vector
  */
  template<int length>
  FAD VectorNorm(LINALG::TMatrix<FAD,length,1> v)
  {
    FAD norm = 0.0;
    for (int i=0;i<length;i++)
    {
      norm+=v(i)*v(i);
    }
    norm = std::pow(norm,0.5);
    return norm;
  }

  //! Calculates the Norm of a FAD vector, since .Norm2() is not available for FAD vectors
  template<typename T>
  T Norm(LINALG::TMatrix<T,3,1> v)
  {
    T norm = 0.0;
    for (int i=0;i<3;i++)
    {
      norm+=v(i)*v(i);
    }
    norm = std::sqrt(norm);

    return norm;
  }

  /*!
  \brief Calculate inner product of two FAD or double vectors
  */
  template<typename type>
  type ScalarProduct(LINALG::TMatrix<type,3,1> a, LINALG::TMatrix<type,3,1> b)
  {
    return a(0)*b(0)+a(1)*b(1)+a(2)*b(2);
  }

  /*!
  \brief Calculate difference of two FAD or double vectors
  */
  template<typename type>
  LINALG::TMatrix<type,3,1> DiffVector(LINALG::TMatrix<type,3,1> a, LINALG::TMatrix<type,3,1> b)
  {
    LINALG::TMatrix<type,3,1> c(true);
    for (int i=0;i<3;i++)
      c(i)=a(i)-b(i);

    return c;
  }

  /*!
  \brief Calculate vector product of two FAD or double vectors
  */
  template<typename type>
  LINALG::TMatrix<type,3,1> VectorProduct(LINALG::TMatrix<type,3,1> first_vector, LINALG::TMatrix<type,3,1> second_vector)
  {
    LINALG::TMatrix<type,3,1> result_vector;
    result_vector.Clear();
    LINALG::TMatrix<type,3,3> S_first_vector;
    S_first_vector.Clear();

    S_first_vector(0,0) = 0.0;
    S_first_vector(0,1) = -first_vector(2);
    S_first_vector(0,2) = first_vector(1);
    S_first_vector(1,0) = first_vector(2);
    S_first_vector(1,1) = 0.0;
    S_first_vector(1,2) = -first_vector(0);
    S_first_vector(2,0) = -first_vector(1);
    S_first_vector(2,1) = first_vector(0);
    S_first_vector(2,2) = 0.0;

    result_vector.Multiply(S_first_vector, second_vector);

    return result_vector;
  }

} // FADUTILS

#endif  // #ifndef FAD_UTILS_H
