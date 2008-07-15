/*----------------------------------------------------------------------*/
/*!
\file blitz_tiny_operation.H

\brief routines for the tiny blitz arrays and matrix vector routines to combine various array formats

\note: such mixing is not optimal in the sense, that matrix vector routines are tuned for one array library and
it is hard to write similarily performing code. So this is more a debug solution until one settles for a array type
in ones personal part of BACI 

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
/*----------------------------------------------------------------------*/
#ifndef BLITZ_TINY_OPERATION_H
#define BLITZ_TINY_OPERATION_H

namespace XFEM
{
  
  /*!
   * \brief routines to work with epetra dense objects and variable + fixed sized blitz arrays
   * 
   * These routines are primarily intended for debugging in the XFEM and XFluid3 code.
   * Their performance for various array sizes is not garantiert.
   */
  namespace BLITZTINY
  {    
    
    //! matrix matrix^T product
    template<int ni, int nj, int nk, class M1, class M2, class M3>
    inline void MMt_product(const M1& m1, const M2& m2, M3& m_out)
    {
      for (int i = 0; i < ni; ++i)
      {
        for (int j = 0; j < nj; ++j)
        {
          m_out(i,j) = 0.0;
          for (int k = 0; k < nk; ++k)
          {
            m_out(i,j) += m1(i,k)*m2(j,k);
          }
        }
      }
    }
    //! matrix^T matrix product
    template<int ni, int nj, int nk, class M1, class M2, class M3>
    inline void MtM_product(const M1& m1, const M2& m2, M3& m_out)
    {
      for (int i = 0; i < ni; ++i)
      {
        for (int j = 0; j < nj; ++j)
        {
          m_out(i,j) = 0.0;
          for (int k = 0; k < nk; ++k)
          {
            m_out(i,j) += m1(k,i)*m2(k,j);
          }
        }
      }
    }
    //! matrix matrix product
    template<int ni, int nj, int nk, class M1, class M2, class M3>
    inline void MM_product(const M1& m1, const M2& m2, M3& m_out)
    {
      for (int i = 0; i < ni; ++i)
      {
        for (int j = 0; j < nj; ++j)
        {
          m_out(i,j) = 0.0;
          for (int k = 0; k < nk; ++k)
          {
            m_out(i,j) += m1(i,k)*m2(k,j);
          }
        }
      }
    }
    //! matrix vector product
    template<int ni, int nj, class M1, class V1, class V2>
    inline void MV_product(const M1& m, const V1& v, V2& rv)
    {
      for (int i = 0; i < ni; ++i)
      {
        rv(i) = 0.0;
        for (int j = 0; j < nj; ++j)
        {
          rv(i) += m(i,j)*v(j);
        }
      }
    }
    //! matrix scaling with scalar f
    template<int ni, int nj, class M1, class F1>
    inline void M_scale(M1& m, const F1& f)
    {
      for (int i = 0; i < ni; ++i)
      {
        for (int j = 0; j < nj; ++j)
        {
          m(i,j) *= f;
        }
      }
    }
    //! vector scaling with scalar f
    template<int ni, class V1, class F1>
    inline void V_scale(V1& v, const F1& f)
    {
      for (int i = 0; i < ni; ++i)
      {
          v(i) *= f;
      }
    }
    //! matrix^T vector product
    template<int ni, int nj, class M1, class V1, class V2>
    inline void MtV_product(const M1& m, const V1& v, V2& rv)
    {
      for (int i = 0; i < ni; ++i)
      {
        rv(i) = 0.0;
        for (int j = 0; j < nj; ++j)
        {
          rv(i) += m(j,i)*v(j);
        }
      }
    }    
  }
}
#endif /*BLITZ_TINY_OPERATION_H_*/
