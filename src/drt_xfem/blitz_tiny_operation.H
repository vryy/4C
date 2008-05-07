#ifndef BLITZ_TINY_OPERATION_H
#define BLITZ_TINY_OPERATION_H

// to be moved to linalg in the future
namespace LINALG {
    /*!
     * \brief routines to work with epetra dense objects and variable + fixed sized blitz arrays
     */
    namespace DENSE { 
    
        //! matrix matrix^T product
        template<int ni, int nj, int nk, class M1, class M2, class M3>
        inline void MMt_product(const M1& deriv,const M2& xyze,M3& xjm)
        {
            for (int i = 0; i < ni; ++i)
            {
                for (int j = 0; j < nj; ++j)
                {
                    xjm(i,j) = 0.0;
                    for (int k = 0; k < nk; ++k)
                    {
                        xjm(i,j) += deriv(i,k)*xyze(j,k);
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
    
    }
}
#endif /*BLITZ_TINY_OPERATION_H_*/
