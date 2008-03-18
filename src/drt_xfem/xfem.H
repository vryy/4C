/*!
\file xfem.H

\brief definitions used in xfem context

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFEM_H
#define XFEM_H

#include <blitz/array.h>

//! namespace for all xfem related stuff
namespace XFEM
{

    typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
    typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix
    
    typedef BlitzVec PosX;             ///< position in physical coordinates (x,y,z)
    typedef BlitzVec PosXiDomain;      ///< position in element coordinates of domain element (r,s,t)
    typedef BlitzVec PosXiBoundary;    ///< position in element coordinates of boundary element (r,s)
    typedef BlitzVec PosEtaDomain;     ///< position in local coordinates of a domain integrationcell (r,s,t)
    typedef BlitzVec PosEtaBoundary;   ///< position in local coordinates of a boundary integrationcell (r,s)

    enum IntegrationType
    {
        sub_integration,
        standard_integration
    };
    
    inline BlitzVec toBlitzArray(const double* x)
    {
        static BlitzVec blitz_x(3);
        blitz_x(0) = x[0];
        blitz_x(1) = x[1];
        blitz_x(2) = x[2];
        
        // doesn't work... why?
        // BlitzVec blitz_x(x, blitz::shape(3), blitz::duplicateData);
        return blitz_x;
    }

}

#endif  // #ifndef XFEM_H
#endif  // #ifdef CCADISCRET
