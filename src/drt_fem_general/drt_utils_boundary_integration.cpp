/*!----------------------------------------------------------------------
\file drt_utils_boundary_integration.cpp

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#include "drt_utils_boundary_integration.H"
#include "drt_utils_integration.H"


/* compute kovariant metric tensor G for surface element     gammi 04/07

                        +-       -+
                        | g11 g12 |
                    G = |         |
                        | g12 g22 |
                        +-       -+

 where (o denotes the inner product, xyz a vector)


                            dxyz   dxyz
                    g11 =   ---- o ----
                             dr     dr

                            dxyz   dxyz
                    g12 =   ---- o ----
                             dr     ds

                            dxyz   dxyz
                    g22 =   ---- o ----
                             ds     ds


 and the square root of the first fundamental form


                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/

 they are needed for the integration over the surface element

*/
void  DRT::UTILS::ComputeMetricTensorForSurface(
  const Epetra_SerialDenseMatrix& xyze,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix&       metrictensor,
  double                         *sqrtdetg)
{
  /*
  |                                              0 1 2
  |                                             +-+-+-+
  |       0 1 2              0...iel-1          | | | | 0
  |      +-+-+-+             +-+-+-+-+          +-+-+-+
  |      | | | | 1           | | | | | 0        | | | | .
  |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
  |      | | | | 2           | | | | | 1        | | | | .
  |      +-+-+-+             +-+-+-+-+          +-+-+-+
  |                                             | | | | iel-1
  |                                             +-+-+-+
  |
  |       dxyzdrs             deriv              xyze^T
  |
  |
  |                                 +-            -+
  |                                 | dx   dy   dz |
  |                                 | --   --   -- |
  |                                 | dr   dr   dr |
  |     yields           dxyzdrs =  |              |
  |                                 | dx   dy   dz |
  |                                 | --   --   -- |
  |                                 | ds   ds   ds |
  |                                 +-            -+
  |
  */
  Epetra_SerialDenseMatrix dxyzdrs (2,3);

  dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);

  /*
  |
  |      +-           -+    +-            -+   +-            -+ T
  |      |             |    | dx   dy   dz |   | dx   dy   dz |
  |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
  |      |             |    | dr   dr   dr |   | dr   dr   dr |
  |      |             |  = |              | * |              |
  |      |             |    | dx   dy   dz |   | dx   dy   dz |
  |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
  |      |             |    | ds   ds   ds |   | ds   ds   ds |
  |      +-           -+    +-            -+   +-            -+
  |
  | the calculation of g21 is redundant since g21=g12
  */
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);

/*
                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/
*/

  sqrtdetg[0]= sqrt(metrictensor(0,0)*metrictensor(1,1)
                    -
                    metrictensor(0,1)*metrictensor(1,0));

  return;
}



/*-----------------------------------------------------------------

\brief Transform Gausspoints on surface element to 3d space of
       parent element (required for integrations of parent-element
       shape functions over boundary elements, required for example
       in weak dirichlet boundary conditions).

  -----------------------------------------------------------------*/
void DRT::UTILS::SurfaceGPToParentGP(
  Epetra_SerialDenseMatrix              & pqxg     ,
  Epetra_SerialDenseMatrix              & derivtrafo,
  const DRT::UTILS::IntegrationPoints2D & intpoints,
  const DRT::Element::DiscretizationType  pdistype ,
  const DRT::Element::DiscretizationType  distype  ,
  const int                               surfaceid)
{
  // resize output array
  pqxg.Shape(intpoints.nquad,3);
  derivtrafo.Shape(3,3);

  if(distype==DRT::Element::quad4 && pdistype==DRT::Element::hex8)
  {
    switch(surfaceid)
    {
    case 0:
    {
      // t=-1
      /*
                parent               surface

                 r|                    s|
                  |                     |
             1         2           3         2
              +-------+             +-------+
              |   |   |      s      |   |   |      r
              |   +-- |  -----      |   +-- |  -----
              |       |             |       |
              +-------+             +-------+
             0         3           0         1
      */

      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][1];
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
        pqxg(iquad,2)=-1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=0.0;
      break;
    }
    case 1:
    {
      // s=-1
      /*
                parent               surface

                 t|                    s|
                  |                     |
             4         5           3         2
              +-------+             +-------+
              |   |   |      r      |   |   |      r
              |   +-- |  -----      |   +-- |  -----
              |       |             |       |
              +-------+             +-------+
             0         1           0         1
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)=-1.0;
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 2:
    {
      // r= 1
      /*
                parent               surface

                 t|                    s|
                  |                     |
             5         6           3         2
              +-------+             +-------+
              |   |   |      s      |   |   |      r
              |   +-- |  -----      |   +-- |  -----
              |       |             |       |
              +-------+             +-------+
             1         2           0         1
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= 1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 3:
    {
      // s= 1
      /*
                parent               surface

                 t|                    s|
                  |                     |
             6         7           3         2
              +-------+             +-------+
        r     |   |   |             |   |   |      r
        ----  |   +-- |             |   +-- |  -----
              |       |             |       |
              +-------+             +-------+
             2         3           0         1
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)=-intpoints.qxg[iquad][0];
        pqxg(iquad,1)= 1.0;
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=-1.0;
      derivtrafo(1,1)=0.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 4:
    {
      // r=-1
      /*
                parent               surface

                 s|                    s|
                  |                     |
             3         7           3         2
              +-------+             +-------+
              |   |   |      t      |   |   |      r
              |   +-- |  -----      |   +-- |  -----
              |       |             |       |
              +-------+             +-------+
             0         4           0         1
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)=-1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][1];
        pqxg(iquad,2)= intpoints.qxg[iquad][0];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 5:
    {
      // t=1
      /*
                parent               surface

                 s|                    s|
                  |                     |
             7         6           3         2
              +-------+             +-------+
              |   |   |      r      |   |   |      r
              |   +-- |  -----      |   +-- |  -----
              |       |             |       |
              +-------+             +-------+
             4         5           0         1
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)= intpoints.qxg[iquad][1];
        pqxg(iquad,2)= 1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=0.0;
      break;
    }
    default:
      dserror("invalid number of surfaces, unable to determine intpoint in parent");
    }
  }
  else if(distype==DRT::Element::nurbs9 && pdistype==DRT::Element::nurbs27)
  {
    switch(surfaceid)
    {
    case 0:
    {
      // t=-1
      /*
                parent               surface

                 s|                    s|
                  |                     |
              +---+---+             +---+---+
             6|  7|  8|      r     6|  7|  8|      r
              +   +-- +  -----      +   +-- +  -----
             3|  4   5|            3|  4   5|
              +---+---+             +---+---+
             0   1   2             0   1   2
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)= intpoints.qxg[iquad][1];
        pqxg(iquad,2)=-1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=0.0;
      break;
    }
    case 1:
    {
      // t=+1
      /*
                parent               surface

                 s|                    s|
                  |                     |
              +---+---+             +---+---+
            24| 25| 26|      r     6|  7|  8|      r
              +   +-- +  -----      +   +-- +  -----
            21| 22  23|            3|  4   5|
              +---+---+             +---+---+
            18  19  20             0   1   2
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)= intpoints.qxg[iquad][1];
        pqxg(iquad,2)= 1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=0.0;
      break;
    }
    case 2:
    {
      // s=-1
      /*
                parent               surface

                 t|                    s|
                  |                     |
              +---+---+             +---+---+
            18| 19| 20|      r     6|  7|  8|      r
              +   +-- +  -----      +   +-- +  -----
             9| 10  11|            3|  4   5|
              +---+---+             +---+---+
             0   1   2             0   1   2
      */

      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)=-1.0;
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 3:
    {
      // s=+1
      /*
                parent               surface

                 t|                    s|
                  |                     |
              +---+---+             +---+---+
            24| 25| 26|    r       6|  7|  8|      r
              +   +-- + ----        +   +-- +  -----
            15| 16  17|            3|  4   5|
              +---+---+             +---+---+
             6   7   8             0   1   2
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)= 1.0;
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 4:
    {
      // r=+1
      /*
                parent               surface

                 t|                    s|
                  |                     |
              +---+---+             +---+---+
            20| 23| 26|      s     6|  7|  8|      r
              +   +-- +  -----      +   +-- +  -----
            11| 14  17|            3|  4   5|
              +---+---+             +---+---+
             2   5   8             0   1   2
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= 1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    case 5:
    {
      // r=-1
      /*
                parent               surface

                 t|                    s|
                  |                     |
              +---+---+             +---+---+
            18| 21| 24|      s     6|  7|  8|      r
              +   +-- +  -----      +   +-- +  -----
             9| 12  15|            3|  4   5|
              +---+---+             +---+---+
             0   3   6             0   1   2
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)=-1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
        pqxg(iquad,2)= intpoints.qxg[iquad][1];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=1.0;
      derivtrafo(2,2)=1.0;
      break;
    }
    default:
      dserror("invalid number of surfaces, unable to determine intpoint in parent");
    }
  }
  else
  {
      dserror("only quad4/hex8 and nurbs9/nurbs27 mappings of surface gausspoint to parent element implemented up to now\n");
  }

  return;
}


/*-----------------------------------------------------------------

\brief Transform Gausspoints on line element to 2d space of
       parent element (required for integrations of parent-element
       shape functions over boundary elements, for example
       in weak dirichlet boundary conditions).

  -----------------------------------------------------------------*/
void DRT::UTILS::LineGPToParentGP(
  Epetra_SerialDenseMatrix              & pqxg     ,
  Epetra_SerialDenseMatrix              & derivtrafo,
  const DRT::UTILS::IntegrationPoints1D & intpoints,
  const DRT::Element::DiscretizationType  pdistype ,
  const DRT::Element::DiscretizationType  distype  ,
  const int                               lineid   )
{
  // resize output array
  pqxg.Shape(intpoints.nquad,2);
  derivtrafo.Shape(2,2);

  if(distype==DRT::Element::line2 && pdistype==DRT::Element::quad4)
  {
    switch(lineid)
    {
    case 0:
    {
    /*                s|
                       |

             3                   2
              +-----------------+
              |                 |
              |                 |
              |                 |
              |        |        |             r
              |        +--      |         -----
              |                 |
              |                 |
              |                 |
              |                 |
              +-----------*-----+
             0                   1
                    -->|gp|<--               */


      // s=-1
      /*

                parent                line

                             r                     r
              +---+---+  -----      +---+---+ ------
             0   1   2             0   1   2

      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)=-1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      break;
    }
    case 1:
    {
    /*                s|
                       |

             3                   2
              +-----------------+
              |                 | |
              |                 | v
              |                 *---
              |        |        | gp          r
              |        +--      |---      -----
              |                 | ^
              |                 | |
              |                 |
              |                 |
              +-----------------+
             0                   1
                                             */

      // r=+1
      /*
                parent               surface

                 s|                        r|
                  |                         |
                      +                     +
                     8|                    2|
                      +                     +
                     5|                    1|
                      +                     +
                     2                     0
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= 1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      break;
    }
    case 2:
    {
    /*                s|
                       |

             3   -->|gp|<--
              +-----*-----------+
              |                 |
              |                 |
              |                 |
              |        |        |             r
              |        +--      |         -----
              |                 |
              |                 |
              |                 |
              |                 |
              +-----------------+
             0                   1
                                             */

      // s=+1
      /*

                parent                line

                             r                           r
              +---+---+  -----             +---+---+ -----
             6   7   8                    0   1   2

      */

      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)=-intpoints.qxg[iquad][0];
        pqxg(iquad,1)= 1.0;
      }
      derivtrafo(0,0)=-1.0;
      derivtrafo(1,1)=0.0;
      break;
    }
    case 3:
    {
    /*                s|
                       |

             3
              +-----*-----------+
              |                 |
              |                 |
            | |                 |
            v |        |        |             r
           ---|        +--      |         -----
            gp|                 |
           ---*                 |
            ^ |                 |
            | |                 |
              +-----------------+
             0                   1
                                             */

      // r=-1
      /*
                parent               surface

                 s|                        r|
                  |                         |
               +                            +
	      6|                           2|
               +                            +
              3|                           1|
               +                            +
              0                            0
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)=-1.0;
        pqxg(iquad,1)=-intpoints.qxg[iquad][0];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=-1.0;
      break;
    }
    default:
      dserror("invalid number of lines, unable to determine intpoint in parent");
    }

  }
  else if(distype==DRT::Element::nurbs3 && pdistype==DRT::Element::nurbs9)
  {
    switch(lineid)
    {
    case 0:
    {
      // s=-1
      /*

                parent                line

                             r                     r
              +---+---+  -----      +---+---+ ------
             0   1   2             0   1   2

      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)=-1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      break;
    }
    case 1:
    {
      // r=+1
      /*
                parent               surface

                 s|                        r|
                  |                         |
                      +                     +
                     8|                    2|
                      +                     +
                     5|                    1|
                      +                     +
                     2                     0
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= 1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=1.0;
      break;
    }
    case 2:
    {
      // s=+1
      /*

                parent                line

                             r                           r
              +---+---+  -----             +---+---+ -----
             6   7   8                    0   1   2

      */

      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)= intpoints.qxg[iquad][0];
        pqxg(iquad,1)= 1.0;
      }
      derivtrafo(0,0)=1.0;
      derivtrafo(1,1)=0.0;
      break;
    }
    case 3:
    {
      // r=-1
      /*
                parent               surface

                 s|                        r|
                  |                         |
               +                            +
              6|                           2|
               +                            +
              3|                           1|
               +                            +
              0                            0
      */
      for (int iquad=0;iquad<intpoints.nquad;++iquad)
      {
        pqxg(iquad,0)=-1.0;
        pqxg(iquad,1)= intpoints.qxg[iquad][0];
      }
      derivtrafo(0,0)=0.0;
      derivtrafo(1,1)=1.0;
      break;
    }
    default:
      dserror("invalid number of lines, unable to determine intpoint in parent");
    }
  }
  else
  {
      dserror("only line2/quad4 and nurbs3/nurbs9 mappings of surface gausspoint to parent element implemented up to now\n");
  }

  return;
}

/*-----------------------------------------------------------------
\brief Template version of transformation of Gausspoints on boundary element to space of
       parent element
 ----------------------------------------------------------------------------------*/

//! specialization for 3D
template<>
void DRT::UTILS::BoundaryGPToParentGP<3>(
Epetra_SerialDenseMatrix                                                    & pqxg     ,
Epetra_SerialDenseMatrix                                                    & derivtrafo,
const DRT::UTILS::IntPointsAndWeights<2> &                                    intpoints,
const DRT::Element::DiscretizationType                                        pdistype ,
const DRT::Element::DiscretizationType                                        distype  ,
const int                                                                     surfaceid)
{
  DRT::UTILS::SurfaceGPToParentGP(
    pqxg     ,
    derivtrafo,
    intpoints.IP(),
    pdistype ,
    distype  ,
    surfaceid);
  return;
}

//! specialization for 2D
template<>
void DRT::UTILS::BoundaryGPToParentGP<2>(
Epetra_SerialDenseMatrix                                                    & pqxg     ,
Epetra_SerialDenseMatrix                                                    & derivtrafo,
const DRT::UTILS::IntPointsAndWeights<1> &                                    intpoints,
const DRT::Element::DiscretizationType                                        pdistype ,
const DRT::Element::DiscretizationType                                        distype  ,
const int                                                                     surfaceid)
{
  DRT::UTILS::LineGPToParentGP(
    pqxg     ,
    derivtrafo,
    intpoints.IP(),
    pdistype ,
    distype  ,
    surfaceid);
  return;
}
