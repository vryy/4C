/*----------------------------------------------------------------------*/
/*!
\file xfluid3_stationary.H

\brief Internal implementation of Fluid3 element -- stationary formulation

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef XFLUID3_STATIONARY_H
#define XFLUID3_STATIONARY_H

#include <math.h>
#include "xfluid3.H"
#include <blitz/array.h>
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/xfem.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_xfem/coordinate_transformation.H"
#include "../drt_lib/drt_utils.H"
#include "fluid3_stabilization.H"
#include "xfluid3_local_assembler.H"
#include "xfluid3_interpolation.H"

#include <Epetra_SerialDenseSolver.h>

using namespace std;
using namespace XFEM::PHYSICS;


typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix


namespace XFLUID
{

/*!
 * \brief calculate the body force from the nodes condition
 */
/*----------------------------------------------------------------------*
 |  get the body force in the nodes of the element (private) gammi 04/07|
 |  the Neumann condition associated with the nodes is stored in the    |
 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType DISTYPE>
BlitzMat BodyForce(DRT::Node** nodes, const double pseudotime)
{
  dserror("not adapted to xfem (, yet)!!!");
    
  std::vector<DRT::Condition*> myneumcond;
  const int numnode = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<numnode;inode++)
  {
    nodes[inode]->GetCondition("VolumeNeumann",myneumcond);

    if (myneumcond.size()>1)
    {
      dserror("more than one VolumeNeumann cond on one node");
    }
    if (myneumcond.size()==1)
    {
      nodecount++;
    }
  }

  BlitzMat edeadng_(3,numnode,blitz::ColumnMajorArray<2>());
  if (nodecount == numnode)
  {
    // find out whether we will use a (pseudo-)time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a (pseudo-)timecurve
    {
      // factor for the intermediate step
      if(pseudotime >= 0.0)
      {
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(pseudotime);
      }
      else
      {
    // do not compute an "alternative" curvefac here since a negative pseudotime value
    // indicates an error.
        dserror("Negative pseudotime value in body force calculation: time = %f",pseudotime);
        //curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
      }
    }
    else // we do not have a (pseudo-)timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // set this condition to the edeadng array
    for (int jnode=0; jnode<numnode; jnode++)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<3;isd++)
      {
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_ = 0.;
  }
  
  return edeadng_;
}

  /*!
   * \brief calculate second global derivatives w.r.t x/y/z at point r,s,t
   */
  /*----------------------------------------------------------------------*
   |  calculate second global derivatives w.r.t. x,y,z at point r,s,t
   |                                            (private)      gammi 07/07
   |
   | From the six equations
   |
   |              +-                     -+
   |  d^2N     d  | dx dN   dy dN   dz dN |
   |  ----   = -- | --*-- + --*-- + --*-- |
   |  dr^2     dr | dr dx   dr dy   dr dz |
   |              +-                     -+
   |
   |              +-                     -+
   |  d^2N     d  | dx dN   dy dN   dz dN |
   |  ------ = -- | --*-- + --*-- + --*-- |
   |  ds^2     ds | ds dx   ds dy   ds dz |
   |              +-                     -+
   |
   |              +-                     -+
   |  d^2N     d  | dx dN   dy dN   dz dN |
   |  ----   = -- | --*-- + --*-- + --*-- |
   |  dt^2     dt | dt dx   dt dy   dt dz |
   |              +-                     -+
   |
   |              +-                     -+
   |  d^2N     d  | dx dN   dy dN   dz dN |
   | -----   = -- | --*-- + --*-- + --*-- |
   | ds dr     ds | dr dx   dr dy   dr dz |
   |              +-                     -+
   |
   |              +-                     -+
   |  d^2N     d  | dx dN   dy dN   dz dN |
   | -----   = -- | --*-- + --*-- + --*-- |
   | dt dr     dt | dr dx   dr dy   dr dz |
   |              +-                     -+
   |
   |              +-                     -+
   |  d^2N     d  | dx dN   dy dN   dz dN |
   | -----   = -- | --*-- + --*-- + --*-- |
   | ds dt     ds | dt dx   dt dy   dt dz |
   |              +-                     -+
   |
   | the matrix (jacobian-bar matrix) system
   |
   | +-                                                                                         -+   +-    -+
   | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
   | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
   | |   \dr/          \dr/           \dr/             dr dr           dr dr           dr dr     |   | dx^2 |
   | |                                                                                           |   |      |
   | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
   | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
   | |   \ds/          \ds/           \ds/             ds ds           ds ds           ds ds     |   | dy^2 |
   | |                                                                                           |   |      |
   | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
   | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
   | |   \dt/          \dt/           \dt/             dt dt           dt dt           dt dt     |   | dz^2 |
   | |                                                                                           | * |      |
   | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
   | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
   | |   dr ds         dr ds          dr ds        dr ds   ds dr   dr ds   ds dr  dr ds   ds dr  |   | dxdy |
   | |                                                                                           |   |      |
   | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
   | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
   | |   dr dt         dr dt          dr dt        dr dt   dt dr   dr dt   dt dr  dr dt   dt dr  |   | dxdz |
   | |                                                                                           |   |      |
   | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
   | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
   | |   dt ds         dt ds          dt ds        dt ds   ds dt   dt ds   ds dt  dt ds   ds dt  |   | dydz |
   | +-                                                                                         -+   +-    -+
   |
   |                  +-    -+     +-                           -+
   |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
   |                  | ---- |     | ----*-- + ----*-- + ----*-- |
   |                  | dr^2 |     | dr^2 dx   dr^2 dy   dr^2 dz |
   |                  |      |     |                             |
   |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
   |                  | ---- |     | ----*-- + ----*-- + ----*-- |
   |                  | ds^2 |     | ds^2 dx   ds^2 dy   ds^2 dz |
   |                  |      |     |                             |
   |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
   |                  | ---- |     | ----*-- + ----*-- + ----*-- |
   |                  | dt^2 |     | dt^2 dx   dt^2 dy   dt^2 dz |
   |              =   |      |  -  |                             |
   |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
   |                  | ---- |     | ----*-- + ----*-- + ----*-- |
   |                  | drds |     | drds dx   drds dy   drds dz |
   |                  |      |     |                             |
   |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
   |                  | ---- |     | ----*-- + ----*-- + ----*-- |
   |                  | drdt |     | drdt dx   drdt dy   drdt dz |
   |                  |      |     |                             |
   |                  | d^2N |     | d^2x dN   d^2y dN   d^2z dN |
   |                  | ---- |     | ----*-- + ----*-- + ----*-- |
   |                  | dtds |     | dtds dx   dtds dy   dtds dz |
   |                  +-    -+     +-                           -+
   |
   |
   | is derived. This is solved for the unknown global derivatives.
   |
   |
   |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
   |                                              |           |
   |                                              +-----------+
   |                                              'chainrulerhs'
   |                                     |                    |
   |                                     +--------------------+
   |                                          'chainrulerhs'
   |
   *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE>
BlitzMat gder2(
        const BlitzMat& xjm_,
        const BlitzMat& derxy_,
        const BlitzMat& deriv2_,
        const BlitzMat& xyze_)
{
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index

    // initialize and zero out everything
    Epetra_SerialDenseMatrix bm(6,6);

    // calculate elements of jacobian_bar matrix
    bm(0,0) = xjm_(0,0)*xjm_(0,0);
    bm(1,0) = xjm_(1,0)*xjm_(1,0);
    bm(2,0) = xjm_(2,0)*xjm_(2,0);
    bm(3,0) = xjm_(0,0)*xjm_(1,0);
    bm(4,0) = xjm_(0,0)*xjm_(2,0);
    bm(5,0) = xjm_(2,0)*xjm_(1,0);

    bm(0,1) = xjm_(0,1)*xjm_(0,1);
    bm(1,1) = xjm_(1,1)*xjm_(1,1);
    bm(2,1) = xjm_(2,1)*xjm_(2,1);
    bm(3,1) = xjm_(0,1)*xjm_(1,1);
    bm(4,1) = xjm_(0,1)*xjm_(2,1);
    bm(5,1) = xjm_(2,1)*xjm_(1,1);

    bm(0,2) = xjm_(0,2)*xjm_(0,2);
    bm(1,2) = xjm_(1,2)*xjm_(1,2);
    bm(2,2) = xjm_(2,2)*xjm_(2,2);
    bm(3,2) = xjm_(0,2)*xjm_(1,2);
    bm(4,2) = xjm_(0,2)*xjm_(2,2);
    bm(5,2) = xjm_(2,2)*xjm_(1,2);

    bm(0,3) = 2.*xjm_(0,0)*xjm_(0,1);
    bm(1,3) = 2.*xjm_(1,0)*xjm_(1,1);
    bm(2,3) = 2.*xjm_(2,0)*xjm_(2,1);
    bm(3,3) = xjm_(0,0)*xjm_(1,1)+xjm_(1,0)*xjm_(0,1);
    bm(4,3) = xjm_(0,0)*xjm_(2,1)+xjm_(2,0)*xjm_(0,1);
    bm(5,3) = xjm_(1,0)*xjm_(2,1)+xjm_(2,0)*xjm_(1,1);

    bm(0,4) = 2.*xjm_(0,0)*xjm_(0,2);
    bm(1,4) = 2.*xjm_(1,0)*xjm_(1,2);
    bm(2,4) = 2.*xjm_(2,0)*xjm_(2,2);
    bm(3,4) = xjm_(0,0)*xjm_(1,2)+xjm_(1,0)*xjm_(0,2);
    bm(4,4) = xjm_(0,0)*xjm_(2,2)+xjm_(2,0)*xjm_(0,2);
    bm(5,4) = xjm_(1,0)*xjm_(2,2)+xjm_(2,0)*xjm_(1,2);

    bm(0,5) = 2.*xjm_(0,1)*xjm_(0,2);
    bm(1,5) = 2.*xjm_(1,1)*xjm_(1,2);
    bm(2,5) = 2.*xjm_(2,1)*xjm_(2,2);
    bm(3,5) = xjm_(0,1)*xjm_(1,2)+xjm_(1,1)*xjm_(0,2);
    bm(4,5) = xjm_(0,1)*xjm_(2,2)+xjm_(2,1)*xjm_(0,2);
    bm(5,5) = xjm_(1,1)*xjm_(2,2)+xjm_(2,1)*xjm_(1,2);

    /*------------------ determine 2nd derivatives of coord.-functions */

    /*
    |
    |         0 1 2              0...iel-1
    |        +-+-+-+             +-+-+-+-+        0 1 2
    |        | | | | 0           | | | | | 0     +-+-+-+
    |        +-+-+-+             +-+-+-+-+       | | | | 0
    |        | | | | 1           | | | | | 1   * +-+-+-+ .
    |        +-+-+-+             +-+-+-+-+       | | | | .
    |        | | | | 2           | | | | | 2     +-+-+-+
    |        +-+-+-+       =     +-+-+-+-+       | | | | .
    |        | | | | 3           | | | | | 3     +-+-+-+ .
    |        +-+-+-+             +-+-+-+-+       | | | | .
    |        | | | | 4           | | | | | 4   * +-+-+-+ .
    |        +-+-+-+             +-+-+-+-+       | | | | .
    |        | | | | 5           | | | | | 5     +-+-+-+
    |        +-+-+-+             +-+-+-+-+       | | | | iel-1
    |                                    +-+-+-+
    |
    |        xder2               deriv2          xyze^T
    |
    |
    |                                     +-                  -+
    |                                 | d^2x   d^2y   d^2z |
    |                                 | ----   ----   ---- |
    |                         | dr^2   dr^2   dr^2 |
    |                         |                    |
    |                         | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                         | ds^2   ds^2   ds^2 |
    |                         |                    |
    |                         | d^2x   d^2y   d^2z |
    |                         | ----   ----   ---- |
    |                         | dt^2   dt^2   dt^2 |
    |               yields    xder2  =    |                    |
    |                                     | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                                     | drds   drds   drds |
    |                                     |                    |
    |                                     | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                                     | drdt   drdt   drdt |
    |                                     |                    |
    |                                     | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                                     | dsdt   dsdt   dsdt |
    |                         +-                  -+
    |
    |
    */

    const BlitzMat xder2_(blitz::sum(deriv2_(i,k)*xyze_(j,k),k));

    /*
    |        0...iel-1             0 1 2
    |        +-+-+-+-+            +-+-+-+
    |        | | | | | 0          | | | | 0
    |        +-+-+-+-+            +-+-+-+            0...iel-1
    |        | | | | | 1          | | | | 1         +-+-+-+-+
    |        +-+-+-+-+            +-+-+-+           | | | | | 0
    |        | | | | | 2          | | | | 2         +-+-+-+-+
    |        +-+-+-+-+       =    +-+-+-+       *   | | | | | 1 * (-1)
    |        | | | | | 3          | | | | 3         +-+-+-+-+
    |        +-+-+-+-+            +-+-+-+           | | | | | 2
    |        | | | | | 4          | | | | 4         +-+-+-+-+
    |        +-+-+-+-+            +-+-+-+
    |        | | | | | 5          | | | | 5          derxy
    |        +-+-+-+-+            +-+-+-+
    |
    |       chainrulerhs          xder2
    */

    BlitzMat derxy2_(deriv2_-blitz::sum(xder2_(i,k)*derxy_(k,j),k));

    /*
    |        0...iel-1            0...iel-1         0...iel-1
    |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
    |        | | | | | 0          | | | | | 0       | | | | | 0
    |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
    |        | | | | | 1          | | | | | 1       | | | | | 1
    |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
    |        | | | | | 2          | | | | | 2       | | | | | 2
    |        +-+-+-+-+       =    +-+-+-+-+    +    +-+-+-+-+
    |        | | | | | 3          | | | | | 3       | | | | | 3
    |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
    |        | | | | | 4          | | | | | 4       | | | | | 4
    |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
    |        | | | | | 5          | | | | | 5       | | | | | 5
    |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
    |
    |       chainrulerhs         chainrulerhs        deriv2
    */

    //derxy2_ += deriv2_;

    /* make LR decomposition and solve system for all right hand sides
     * (i.e. the components of chainrulerhs)
    |
    |          0  1  2  3  4  5         i        i
    |        +--+--+--+--+--+--+       +-+      +-+
    |        |  |  |  |  |  |  | 0     | | 0    | | 0
    |        +--+--+--+--+--+--+       +-+      +-+
    |        |  |  |  |  |  |  | 1     | | 1    | | 1
    |        +--+--+--+--+--+--+       +-+      +-+
    |        |  |  |  |  |  |  | 2     | | 2    | | 2
    |        +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
    |        |  |  |  |  |  |  | 3     | | 3    | | 3
    |        +--+--+--+--+--+--+       +-+      +-+
    |        |  |  |  |  |  |  | 4     | | 4    | | 4
    |        +--+--+--+--+--+--+       +-+      +-+
    |        |  |  |  |  |  |  | 5     | | 5    | | 5
    |        +--+--+--+--+--+--+       +-+      +-+
    |                                   |        |
    |                                   |        |
    |                                   derxy2[i]|
    |                                    |
    |                                    chainrulerhs[i]
    |
    |   yields
    |
    |                      0...iel-1
    |                      +-+-+-+-+
    |                      | | | | | 0 = drdr
    |                      +-+-+-+-+
    |                      | | | | | 1 = dsds
    |                      +-+-+-+-+
    |                      | | | | | 2 = dtdt
    |            derxy2 =  +-+-+-+-+
    |                      | | | | | 3 = drds
    |                      +-+-+-+-+
    |                      | | | | | 4 = drdt
    |                      +-+-+-+-+
    |                      | | | | | 5 = dsdt
    |                  +-+-+-+-+
    */
    const int numnode = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
    Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,numnode);

    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(bm);

    // No need for a separate rhs. We assemble the rhs to the solution
    // vector. The solver will destroy the rhs and return the solution.
    solver.SetVectors(ederxy2,ederxy2);
    solver.Solve();

    return derxy2_;
}
  

//! calculate element metric tensor in gausspoint (r,s)
BlitzMat computeMetricTensor(
        const BlitzMat xs,                    ///< element nodal coordinates in x space
        const BlitzMat deriv                  ///< shape functions derivative with respect to r and s
        )
{
    // init blitz indices
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index
      
    // compute dXYZ / drs
    const BlitzMat dxyzdrs(blitz::sum(deriv(i,k)*xs(j,k),k));

    /* compute covariant metric tensor G for surface element
    **                        | g11   g12 |
    **                    G = |           |
    **                        | g12   g22 |
    ** where (o denotes the inner product, xyz a vector)
    **
    **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
    ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
    **        dr     dr            dr     ds            ds     ds
    */
    const BlitzMat metrictensor(blitz::sum(dxyzdrs(i,k)*dxyzdrs(j,k),k));

    return metrictensor;
}



//! calculate element normal in gausspoint (r,s)
BlitzVec computeNormalToBoundaryElement(
        const BlitzMat deriv,                  ///< shape functions derivative
        const BlitzMat xs,                     ///< element nodal coordinates in x space
        const double r,                        ///< coord in r-direction
        const double s                         ///< coord in s-direction
        )
{
    // init blitz indices
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index
      
    // compute dXYZ / drs
    const BlitzMat dxyzdrs(blitz::sum(deriv(i,k)*xs(j,k),k));

    /* compute covariant metric tensor G for surface element
    **                        | g11   g12 |
    **                    G = |           |
    **                        | g12   g22 |
    ** where (o denotes the inner product, xyz a vector)
    **
    **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
    ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
    **        dr     dr            dr     ds            ds     ds
    */
    const BlitzMat metrictensor(blitz::sum(dxyzdrs(i,k)*dxyzdrs(j,k),k));
    
    BlitzVec unrm(3);
    unrm(0) = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    unrm(1) = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    unrm(2) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);
    
    // scale
    const double length = sqrt(unrm(0)*unrm(0) + unrm(1)*unrm(1) + unrm(2)*unrm(2));
    
    BlitzVec unit_normal(3);
    unit_normal = 1.0/length * unrm;

    return unit_normal;
}
  
  
  
    
template<XFEM::AssemblyType>
int getNumParam(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        int numnode)
{
    return dofman.NumDofPerField(field);
}
template<> 
int getNumParam<XFEM::standard_assembly>(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        int numnode)
{
    return numnode;
}
  
  /*!
    Calculate inverse of a matrix

   */
BlitzMat Inverse3x3(const BlitzMat& xjm_, const double det)
{
    BlitzMat xji_(3,3,blitz::ColumnMajorArray<2>());
    xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
    xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
    xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
    xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
    xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
    xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
    xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
    xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
    xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;
    return xji_;
}
  
/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void Sysmat(
        const DRT::ELEMENTS::XFluid3*     ele,           ///< the element those matrix is calculated
        const RCP<XFEM::InterfaceHandle>  ih,            ///< connection to the interface handler
        const XFEM::ElementDofManager&    dofman,        ///< dofmanager of the current element
        const std::vector<double>&        locval,        ///< nodal unknowns at n+1, i
        BlitzMat&                         estif,         ///< element matrix to calculate
        BlitzVec&                         eforce,        ///< element rhs to calculate
        struct _MATERIAL*                 material,      ///< fluid material
        double                            pseudotime,    ///< current simulation pseudotime
        bool                              newton,        ///< full Newton or fixed-point-like
        bool                              pstab,         ///< flag for stabilisation
        bool                              supg,          ///< flag for stabilisation
        bool                              vstab,         ///< flag for stabilisation
        bool                              cstab          ///< flag for stabilisation
        )
{
    // number of nodes for element
    const int numnode = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
      
    // get node coordinates
    const BlitzMat xyze_(DRT::UTILS::PositionArray(ele));

    // dead load in element nodes
    //////////////////////////////////////////////////// , BlitzMat edeadng_(BodyForce(ele->Nodes(),pseudotime));

    // get viscosity
    // check here, if we really have a fluid !!
    dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
    const double visc = material->m.fluid->viscosity;
    const double reciproke_visc = 1.0 / visc;

    // flag for higher order elements
    const bool higher_order_ele = ele->isHigherOrderElement<DISTYPE>();

    DRT::UTILS::GaussRule3D gaussrule;
    if (ih->ElementIntersected(ele->Id()))
    {
        // attention: don't use gausspoint rules that have gausspoints on the element boundary 
        // (so called closed rules)
        // -> enrichment values might be undefined 
        gaussrule = DRT::UTILS::intrule_tet_11point;
    }
    else
    {
        gaussrule = ele->GetGaussrule();
    }
    
    LocalAssembler<DISTYPE, ASSTYPE> assembler(dofman, estif, eforce);
    
    const blitz::Range _  = blitz::Range::all();
    
    // We define the variables i,j,k to be indices to blitz arrays.
    // These are used for array expressions, that is matrix-vector
    // products in the following.
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index

    const DRT::Element::DiscretizationType stressdistype = XFLUID::getStressInterpolationType3D(DISTYPE);
    
    // figure out whether we have stress unknowns at all
    const bool stress_unknowns_present = (getNumParam<ASSTYPE>(dofman, Tauxx, 0) > 0);
//    const bool velocity_unknowns_present = (getNumParam<ASSTYPE>(dofman, Velx, 1) > 0);
//    const bool pressure_unknowns_present = (getNumParam<ASSTYPE>(dofman, Pres, 1) > 0);
//    cout << endl;
//    if (ASSTYPE == XFEM::standard_assembly)
//        cout << "standard assembly" << endl;
//    else
//        cout << "xfem assembly" << endl;
//            
//    cout << "stress unknowns present  : " << stress_unknowns_present << endl;
//    cout << "velocity unknowns present: " << velocity_unknowns_present << endl;
//    cout << "pressure unknowns present: " << pressure_unknowns_present << endl;
    
    
    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = getNumParam<ASSTYPE>(dofman, Velx, numnode);
    const int numparamvely = getNumParam<ASSTYPE>(dofman, Vely, numnode);
    const int numparamvelz = getNumParam<ASSTYPE>(dofman, Velz, numnode);
    const int numparampres = getNumParam<ASSTYPE>(dofman, Pres, numnode);
    const int numparamtauxx = getNumParam<ASSTYPE>(dofman, Tauxx, 1);
    const int numparamtauyy = getNumParam<ASSTYPE>(dofman, Tauyy, 1);
    const int numparamtauzz = getNumParam<ASSTYPE>(dofman, Tauzz, 1);
    const int numparamtauxy = getNumParam<ASSTYPE>(dofman, Tauxy, 1);
    const int numparamtauxz = getNumParam<ASSTYPE>(dofman, Tauxz, 1);
    const int numparamtauyz = getNumParam<ASSTYPE>(dofman, Tauyz, 1);
    

    
    const vector<int> velxdof = dofman.LocalDofPosPerField(Velx);
    const vector<int> velydof = dofman.LocalDofPosPerField(Vely);
    const vector<int> velzdof = dofman.LocalDofPosPerField(Velz);
    const vector<int> presdof = dofman.LocalDofPosPerField(Pres);
    const vector<int> tauxxdof = dofman.LocalDofPosPerField(Tauxx);
    const vector<int> tauyydof = dofman.LocalDofPosPerField(Tauyy);
    const vector<int> tauzzdof = dofman.LocalDofPosPerField(Tauzz);
    const vector<int> tauxydof = dofman.LocalDofPosPerField(Tauxy);
    const vector<int> tauxzdof = dofman.LocalDofPosPerField(Tauxz);
    const vector<int> tauyzdof = dofman.LocalDofPosPerField(Tauyz);

//    cout << "numparamvelx  " << numparamvelx << "  " <<  velxdof.size() << endl;
//    cout << "numparamtauxx " << numparamtauxx << "  " <<  tauxxdof.size() << endl;
//    cout << "numparamtauyy " << numparamtauyy << "  " <<  tauyydof.size() << endl;
//    cout << "numparamtauzz " << numparamtauzz << "  " <<  tauzzdof.size() << endl;
//    cout << "numparamtauxy " << numparamtauxy << "  " <<  tauxydof.size() << endl;
//    cout << "numparamtauxz " << numparamtauxz << "  " <<  tauxzdof.size() << endl;
//    cout << "numparamtauyz " << numparamtauyz << "  " <<  tauyzdof.size() << endl;

//    if (numparamtauxx > 0)
//        cout << tauxxdof[0] << " " << tauxxdof[1] << " " << tauxxdof[2] << " " << tauxxdof[3] << " " << endl;

    //    if (numparamvelx != numparamvely or numparamvelx != numparamvelz or numparamvelx != numparampres)
    //        dserror("for now, we enrich velocity and pressure together");
    
    // split velocity and pressure (and stress)
    blitz::Array<double, 1> eprenp(numparampres);
    blitz::Array<double, 2> evelnp(3,numparamvelx,blitz::ColumnMajorArray<2>());
    blitz::Array<double, 2> etau(6,numparamtauxx,blitz::ColumnMajorArray<2>());
    
    for (int iparam=0; iparam<numparamvelx; ++iparam)   evelnp(0,iparam) = locval[velxdof[iparam]];
    for (int iparam=0; iparam<numparamvely; ++iparam)   evelnp(1,iparam) = locval[velydof[iparam]];
    for (int iparam=0; iparam<numparamvelz; ++iparam)   evelnp(2,iparam) = locval[velzdof[iparam]];
    for (int iparam=0; iparam<numparampres; ++iparam)   eprenp(iparam) = locval[presdof[iparam]];
    if (stress_unknowns_present)
    {
        for (int iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = locval[tauxxdof[iparam]];
        for (int iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = locval[tauyydof[iparam]];
        for (int iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = locval[tauzzdof[iparam]];
        for (int iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = locval[tauxydof[iparam]];
        for (int iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = locval[tauxzdof[iparam]];
        for (int iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = locval[tauyzdof[iparam]];
    }
    
    // stabilization parameter
    const double hk = XFLUID::HK<DISTYPE>(evelnp,visc,xyze_);
    const double mk = XFLUID::MK<DISTYPE>();
    
    // information about domain integration cells
    const XFEM::DomainIntCells   domainIntCells   = ih->GetDomainIntCells(ele->Id(),DISTYPE);
    // loop over integration cells
    for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {

        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            XFEM::PosEtaDomain pos_eta(3);
            pos_eta(0) = intpoints.qxg[iquad][0];
            pos_eta(1) = intpoints.qxg[iquad][1];
            pos_eta(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            const XFEM::PosXiDomain posXiDomain(XFEM::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta));
            const double detcell = XFEM::detEtaToXi3D<ASSTYPE>(*cell, pos_eta);
            
            // shape functions and their first derivatives
            const BlitzVec funct_(DRT::UTILS::shape_function_3D(posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE));
            const BlitzMat deriv_(DRT::UTILS::shape_function_3D_deriv1(posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE));
      
            BlitzVec funct_stress(DRT::UTILS::getNumberOfElementNodes(stressdistype));
            BlitzMat deriv_stress(3, DRT::UTILS::getNumberOfElementNodes(stressdistype));
            if (stress_unknowns_present)
            {
                DRT::UTILS::shape_function_3D(funct_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
                DRT::UTILS::shape_function_3D_deriv1(deriv_stress,posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype);
                //cout << "funct_stress" << funct_stress << endl;
                //cout << "deriv_stress" << deriv_stress << endl;
            }

            // position of the gausspoint in physical coordinates
            const BlitzVec gauss_pos_xyz(blitz::sum(funct_(j)*xyze_(i,j),j));
      
            // get transposed of the jacobian matrix d x / d \xi
            const BlitzMat xjm_(blitz::sum(deriv_(i,k)*xyze_(j,k),k));
            const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
                               xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
                               xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
                               xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
                               xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
                               xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
            const double fac = intpoints.qwgt[iquad]*det*detcell;

            if (det < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %lf", ele->Id(), det);
            }

            // inverse of jacobian
            const BlitzMat xji_(Inverse3x3(xjm_, det));

            // compute global derivates
            const BlitzMat derxy_(blitz::sum(xji_(i,k)*deriv_(k,j),k));
            const BlitzMat derxy_stress(blitz::sum(xji_(i,k)*deriv_stress(k,j),k));

            // compute second global derivative
            BlitzMat derxy2_(6,numnode,blitz::ColumnMajorArray<2>());
            if (higher_order_ele)
            {
                BlitzMat deriv2_(6,numnode,blitz::ColumnMajorArray<2>());
                DRT::UTILS::shape_function_3D_deriv2(deriv2_,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
                derxy2_ = gder2<DISTYPE>(xjm_, derxy_, deriv2_, xyze_);
            }
            else
            {
                derxy2_ = 0.;
            }

            // after this call, one should only use the enriched shape functions and derivatives!
            BlitzVec enr_funct_(numparamvelx);
            BlitzMat enr_derxy_(3,numparamvelx,blitz::ColumnMajorArray<2>());
            BlitzMat enr_derxy2_(6,numparamvelx,blitz::ColumnMajorArray<2>());
      
            BlitzVec enr_funct_stress(numparamtauxx);
            BlitzMat enr_derxy_stress(3,numparamtauxx,blitz::ColumnMajorArray<2>());
            if (ASSTYPE == XFEM::xfem_assembly)
            {
                // shape function for nodal dofs
                XFEM::ComputeEnrichedShapefunction(
                        *ele,
                        ih,
                        dofman,
                        Velx,
                        gauss_pos_xyz,
                        XFEM::Enrichment::approachUnknown,
                        funct_,
                        derxy_,
                        derxy2_, 
                        enr_funct_,
                        enr_derxy_,
                        enr_derxy2_);
          
                if (stress_unknowns_present)
                {
                    // shape functions for element dofs
                    XFEM::ComputeEnrichedStressShapefunction(
                            *ele,
                            ih,
                            dofman,
                            Tauxx,
                            gauss_pos_xyz,
                            XFEM::Enrichment::approachUnknown,
                            funct_stress,
                            derxy_stress,
                            enr_funct_stress,
                            enr_derxy_stress);
                }
            }
            else
            {
                enr_funct_  = funct_;
                enr_derxy_  = derxy_;
                enr_derxy2_ = derxy2_;
          
                enr_funct_stress  = funct_stress;
                enr_derxy_stress  = derxy_stress;
            }
      
            // create vievs on the shape function arrays for easy handling in the assembly process
            const BlitzVec shp(enr_funct_);
            const BlitzVec shp_dx(enr_derxy_(0,_));
            const BlitzVec shp_dy(enr_derxy_(1,_));
            const BlitzVec shp_dz(enr_derxy_(2,_));
            const BlitzVec shp_dxdx(enr_derxy2_(0,_)); const BlitzVec shp_dxdy(enr_derxy2_(3,_)); const BlitzVec shp_dxdz(enr_derxy2_(4,_));
            const BlitzVec shp_dydx(shp_dxdy);         const BlitzVec shp_dydy(enr_derxy2_(1,_)); const BlitzVec shp_dydz(enr_derxy2_(5,_));
            const BlitzVec shp_dzdx(shp_dxdz);         const BlitzVec shp_dzdy(shp_dydz);         const BlitzVec shp_dzdz(enr_derxy2_(2,_));
            
            const BlitzVec shp_tau(enr_funct_stress);
            const BlitzVec shp_tau_dx(enr_derxy_stress(0,_));
            const BlitzVec shp_tau_dy(enr_derxy_stress(1,_));
            const BlitzVec shp_tau_dz(enr_derxy_stress(2,_));
      
            // get velocities (n+g,i) at integration point
            const BlitzVec velint_(blitz::sum(enr_funct_(j)*evelnp(i,j),j));

            // get velocity (np,i) derivatives at integration point
            const BlitzMat vderxy_(blitz::sum(enr_derxy_(j,k)*evelnp(i,k),k));
      
            // calculate 2nd velocity derivatives at integration point
            BlitzMat vderxy2_(3,6,blitz::ColumnMajorArray<2>());
            if (higher_order_ele)
            {
                vderxy2_ = blitz::sum(enr_derxy2_(j,k)*evelnp(i,k),k);
            }
            else
            {
                vderxy2_ = 0.;
            }

            // get pressure gradients
            const BlitzVec gradp_(blitz::sum(enr_derxy_(i,j)*eprenp(j),j));
    
            // get pressure
            const double press = blitz::sum(enr_funct_*eprenp);
    
            // get viscous stress unknowns
            BlitzMat tau(3,3,blitz::ColumnMajorArray<2>());
            if (stress_unknowns_present)
            {
                tau(0,0) = blitz::sum(shp_tau(_)*etau(0,_)); tau(0,1) = blitz::sum(shp_tau(_)*etau(3,_)); tau(0,2) = blitz::sum(shp_tau(_)*etau(4,_));
                tau(1,0) = tau(0,1);                         tau(1,1) = blitz::sum(shp_tau(_)*etau(1,_)); tau(1,2) = blitz::sum(shp_tau(_)*etau(5,_));
                tau(2,0) = tau(0,2);                         tau(2,1) = tau(1,2);                         tau(2,2) = blitz::sum(shp_tau(_)*etau(2,_));
            }
          
//            BlitzVec nabla_dot_tau(3);
//            if (stress_unknowns_present)
//            {
//                nabla_dot_tau(0) = blitz::sum(shp_tau_dx(_)*etau(0,_)) + blitz::sum(shp_tau_dy(_)*etau(3,_)) + blitz::sum(shp_tau_dz(_)*etau(4,_));
//                nabla_dot_tau(1) = blitz::sum(shp_tau_dx(_)*etau(3,_)) + blitz::sum(shp_tau_dy(_)*etau(1,_)) + blitz::sum(shp_tau_dz(_)*etau(5,_));
//                nabla_dot_tau(2) = blitz::sum(shp_tau_dx(_)*etau(4,_)) + blitz::sum(shp_tau_dy(_)*etau(5,_)) + blitz::sum(shp_tau_dz(_)*etau(2,_));
//            }
            
            //const double trace_tau = (tau(0, 0) + tau(1, 1) + tau(2, 2));
      
      
//            const BlitzMat eps(0.5*(vderxy_(i,j) + vderxy_(j,i));
//            const BlitzMat epstau(tau*reciproke_visc);
//      
//            if (ele->Id() == 0)
//            {
//                cout << endl;
//                cout << "eps^tau: " << epstau << endl;
//                cout << "eps^u:   " << eps << endl;
//                cout << "pressure = " << press << endl;
//            }

            
            // get bodyforce in gausspoint
            BlitzVec bodyforce_(3);
            bodyforce_ = 0.0;
            //////////////////////////////////////////BlitzVec bodyforce_(blitz::sum(enr_edeadng_(i,j)*enr_funct_(j),j));

            

            // get velocity norm
            const double vel_norm = sqrt(blitz::sum(velint_*velint_));

            // normed velocity at element centre
            BlitzVec velino_(3);
            if (vel_norm>=1e-6)
            {
                velino_ = velint_/vel_norm;
            }
            else
            {
                velino_ = 0.;
                velino_(0) = 1;
            }

            // get streamlength
            const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
            const double strle = 2.0/val;

            // calculate tau: stabilization parameters for stationary case

            // compute tau_Mu
            const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
            const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
            const double tau_stab_0 = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);

            // compute tau_Mp
            const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
            const double xi_tau_mp = DMAX(re_tau_mp,1.0);
            const double tau_stab_1 = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);

            // compute tau_C
            const double xi_tau_c = DMIN(re_tau_mp, 1.0);
            const double tau_stab_2 = 0.5*vel_norm*hk*xi_tau_c;
      
            // stabilisation parameter
            const double tau_M  = tau_stab_0*fac;
            const double tau_Mp = tau_stab_1*fac;
            const double tau_C  = tau_stab_2*fac;

            /*------------------------- evaluate rhs vector at integration point ---*/
            //   rhsint_ = histvec_(i) + bodyforce_(i);
            // histvec is always zero in stationary case (!):
            const BlitzVec rhsint_(bodyforce_(i));

            /*----------------- get numerical representation of single operators ---*/

            /* Convective term  u_old * grad u_old: */
            const BlitzVec conv_old_(blitz::sum(vderxy_(i, j)*velint_(j), j));
            
            /* Viscous term  div epsilon(u_old) */
            BlitzVec visc_old_(3);
            visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,3) + vderxy2_(0,2) + vderxy2_(2,4));
            visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,3) + vderxy2_(1,2) + vderxy2_(2,5));
            visc_old_(2) = vderxy2_(2,2) + 0.5 * (vderxy2_(2,0) + vderxy2_(0,4) + vderxy2_(2,1) + vderxy2_(1,5));
            
            // evaluate residual once for all stabilisation right hand sides
            BlitzVec res_old(-rhsint_+(conv_old_+gradp_-2.0*visc*visc_old_));
      
            /* Reactive term  u:  funct */
            /* linearise convective term */

            /*--- convective part u_old * grad (funct) --------------------------*/
            /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
             with  N .. form function matrix                                   */
            const BlitzVec enr_conv_c_(blitz::sum(enr_derxy_(j,i)*velint_(j), j));

//              /*--- convective grid part u_G * grad (funct) -----------------------*/
//              /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
//              enr_conv_g_ = 0.0;


          /*--- viscous term  - grad * epsilon(u): ----------------------------*/
          /*   /                                                \
               |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
             1 |                                                |
             - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
             2 |                                                |
               |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
               \                                                /
    
               with N_x .. x-line of N
               N_y .. y-line of N                                             */
            blitz::Array<double,3> enr_viscs2_(3,3,numparamvelx,blitz::ColumnMajorArray<3>());
            enr_viscs2_(0,0,_) = 0.5 * (2.0 * shp_dxdx + shp_dydy + shp_dzdz);
            enr_viscs2_(0,1,_) = 0.5 *  shp_dxdy;
            enr_viscs2_(0,2,_) = 0.5 *  shp_dxdz;
            enr_viscs2_(1,0,_) = 0.5 *  shp_dydx;
            enr_viscs2_(1,1,_) = 0.5 * (shp_dxdx + 2.0 * shp_dydy + shp_dzdz);
            enr_viscs2_(1,2,_) = 0.5 *  shp_dydz;
            enr_viscs2_(2,0,_) = 0.5 *  shp_dzdx;
            enr_viscs2_(2,1,_) = 0.5 *  shp_dzdy;
            enr_viscs2_(2,2,_) = 0.5 * (shp_dxdx + shp_dydy + 2.0 * shp_dzdz);


            /*--------------------------------- now build single stiffness terms ---*/

            //----------------------------------------------------------------------
            //                            GALERKIN PART

            /* convection, convective part */
            /*
                         /                       \
                        |      / n+1       \      |
                        | v , | u   o nabla | Du  |
                        |      \ (i)       /      |
                         \                       /
            */
            assembler.template Matrix<Velx,Velx>(shp, fac, enr_conv_c_);
            assembler.template Matrix<Vely,Vely>(shp, fac, enr_conv_c_);
            assembler.template Matrix<Velz,Velz>(shp, fac, enr_conv_c_);
            
            assembler.template Vector<Velx>(shp, -fac*velint_(0)*vderxy_(0,0));
            assembler.template Vector<Velx>(shp, -fac*velint_(1)*vderxy_(0,1));
            assembler.template Vector<Velx>(shp, -fac*velint_(2)*vderxy_(0,2));
            
            assembler.template Vector<Vely>(shp, -fac*velint_(0)*vderxy_(1,0));
            assembler.template Vector<Vely>(shp, -fac*velint_(1)*vderxy_(1,1));
            assembler.template Vector<Vely>(shp, -fac*velint_(2)*vderxy_(1,2));
            
            assembler.template Vector<Velz>(shp, -fac*velint_(0)*vderxy_(2,0));
            assembler.template Vector<Velz>(shp, -fac*velint_(1)*vderxy_(2,1));
            assembler.template Vector<Velz>(shp, -fac*velint_(2)*vderxy_(2,2));
            
            if (newton)
            {
                /*  convection, reactive part */
                /*
                       /                         \
                      |      /          \   n+1   |
                      | v , | Du o nabla | u      |
                      |      \          /   (i)   |
                       \                         /
                */
                assembler.template Matrix<Velx,Velx>(shp, fac*vderxy_(0,0), shp);
                assembler.template Matrix<Velx,Vely>(shp, fac*vderxy_(0,1), shp);
                assembler.template Matrix<Velx,Velz>(shp, fac*vderxy_(0,2), shp);
                assembler.template Matrix<Vely,Velx>(shp, fac*vderxy_(1,0), shp);
                assembler.template Matrix<Vely,Vely>(shp, fac*vderxy_(1,1), shp);
                assembler.template Matrix<Vely,Velz>(shp, fac*vderxy_(1,2), shp);
                assembler.template Matrix<Velz,Velx>(shp, fac*vderxy_(2,0), shp);
                assembler.template Matrix<Velz,Vely>(shp, fac*vderxy_(2,1), shp);
                assembler.template Matrix<Velz,Velz>(shp, fac*vderxy_(2,2), shp);
            }
            
            /* Viskositaetsterm */
            /*
                          /                        \
                         |       / \         /  \   |
                         |  eps | v | , tau | Du |  |
                         |       \ /         \  /   |
                          \                        /
            */
            assembler.template Matrix<Velx,Velx>(shp_dx, 2.0*visc*fac, shp_dx);
            assembler.template Matrix<Velx,Velx>(shp_dy,     visc*fac, shp_dy);
            assembler.template Matrix<Velx,Velx>(shp_dz,     visc*fac, shp_dz);
            assembler.template Matrix<Velx,Vely>(shp_dy,     visc*fac, shp_dx);
            assembler.template Matrix<Velx,Velz>(shp_dz,     visc*fac, shp_dx);
            
            assembler.template Matrix<Vely,Velx>(shp_dx,     visc*fac, shp_dy);
            assembler.template Matrix<Vely,Vely>(shp_dx,     visc*fac, shp_dx);
            assembler.template Matrix<Vely,Vely>(shp_dy, 2.0*visc*fac, shp_dy);
            assembler.template Matrix<Vely,Vely>(shp_dz,     visc*fac, shp_dz);
            assembler.template Matrix<Vely,Velz>(shp_dz,     visc*fac, shp_dy);
            
            assembler.template Matrix<Velz,Velx>(shp_dx,     visc*fac, shp_dz);
            assembler.template Matrix<Velz,Vely>(shp_dy,     visc*fac, shp_dz);
            assembler.template Matrix<Velz,Velz>(shp_dx,     visc*fac, shp_dx);
            assembler.template Matrix<Velz,Velz>(shp_dy,     visc*fac, shp_dy);
            assembler.template Matrix<Velz,Velz>(shp_dz, 2.0*visc*fac, shp_dz);
            
            assembler.template Vector<Velx>(shp_dx, -2.0*visc*fac*vderxy_(0, 0));
            assembler.template Vector<Velx>(shp_dy,     -visc*fac*vderxy_(0, 1));
            assembler.template Vector<Velx>(shp_dy,     -visc*fac*vderxy_(1, 0));
            assembler.template Vector<Velx>(shp_dz,     -visc*fac*vderxy_(0, 2));
            assembler.template Vector<Velx>(shp_dz,     -visc*fac*vderxy_(2, 0));
            
            assembler.template Vector<Vely>(shp_dx,     -visc*fac*vderxy_(1, 0));
            assembler.template Vector<Vely>(shp_dx,     -visc*fac*vderxy_(0, 1));
            assembler.template Vector<Vely>(shp_dy, -2.0*visc*fac*vderxy_(1, 1));
            assembler.template Vector<Vely>(shp_dz,     -visc*fac*vderxy_(1, 2));
            assembler.template Vector<Vely>(shp_dz,     -visc*fac*vderxy_(2, 1));
            
            assembler.template Vector<Velz>(shp_dx,     -visc*fac*vderxy_(2, 0));
            assembler.template Vector<Velz>(shp_dx,     -visc*fac*vderxy_(0, 2));
            assembler.template Vector<Velz>(shp_dy,     -visc*fac*vderxy_(2, 1));
            assembler.template Vector<Velz>(shp_dy,     -visc*fac*vderxy_(1, 2));
            assembler.template Vector<Velz>(shp_dz, -2.0*visc*fac*vderxy_(2, 2));
            
            /* Druckterm */
            /*
                            /                \
                           |                  |
                           |  nabla o v , Dp  |
                           |                  |
                            \                /
            */
            assembler.template Matrix<Velx,Pres>(shp_dx, -fac, shp);
            assembler.template Matrix<Vely,Pres>(shp_dy, -fac, shp);
            assembler.template Matrix<Velz,Pres>(shp_dz, -fac, shp);
            
            assembler.template Vector<Velx>(shp_dx, press*fac);
            assembler.template Vector<Vely>(shp_dy, press*fac);
            assembler.template Vector<Velz>(shp_dz, press*fac);
            
            /* Divergenzfreiheit - continuity equation*/
            /*
                           /              \
                          |                |
                          | q , nabla o Du |
                          |                |
                           \              /
            */
            assembler.template Matrix<Pres,Velx>(shp, fac, shp_dx);
            assembler.template Matrix<Pres,Vely>(shp, fac, shp_dy);
            assembler.template Matrix<Pres,Velz>(shp, fac, shp_dz);
            
            assembler.template Vector<Pres>(shp, -fac*vderxy_(0,0));
            assembler.template Vector<Pres>(shp, -fac*vderxy_(1,1));
            assembler.template Vector<Pres>(shp, -fac*vderxy_(2,2));
            
            // source term of the right hand side
            assembler.template Vector<Velx>(shp, fac*rhsint_(0));
            assembler.template Vector<Vely>(shp, fac*rhsint_(1));
            assembler.template Vector<Velz>(shp, fac*rhsint_(2));
            
            
            // Hellinger-Reissner terms
            if (stress_unknowns_present)
            {
                /* quadratic part */
                /*
                               /                    \
                              |                      |
                              | virt tau , eps(Dtau) |
                              |                      |
                               \                    /
                */
                assembler.template Matrix<Tauxx,Tauxx>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauxy,Tauxy>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauxz,Tauxz>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauyx,Tauyx>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauyy,Tauyy>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauyz,Tauyz>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauzx,Tauzx>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauzy,Tauzy>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                assembler.template Matrix<Tauzz,Tauzz>(shp_tau, -0.5*reciproke_visc*fac, shp_tau);
                
                assembler.template Vector<Tauxx>(shp_tau, 0.5*reciproke_visc*fac*tau(0,0));
                assembler.template Vector<Tauxy>(shp_tau, 0.5*reciproke_visc*fac*tau(0,1));
                assembler.template Vector<Tauxz>(shp_tau, 0.5*reciproke_visc*fac*tau(0,2));
                assembler.template Vector<Tauyx>(shp_tau, 0.5*reciproke_visc*fac*tau(1,0));
                assembler.template Vector<Tauyy>(shp_tau, 0.5*reciproke_visc*fac*tau(1,1));
                assembler.template Vector<Tauyz>(shp_tau, 0.5*reciproke_visc*fac*tau(1,2));
                assembler.template Vector<Tauzx>(shp_tau, 0.5*reciproke_visc*fac*tau(2,0));
                assembler.template Vector<Tauzy>(shp_tau, 0.5*reciproke_visc*fac*tau(2,1));
                assembler.template Vector<Tauzz>(shp_tau, 0.5*reciproke_visc*fac*tau(2,2));
                
                /* rectangular part */
                /*
                               /                  \
                              |                    |
                              | virt tau , eps(Du) |
                              |                    |
                               \                  /
                */
                assembler.template Matrix<Tauxx,Velx>(shp_tau, 2.0*fac*0.5, shp_dx);
                assembler.template Matrix<Tauxy,Velx>(shp_tau,     fac*0.5, shp_dy);
                assembler.template Matrix<Tauxy,Vely>(shp_tau,     fac*0.5, shp_dx);
                assembler.template Matrix<Tauxz,Velx>(shp_tau,     fac*0.5, shp_dz);
                assembler.template Matrix<Tauxz,Velz>(shp_tau,     fac*0.5, shp_dx);
                
                assembler.template Matrix<Tauyx,Velx>(shp_tau,     fac*0.5, shp_dy);
                assembler.template Matrix<Tauyx,Vely>(shp_tau,     fac*0.5, shp_dx);
                assembler.template Matrix<Tauyy,Vely>(shp_tau, 2.0*fac*0.5, shp_dy);
                assembler.template Matrix<Tauyz,Velz>(shp_tau,     fac*0.5, shp_dz);
                assembler.template Matrix<Tauyz,Vely>(shp_tau,     fac*0.5, shp_dx);
                
                assembler.template Matrix<Tauzx,Velx>(shp_tau,     fac*0.5, shp_dz);
                assembler.template Matrix<Tauzx,Velz>(shp_tau,     fac*0.5, shp_dx);
                assembler.template Matrix<Tauzy,Vely>(shp_tau,     fac*0.5, shp_dz);
                assembler.template Matrix<Tauzy,Velz>(shp_tau,     fac*0.5, shp_dy);
                assembler.template Matrix<Tauzz,Velz>(shp_tau, 2.0*fac*0.5, shp_dz);
                
                assembler.template Vector<Tauxx>(shp_tau, -2.0*fac*0.5*vderxy_(0, 0));
                assembler.template Vector<Tauxy>(shp_tau, -    fac*0.5*vderxy_(0, 1));
                assembler.template Vector<Tauxy>(shp_tau, -    fac*0.5*vderxy_(1, 0));
                assembler.template Vector<Tauxz>(shp_tau, -    fac*0.5*vderxy_(0, 2));
                assembler.template Vector<Tauxz>(shp_tau, -    fac*0.5*vderxy_(2, 0));
                
                assembler.template Vector<Tauyx>(shp_tau, -    fac*0.5*vderxy_(1, 0));
                assembler.template Vector<Tauyx>(shp_tau, -    fac*0.5*vderxy_(0, 1));
                assembler.template Vector<Tauyy>(shp_tau, -2.0*fac*0.5*vderxy_(1, 1));
                assembler.template Vector<Tauyz>(shp_tau, -    fac*0.5*vderxy_(1, 2));
                assembler.template Vector<Tauyz>(shp_tau, -    fac*0.5*vderxy_(2, 1));
                
                assembler.template Vector<Tauzx>(shp_tau, -    fac*0.5*vderxy_(2, 0));
                assembler.template Vector<Tauzx>(shp_tau, -    fac*0.5*vderxy_(0, 2));
                assembler.template Vector<Tauzy>(shp_tau, -    fac*0.5*vderxy_(2, 1));
                assembler.template Vector<Tauzy>(shp_tau, -    fac*0.5*vderxy_(1, 2));
                assembler.template Vector<Tauzz>(shp_tau, -2.0*fac*0.5*vderxy_(2, 2));
                
                /*
                               /              \
                              |                |
                              | virt tau , p I |
                              |                |
                               \              /
                */
                assembler.template Matrix<Tauxx,Pres>(shp_tau, -0.5*reciproke_visc*fac, shp);
                assembler.template Matrix<Tauyy,Pres>(shp_tau, -0.5*reciproke_visc*fac, shp);
                assembler.template Matrix<Tauzz,Pres>(shp_tau, -0.5*reciproke_visc*fac, shp);
                
                assembler.template Vector<Tauxx>(shp_tau, 0.5*reciproke_visc*fac*press);
                assembler.template Vector<Tauyy>(shp_tau, 0.5*reciproke_visc*fac*press);
                assembler.template Vector<Tauzz>(shp_tau, 0.5*reciproke_visc*fac*press);
            }
            
            //----------------------------------------------------------------------
            //                 PRESSURE STABILISATION PART
            if(pstab)
            {
                /* pressure stabilisation: convection, convective part */
                /*
                          /                             \
                         |             / n+1       \     |
                         | nabla q ,  | u   o nabla | Du |
                         |             \ i         /     |
                          \                             /
                */
                assembler.template Matrix<Pres,Velx>(shp_dx, tau_Mp, enr_conv_c_);
                assembler.template Matrix<Pres,Vely>(shp_dy, tau_Mp, enr_conv_c_);
                assembler.template Matrix<Pres,Velz>(shp_dz, tau_Mp, enr_conv_c_);
                
                if (newton)
                {
                    /*  pressure stabilisation: convection, reactive part
                          /                             \
                         |           /          \   n+1  |
                         | grad q , | Du o nabla | u     |
                         |           \          /   (i)  |
                          \                             /
                    */
                    assembler.template Matrix<Pres,Velx>(shp_dx, tau_Mp*vderxy_(0,0), shp);
                    assembler.template Matrix<Pres,Velx>(shp_dy, tau_Mp*vderxy_(1,0), shp);
                    assembler.template Matrix<Pres,Velx>(shp_dz, tau_Mp*vderxy_(2,0), shp);
                    
                    assembler.template Matrix<Pres,Vely>(shp_dx, tau_Mp*vderxy_(0,1), shp);
                    assembler.template Matrix<Pres,Vely>(shp_dy, tau_Mp*vderxy_(1,1), shp);
                    assembler.template Matrix<Pres,Vely>(shp_dz, tau_Mp*vderxy_(2,1), shp);
                    
                    assembler.template Matrix<Pres,Velz>(shp_dx, tau_Mp*vderxy_(0,2), shp);
                    assembler.template Matrix<Pres,Velz>(shp_dy, tau_Mp*vderxy_(1,2), shp);
                    assembler.template Matrix<Pres,Velz>(shp_dz, tau_Mp*vderxy_(2,2), shp);
                }
                
//                assembler.template Vector<Pres>(shp_dx, -tau_Mp*conv_old_(0));
//                assembler.template Vector<Pres>(shp_dy, -tau_Mp*conv_old_(1));
//                assembler.template Vector<Pres>(shp_dz, -tau_Mp*conv_old_(2));
                
                /* pressure stabilisation: viscosity (-L_visc_u) */
                /*
                           /                             \
                          |                         /  \  |
                        - |  nabla q , nabla o tau | Du | |
                          |                         \  /  |
                           \                             /
                */
                assembler.template Matrix<Pres,Velx>(shp_dx, -2.0*visc*tau_Mp, enr_viscs2_(0, 0, _));
                assembler.template Matrix<Pres,Vely>(shp_dx, -2.0*visc*tau_Mp, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Pres,Velz>(shp_dx, -2.0*visc*tau_Mp, enr_viscs2_(0, 2, _));
                
                assembler.template Matrix<Pres,Velx>(shp_dy, -2.0*visc*tau_Mp, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Pres,Vely>(shp_dy, -2.0*visc*tau_Mp, enr_viscs2_(1, 1, _));
                assembler.template Matrix<Pres,Velz>(shp_dy, -2.0*visc*tau_Mp, enr_viscs2_(1, 2, _));
                
                assembler.template Matrix<Pres,Velx>(shp_dz, -2.0*visc*tau_Mp, enr_viscs2_(0, 2, _));
                assembler.template Matrix<Pres,Vely>(shp_dz, -2.0*visc*tau_Mp, enr_viscs2_(1, 2, _));
                assembler.template Matrix<Pres,Velz>(shp_dz, -2.0*visc*tau_Mp, enr_viscs2_(2, 2, _));
                
//                assembler.template Vector<Pres>(shp_dx, tau_Mp*2.0*visc*visc_old_(0));
//                assembler.template Vector<Pres>(shp_dy, tau_Mp*2.0*visc*visc_old_(1));
//                assembler.template Vector<Pres>(shp_dz, tau_Mp*2.0*visc*visc_old_(2));
                      
                /* pressure stabilisation: pressure( L_pres_p) */
                /*
                          /                    \
                         |                      |
                         |  nabla q , nabla Dp  |
                         |                      |
                          \                    /
                */
                assembler.template Matrix<Pres,Pres>(shp_dx, tau_Mp, shp_dx);
                assembler.template Matrix<Pres,Pres>(shp_dy, tau_Mp, shp_dy);
                assembler.template Matrix<Pres,Pres>(shp_dz, tau_Mp, shp_dz);
                
//                assembler.template Vector<Pres>(shp_dx, -tau_Mp*gradp_(0));
//                assembler.template Vector<Pres>(shp_dy, -tau_Mp*gradp_(1));
//                assembler.template Vector<Pres>(shp_dz, -tau_Mp*gradp_(2));
                
                assembler.template Vector<Pres>(shp_dx, -tau_Mp*res_old(0));
                assembler.template Vector<Pres>(shp_dy, -tau_Mp*res_old(1));
                assembler.template Vector<Pres>(shp_dz, -tau_Mp*res_old(2));
            }
            
            //----------------------------------------------------------------------
            //                     SUPG STABILISATION PART
            if(supg)
            {
                /* supg stabilisation: convective part ( L_conv_u) */
                /*
                     /                                          \
                    |  / n+1        \        / n+1        \      |
                    | | u    o nabla | v ,  | u    o nabla | Du  |
                    |  \ (i)        /        \ (i)        /      |
                     \                                          /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, tau_M, enr_conv_c_);
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, tau_M, enr_conv_c_);
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, tau_M, enr_conv_c_);
                /* supg stabilisation: pressure part  ( L_pres_p) */
                /*
                          /                             \
                         |   / n+1       \               |
                         |  | u   o nabla | v , nabla Dp |
                         |   \ (i)       /               |
                          \                             /
                */
                assembler.template Matrix<Velx,Pres>(enr_conv_c_, tau_M, shp_dx);
                assembler.template Matrix<Vely,Pres>(enr_conv_c_, tau_M, shp_dy);
                assembler.template Matrix<Velz,Pres>(enr_conv_c_, tau_M, shp_dz);

                /* supg stabilisation: viscous part  (-L_visc_u) */
                /*
                      /                                        \
                     |               /  \    / n+1        \     |
                     |  nabla o eps | Du |, | u    o nabla | v  |
                     |               \  /    \ (i)        /     |
                      \                                        /
                */
                assembler.template Matrix<Velx,Velx>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(0, 0, _));
                assembler.template Matrix<Velx,Vely>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Velx,Velz>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(0, 2, _));
    
                assembler.template Matrix<Vely,Velx>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(0, 1, _));
                assembler.template Matrix<Vely,Vely>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(1, 1, _));
                assembler.template Matrix<Vely,Velz>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(1, 2, _));
    
                assembler.template Matrix<Velz,Velx>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(0, 2, _));
                assembler.template Matrix<Velz,Vely>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(1, 2, _));
                assembler.template Matrix<Velz,Velz>(enr_conv_c_, -2.0*visc*tau_M, enr_viscs2_(2, 2, _));
                
                if (newton)
                {
                    /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                    /*
                               /                                           \
                              |    /          \   n+1    / n+1        \     |
                              |   | Du o nabla | u    , | u    o nabla | v  |
                              |    \          /   (i)    \ (i)        /     |
                               \                                           /
                    */
                    assembler.template Matrix<Velx,Velx>(enr_conv_c_, tau_M*vderxy_(0,0), shp);
                    assembler.template Matrix<Velx,Vely>(enr_conv_c_, tau_M*vderxy_(0,1), shp);
                    assembler.template Matrix<Velx,Velz>(enr_conv_c_, tau_M*vderxy_(0,2), shp);
                    
                    assembler.template Matrix<Vely,Velx>(enr_conv_c_, tau_M*vderxy_(1,0), shp);
                    assembler.template Matrix<Vely,Vely>(enr_conv_c_, tau_M*vderxy_(1,1), shp);                    
                    assembler.template Matrix<Vely,Velz>(enr_conv_c_, tau_M*vderxy_(1,2), shp);                    

                    assembler.template Matrix<Velz,Velx>(enr_conv_c_, tau_M*vderxy_(2,0), shp);
                    assembler.template Matrix<Velz,Vely>(enr_conv_c_, tau_M*vderxy_(2,1), shp);
                    assembler.template Matrix<Velz,Velz>(enr_conv_c_, tau_M*vderxy_(2,2), shp);
                    
                    /*
                             /                                           \
                            |    / n+1        \   n+1    /          \     |
                            |   | u    o nabla | u    , | Du o nabla | v  |
                            |    \ (i)        /   (i)    \          /     |
                             \                                           /
                    */
                    assembler.template Matrix<Velx,Velx>(shp_dx, tau_M*velint_(0)*vderxy_(0,0), shp);
                    assembler.template Matrix<Velx,Velx>(shp_dx, tau_M*velint_(1)*vderxy_(0,1), shp);
                    assembler.template Matrix<Velx,Velx>(shp_dx, tau_M*velint_(2)*vderxy_(0,2), shp);
                    
                    assembler.template Matrix<Velx,Vely>(shp_dy, tau_M*velint_(0)*vderxy_(0,0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, tau_M*velint_(1)*vderxy_(0,1), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, tau_M*velint_(2)*vderxy_(0,2), shp);                        
                    
                    assembler.template Matrix<Velx,Velz>(shp_dz, tau_M*velint_(0)*vderxy_(0,0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, tau_M*velint_(1)*vderxy_(0,1), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, tau_M*velint_(2)*vderxy_(0,2), shp); 
                    
                    
                    assembler.template Matrix<Vely,Velx>(shp_dx, tau_M*velint_(0)*vderxy_(1,0), shp);
                    assembler.template Matrix<Vely,Velx>(shp_dx, tau_M*velint_(1)*vderxy_(1,1), shp);
                    assembler.template Matrix<Vely,Velx>(shp_dx, tau_M*velint_(2)*vderxy_(1,2), shp);
                    
                    assembler.template Matrix<Vely,Vely>(shp_dy, tau_M*velint_(0)*vderxy_(1,0), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, tau_M*velint_(1)*vderxy_(1,1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, tau_M*velint_(2)*vderxy_(1,2), shp);
                    
                    assembler.template Matrix<Vely,Velz>(shp_dz, tau_M*velint_(0)*vderxy_(1,0), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, tau_M*velint_(1)*vderxy_(1,1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, tau_M*velint_(2)*vderxy_(1,2), shp);
                    
                    
                    assembler.template Matrix<Velz,Velx>(shp_dx, tau_M*velint_(0)*vderxy_(2,0), shp);
                    assembler.template Matrix<Velz,Velx>(shp_dx, tau_M*velint_(1)*vderxy_(2,1), shp);
                    assembler.template Matrix<Velz,Velx>(shp_dx, tau_M*velint_(2)*vderxy_(2,2), shp);
                    
                    assembler.template Matrix<Velz,Vely>(shp_dy, tau_M*velint_(0)*vderxy_(2,0), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, tau_M*velint_(1)*vderxy_(2,1), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, tau_M*velint_(2)*vderxy_(2,2), shp);
                    
                    assembler.template Matrix<Velz,Velz>(shp_dz, tau_M*velint_(0)*vderxy_(2,0), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, tau_M*velint_(1)*vderxy_(2,1), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, tau_M*velint_(2)*vderxy_(2,2), shp);


                    /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
                    /*
                                    /                               \
                                   |         n+1    /          \     |
                                   |  nabla p    , | Du o nabla | v  |
                                   |         (i)    \          /     |
                                    \                               /
                    */
                    assembler.template Matrix<Velx,Velx>(shp_dx, tau_M*gradp_(0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, tau_M*gradp_(0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, tau_M*gradp_(0), shp);

                    assembler.template Matrix<Vely,Velx>(shp_dx, tau_M*gradp_(1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, tau_M*gradp_(1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, tau_M*gradp_(1), shp);
                        
                    assembler.template Matrix<Velz,Velx>(shp_dx, tau_M*gradp_(2), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, tau_M*gradp_(2), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, tau_M*gradp_(2), shp);


                      /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                      /*
                              /                                         \
                             |               / n+1 \    /          \     |
                             |  nabla o eps | u     |, | Du o nabla | v  |
                             |               \ (i) /    \          /     |
                              \                                         /
                      */
                    assembler.template Matrix<Velx,Velx>(shp_dx, -2.0*visc*tau_M*visc_old_(0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, -2.0*visc*tau_M*visc_old_(0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, -2.0*visc*tau_M*visc_old_(0), shp);

                    assembler.template Matrix<Vely,Velx>(shp_dx, -2.0*visc*tau_M*visc_old_(1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, -2.0*visc*tau_M*visc_old_(1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, -2.0*visc*tau_M*visc_old_(1), shp);
                        
                    assembler.template Matrix<Velz,Velx>(shp_dx, -2.0*visc*tau_M*visc_old_(2), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, -2.0*visc*tau_M*visc_old_(2), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, -2.0*visc*tau_M*visc_old_(2), shp);

                    /* supg stabilisation: bodyforce part, linearisation of test function */

                    /*
                                  /                             \
                                 |              /          \     |
                                 |  rhsint   , | Du o nabla | v  |
                                 |              \          /     |
                                  \                             /

                    */
                    assembler.template Matrix<Velx,Velx>(shp_dx, -tau_M*rhsint_(0), shp);
                    assembler.template Matrix<Velx,Vely>(shp_dy, -tau_M*rhsint_(0), shp);
                    assembler.template Matrix<Velx,Velz>(shp_dz, -tau_M*rhsint_(0), shp);
                    
                    assembler.template Matrix<Vely,Velx>(shp_dx, -tau_M*rhsint_(1), shp);
                    assembler.template Matrix<Vely,Vely>(shp_dy, -tau_M*rhsint_(1), shp);
                    assembler.template Matrix<Vely,Velz>(shp_dz, -tau_M*rhsint_(1), shp);
                    
                    assembler.template Matrix<Velz,Velx>(shp_dx, -tau_M*rhsint_(2), shp);
                    assembler.template Matrix<Velz,Vely>(shp_dy, -tau_M*rhsint_(2), shp);
                    assembler.template Matrix<Velz,Velz>(shp_dz, -tau_M*rhsint_(2), shp);
                } // if newton
                
                // supg stabilisation
                assembler.template Vector<Velx>(enr_conv_c_, -tau_M*res_old(0));
                assembler.template Vector<Vely>(enr_conv_c_, -tau_M*res_old(1));
                assembler.template Vector<Velz>(enr_conv_c_, -tau_M*res_old(2));
            }
            
            
            //----------------------------------------------------------------------
            //                       STABILISATION, VISCOUS PART
            if(vstab)
            {
                //
                //  viscous part of stabilisation is switched off!
                //  vstab is set to false within fluid3_evaluate.cpp
                //
                dserror("not translated to XFEM assembler yet!");
            };
            
            //----------------------------------------------------------------------
            //                     STABILISATION, CONTINUITY PART
            if(cstab)
            {
                /* continuity stabilisation */
                /*
                         /                        \
                        |                          |
                        | nabla o Du  , nabla o v  |
                        |                          |
                         \                        /
                */
                assembler.template Matrix<Velx,Velx>(shp_dx, tau_C, shp_dx);
                assembler.template Matrix<Velx,Vely>(shp_dx, tau_C, shp_dy);
                assembler.template Matrix<Velx,Velz>(shp_dx, tau_C, shp_dz);
                
                assembler.template Matrix<Vely,Velx>(shp_dy, tau_C, shp_dx);
                assembler.template Matrix<Vely,Vely>(shp_dy, tau_C, shp_dy);
                assembler.template Matrix<Vely,Velz>(shp_dy, tau_C, shp_dz);
                
                assembler.template Matrix<Velz,Velx>(shp_dz, tau_C, shp_dx);
                assembler.template Matrix<Velz,Vely>(shp_dz, tau_C, shp_dy);
                assembler.template Matrix<Velz,Velz>(shp_dz, tau_C, shp_dz);
                
                const double tau_C_divunp = tau_C*(vderxy_(0, 0)+vderxy_(1, 1)+vderxy_(2, 2));
                assembler.template Vector<Velx>(shp_dx, -tau_C_divunp);
                assembler.template Vector<Vely>(shp_dy, -tau_C_divunp);
                assembler.template Vector<Velz>(shp_dz, -tau_C_divunp);
            } // endif cstab
        } // end loop over gauss points
    } // end loop over integration cells
    
    
    
    //if (false)
    if (ASSTYPE == XFEM::xfem_assembly)
    {
        // for now, I don't try to compare to elements wiout stress unknowns, since they lock anyway
        if (stress_unknowns_present){
    
    // information about boundary integration cells
    const XFEM::BoundaryIntCells boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());
    
    
    
    bool printcell = true;
    
    if (printcell)
    {
        printcell = false;
        std::ofstream f_system("example.pos");
        
        
        {
            // draw elements with associated gid
            stringstream gmshfilecontent;
            gmshfilecontent << "View \" " << "DomainIntegrationCells \" {" << endl;
            XFEM::DomainIntCells::const_iterator cell;
            for (cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
            {
                //gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXYZ(*ele), 0.0, cell->Shape()) << endl;
                gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXiDomain(), 0.0, cell->Shape()) << endl;
            };
            gmshfilecontent << "};" << endl;
            f_system << gmshfilecontent.str();
        }
        {
            // draw elements with associated gid
            stringstream gmshfilecontent;
            gmshfilecontent << "View \" " << "BoundaryIntegrationCells \" {" << endl;
            XFEM::BoundaryIntCells::const_iterator cell;
            for (cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
            {
                //gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXYZ(*ele), 0.0, cell->Shape()) << endl;
                gmshfilecontent << IO::GMSH::cellToString(cell->NodalPosXiDomain(), 0.0, cell->Shape()) << endl;
            };
            gmshfilecontent << "};" << endl;
            f_system << gmshfilecontent.str();
        }
        f_system.close();

    }
    
    
    
    
    
    // loop over boundary integration cells
    for (XFEM::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
    {

        // gaussian points
        const DRT::UTILS::IntegrationPoints2D intpoints(DRT::UTILS::intrule_tri_6point);
        
        // get node coordinates
        //TODO: get the right element
        //const RCP<DRT::Element> boundaryele = ih->GetBoundaryEle(cell->GetSurfaceEleGid());
        const DRT::Element* boundaryele = ih->GetBoundaryEle(0);
        //cout << (*boundaryele) << endl;
        
        const BlitzMat xyze_boundary(DRT::UTILS::PositionArray(boundaryele));

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta^\boundary
            XFEM::PosEtaBoundary pos_eta_boundary(2);
            pos_eta_boundary(0) = intpoints.qxg[iquad][0];
            pos_eta_boundary(1) = intpoints.qxg[iquad][1];
      
            // coordinates of the current integration point in element coordinates \xi
            const XFEM::PosXiBoundary posXiBoundary(mapEtaBToXiB(*cell, pos_eta_boundary));
            const XFEM::PosXiDomain posXiDomain(mapEtaBToXiD(*cell, pos_eta_boundary));
            const double detcell = abs(detEtaBToXiB(*cell, pos_eta_boundary));

            // shape functions and their first derivatives
            const BlitzVec funct_boundary(DRT::UTILS::shape_function_2D(posXiBoundary(0),posXiBoundary(1),boundaryele->Shape()));
            const BlitzMat deriv_boundary(DRT::UTILS::shape_function_2D_deriv1(posXiBoundary(0),posXiBoundary(1),boundaryele->Shape()));
            
            // shape functions and their first derivatives
            const BlitzVec funct_(DRT::UTILS::shape_function_3D(posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE));
            const BlitzMat deriv_(DRT::UTILS::shape_function_3D_deriv1(posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE));
      
            const BlitzVec funct_stress(DRT::UTILS::shape_function_3D(posXiDomain(0),posXiDomain(1),posXiDomain(2),stressdistype));

            // position of the gausspoint in physical coordinates
            const BlitzVec gauss_pos_xyz(blitz::sum(funct_boundary(j)*xyze_boundary(i,j),j));
      
            // get jacobian matrix d x / d \xi
            const BlitzMat metric = computeMetricTensor(xyze_boundary,deriv_boundary);
            
            const double detmetric = sqrt(metric(0,0)*metric(1,1) - metric(0,1)*metric(1,0));
            const double fac = intpoints.qwgt[iquad]*detmetric*detcell;

            if (detmetric < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %lf", ele->Id(), detmetric);
            }

            // after this call, one should only use the enriched shape functions and derivatives!
            BlitzVec enr_funct_(numparamvelx);
            BlitzVec enr_funct_stress(numparamtauxx);
            
            // shape function for nodal dofs
            XFEM::ComputeEnrichedShapefunction(
                    *ele,
                    ih,
                    dofman,
                    Velx,
                    gauss_pos_xyz,
                    XFEM::Enrichment::approachFromPlus,
                    funct_,
                    enr_funct_);

            // shape functions for element dofs
            XFEM::ComputeEnrichedStressShapefunction(
                    *ele,
                    ih,
                    dofman,
                    Tauxx,
                    gauss_pos_xyz,
                    XFEM::Enrichment::approachFromPlus,
                    funct_stress,
                    enr_funct_stress);
                
            // perform integration for entire matrix and rhs
            // create vievs on the shape function arrays for easy handling in the assembly process
            const BlitzVec shp(enr_funct_);
            const BlitzVec shp_tau(enr_funct_stress);
            
            // get normal vector (in x coordinates) to surface element at integration point
            const BlitzVec normalvec(computeNormalToBoundaryElement(deriv_boundary, xyze_boundary, posXiBoundary(0), posXiBoundary(1)));
      
            // get velocities (n+g,i) at integration point
            const BlitzVec velint_(blitz::sum(shp(j)*evelnp(i,j),j));

            // get viscous stress unknowns
            BlitzMat tau(3,3,blitz::ColumnMajorArray<2>());
            tau(0,0) = blitz::sum(shp_tau(_)*etau(0,_)); tau(0,1) = blitz::sum(shp_tau(_)*etau(3,_)); tau(0,2) = blitz::sum(shp_tau(_)*etau(4,_));
            tau(1,0) = tau(0,1);                         tau(1,1) = blitz::sum(shp_tau(_)*etau(1,_)); tau(1,2) = blitz::sum(shp_tau(_)*etau(5,_));
            tau(2,0) = tau(0,2);                         tau(2,1) = tau(1,2);                         tau(2,2) = blitz::sum(shp_tau(_)*etau(2,_));



            /*--------------------------------- now build single stiffness terms ---*/

            /*
                           /                   \
                          |                     |
                          | (virt tau) * n , Du |
                          |                     |
                           \                   /
            */
            // i = x
            assembler.template Matrix<Tauxx,Velx>(shp_tau, fac*normalvec(0), shp);
            assembler.template Matrix<Tauxy,Velx>(shp_tau, fac*normalvec(1), shp);
            assembler.template Matrix<Tauxz,Velx>(shp_tau, fac*normalvec(2), shp);
            // i = y
            assembler.template Matrix<Tauyx,Vely>(shp_tau, fac*normalvec(0), shp);
            assembler.template Matrix<Tauyy,Vely>(shp_tau, fac*normalvec(1), shp);
            assembler.template Matrix<Tauyz,Vely>(shp_tau, fac*normalvec(2), shp);
            // i = z
            assembler.template Matrix<Tauzx,Velz>(shp_tau, fac*normalvec(0), shp);
            assembler.template Matrix<Tauzy,Velz>(shp_tau, fac*normalvec(1), shp);
            assembler.template Matrix<Tauzz,Velz>(shp_tau, fac*normalvec(2), shp);
            
            // i = x
            assembler.template Vector<Tauxx>(shp_tau, -fac*normalvec(0)*velint_(0));
            assembler.template Vector<Tauxy>(shp_tau, -fac*normalvec(1)*velint_(0));
            assembler.template Vector<Tauxz>(shp_tau, -fac*normalvec(2)*velint_(0));
            // i = y
            assembler.template Vector<Tauyx>(shp_tau, -fac*normalvec(0)*velint_(1));
            assembler.template Vector<Tauyy>(shp_tau, -fac*normalvec(1)*velint_(1));
            assembler.template Vector<Tauyz>(shp_tau, -fac*normalvec(2)*velint_(1));
            // i = z            
            assembler.template Vector<Tauzx>(shp_tau, -fac*normalvec(0)*velint_(2));
            assembler.template Vector<Tauzy>(shp_tau, -fac*normalvec(1)*velint_(2));
            assembler.template Vector<Tauzz>(shp_tau, -fac*normalvec(2)*velint_(2));
            
            // symmetric part from partial integration of stress term
            /*
                           /           \
                          |             |
                        - | v , tau * n |
                          |             |
                           \           /
            */
            assembler.template Matrix<Velx,Tauxx>(shp, fac*normalvec(0), shp_tau);
            assembler.template Matrix<Velx,Tauxy>(shp, fac*normalvec(1), shp_tau);
            assembler.template Matrix<Velx,Tauxz>(shp, fac*normalvec(2), shp_tau);
            
            assembler.template Matrix<Vely,Tauyx>(shp, fac*normalvec(0), shp_tau);
            assembler.template Matrix<Vely,Tauyy>(shp, fac*normalvec(1), shp_tau);
            assembler.template Matrix<Vely,Tauyz>(shp, fac*normalvec(2), shp_tau);
            
            assembler.template Matrix<Velz,Tauzx>(shp, fac*normalvec(0), shp_tau);
            assembler.template Matrix<Velz,Tauzy>(shp, fac*normalvec(1), shp_tau);
            assembler.template Matrix<Velz,Tauzz>(shp, fac*normalvec(2), shp_tau);
            
            const BlitzVec tau_times_n(blitz::sum(tau(i,j)*normalvec(j),j));
            assembler.template Vector<Velx>(shp, -fac*tau_times_n(0));
            assembler.template Vector<Vely>(shp, -fac*tau_times_n(1));
            assembler.template Vector<Velz>(shp, -fac*tau_times_n(2));
            
//            /*
//                           /         \
//                          |           |
//                          | v , Dp n  |
//                          |           |
//                           \         /
//            */
//            assembler.template Matrix<Velx,Pres>(shp, -fac*normalvec(0), shp);
//            assembler.template Matrix<Vely,Pres>(shp, -fac*normalvec(1), shp);
//            assembler.template Matrix<Velz,Pres>(shp, -fac*normalvec(2), shp);
//            
//            assembler.template Vector<Velx>(shp, fac*normalvec(0)*press);
//            assembler.template Vector<Vely>(shp, fac*normalvec(1)*press);
//            assembler.template Vector<Velz>(shp, fac*normalvec(2)*press);  

        } // end loop over gauss points
    } // end loop over boundary integration cells
    }
    } // if (ASSTYPE == XFEM::xfem_assembly)
    return;
}
  
void callSysmat(
        const XFEM::AssemblyType assembly_type,
        const DRT::ELEMENTS::XFluid3* ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&    eleDofManager,
        const std::vector<double>&        locval,
        blitz::Array<double,2>&           estif,
        blitz::Array<double,1>&           eforce,
        struct _MATERIAL*       material,
        double                  pseudotime,
        bool                    newton ,
        bool                    pstab  ,
        bool                    supg   ,
        bool                    vstab  ,
        bool                    cstab
        )
{
    using DRT::Element;

    if (assembly_type == XFEM::standard_assembly)
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, estif, eforce,
                        material, pseudotime, newton, pstab, supg, vstab, cstab);
                break;
            case DRT::Element::hex20:
                Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, estif, eforce,
                        material, pseudotime, newton, pstab, supg, vstab, cstab);
                break;
            case DRT::Element::hex27:
                Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
                        ele, ih, eleDofManager, locval, estif, eforce,
                        material, pseudotime, newton, pstab, supg, vstab, cstab);
                break;
            default:
                dserror("not templated yet");
        };
    }
    else
    {
        switch (ele->Shape())
        {
            case DRT::Element::hex8:
                Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, locval, estif, eforce,
                        material, pseudotime, newton, pstab, supg, vstab, cstab);
                break;
            case DRT::Element::hex20:
                Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, locval, estif, eforce,
                        material, pseudotime, newton, pstab, supg, vstab, cstab);
                break;
            case DRT::Element::hex27:
                Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
                        ele, ih, eleDofManager, locval, estif, eforce,
                        material, pseudotime, newton, pstab, supg, vstab, cstab);
                break;
            default:
                dserror("not templated yet");
        };
    }
}
} // end namespace XFLUID

#endif

#endif
#endif
