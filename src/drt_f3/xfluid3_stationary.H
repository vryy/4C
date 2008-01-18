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

#include "xfluid3.H"
#include <blitz/array.h>
#include "../drt_mat/newtonianfluid.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/xfem.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_lib/drt_utils.H"
#include "fluid3_stabilization.H"
#include "xfluid3_local_assembler.H"
#include "xfluid3_interpolation.H"

#include <Epetra_SerialDenseSolver.h>

using namespace std;
using namespace XFEM::PHYSICS;


typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix


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
  const int numnode = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
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
        curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(pseudotime);
      }
      else
      {
    // do not compute an "alternative" curvefac here since a negative pseudotime value
    // indicates an error.
        dserror("Negative pseudotime value in body force calculation: time = %f",pseudotime);
        //curvefac = DRT::Utils::TimeCurveManager::Instance().Curve(curvenum).f(0.0);
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
  template <DRT::Element::DiscretizationType DISTYPE>
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
    const int numnode = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
    Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,numnode);

    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(bm);

    // No need for a separate rhs. We assemble the rhs to the solution
    // vector. The solver will destroy the rhs and return the solution.
    solver.SetVectors(ederxy2,ederxy2);
    solver.Solve();

    return derxy2_;
  }
  

  // get volume in parameter space using Gauss integration
  template <XFEM::AssemblyType ASSTYPE>
  vector<double> modifyGaussRule3D(
          const XFEM::DomainIntCell& cell,
          const double               cell_e0,
          const double               cell_e1,
          const double               cell_e2)
  {   
      // return value
      vector<double> element_e(4);
      if (ASSTYPE == XFEM::standard_assembly) {
          // gauss coordinates of cell in element coordinates
          element_e[0] = cell_e0;
          element_e[1] = cell_e1;
          element_e[2] = cell_e2;
          element_e[3] = 1.0;
      } else {
      
          const DRT::Element::DiscretizationType celldistype = cell.Shape();
          const int numnode = DRT::Utils::getNumberOfElementNodes(celldistype);
          const int nsd = 3;
      
          // get node coordinates
          blitz::Array<double,2> xyze_cell(nsd,numnode,blitz::ColumnMajorArray<2>());
          for (int inode=0; inode<numnode; inode++)
              {
              xyze_cell(0,inode) = cell.GetDomainCoord()[inode][0];
              xyze_cell(1,inode) = cell.GetDomainCoord()[inode][1];
              xyze_cell(2,inode) = cell.GetDomainCoord()[inode][2];
              }    

          // init blitz indices
          blitz::firstIndex i;    // Placeholder for the first index
          blitz::secondIndex j;   // Placeholder for the second index
          blitz::thirdIndex k;    // Placeholder for the third index
      
          // create shape function vectors 
          blitz::Array<double,1> funct(numnode);
          blitz::Array<double,2> deriv(nsd,numnode,blitz::ColumnMajorArray<2>());
          DRT::Utils::shape_function_3D(funct,cell_e0,cell_e1,cell_e2,celldistype);
          DRT::Utils::shape_function_3D_deriv1(deriv,cell_e0,cell_e1,cell_e2,celldistype);

          // translate position into from cell coordinates to element coordinates
          const blitz::Array<double,1> e(blitz::sum(funct(j)*xyze_cell(i,j),j));


          // get Jacobian matrix and determinant
          // actually compute its transpose....
      /*
        +-            -+ T      +-            -+
        | dx   dx   dx |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dr   dr   dr |
        |              |        |              |
        | dy   dy   dy |        | dx   dy   dz |
        | --   --   -- |   =    | --   --   -- |
        | dr   ds   dt |        | ds   ds   ds |
        |              |        |              |
        | dz   dz   dz |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dt   dt   dt |
        +-            -+        +-            -+
      */
          const blitz::Array<double,2> xjm(blitz::sum(deriv(i,k)*xyze_cell(j,k),k));
          const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                             xjm(0,1)*xjm(1,2)*xjm(2,0)+
                             xjm(0,2)*xjm(1,0)*xjm(2,1)-
                             xjm(0,2)*xjm(1,1)*xjm(2,0)-
                             xjm(0,0)*xjm(1,2)*xjm(2,1)-
                             xjm(0,1)*xjm(1,0)*xjm(2,2);
      
    
          // gauss coordinates of cell in element coordinates + gauss weight
          element_e[0] = e(0);
          element_e[1] = e(1);
          element_e[2] = e(2);
          element_e[3] = det;
      }  
      return element_e;
  };
  
  
  template<XFEM::AssemblyType>
  int getNumParam(
          const XFEM::ElementDofManager& dofman,
          const XFEM::PHYSICS::Field     field,
          int numnode)
  {
      return dofman.NumDofPerField(field);
  }
  template<> int getNumParam<XFEM::standard_assembly>(
          const XFEM::ElementDofManager& dofman,
          const XFEM::PHYSICS::Field     field,
          int numnode)
  {
      return numnode;
  }
  

  
  /*!
    Calculate matrix and rhs for stationary problem formulation

   */
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE>
  void Sysmat(
          DRT::Elements::XFluid3*           ele,           ///< the element those matrix is calculated
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
    const int numnode = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
      
    // information about domain integration cells
    const XFEM::DomainIntCells   domainIntCells   = ih->GetDomainIntCells(ele->Id(),DISTYPE);
    // information about boundary integration cells
    const XFEM::BoundaryIntCells boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

    // get node coordinates
    DRT::Node** const nodes = ele->Nodes();
    BlitzMat xyze_(3,numnode,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<numnode; inode++)
    {
      const double* x = nodes[inode]->X();
      xyze_(0,inode) = x[0];
      xyze_(1,inode) = x[1];
      xyze_(2,inode) = x[2];
    }

    // dead load in element nodes
    //////////////////////////////////////////////////// , BlitzMat edeadng_(BodyForce(ele->Nodes(),pseudotime));

    // get viscosity
    // check here, if we really have a fluid !!
    dsassert(material->mattyp == m_fluid, "Material law is not of type m_fluid.");
    const double visc = material->m.fluid->viscosity;
    const double reciproke_visc = 1.0 / visc;

    // flag for higher order elements
    const bool higher_order_ele = ele->isHigherOrderElement<DISTYPE>();

    // integrate of integration cell
    dsassert(domainIntCells.empty() == false, "this is a bug!");

    DRT::Utils::GaussRule3D gaussrule;
    if (domainIntCells.size() == 1)
    {
        gaussrule = ele->GetGaussrule();
    }
    else
    {
        // attention: don't use gausspoint rules that have gausspoints on the element boundary 
        // (so called closed rules)
        // -> enrichment value might be undefined 
        gaussrule = DRT::Utils::intrule_tet_11point;
    }
    
    LocalAssembler<DISTYPE, ASSTYPE> assembler(dofman, estif, eforce);
    
    const blitz::Range _  = blitz::Range::all();
    
    // We define the variables i,j,k to be indices to blitz arrays.
    // These are used for array expressions, that is matrix-vector
    // products in the following.

    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
    const int numparamvely = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Vely, numnode);
    const int numparamvelz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velz, numnode);
    const int numparampres = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Pres, numnode);
    const int numparamtauxx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 0);
    const int numparamtauyy = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyy, 0);
    const int numparamtauzz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauzz, 0);
    const int numparamtauxy = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxy, 0);
    const int numparamtauxz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxz, 0);
    const int numparamtauyz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyz, 0);
//          cout << "numparampres " << numparampres << endl;
//          cout << "numparamvelx " << numparamvelx << endl;
//          cout << "numparamvely " << numparamvely << endl;
//          cout << "numparamvelz " << numparamvelz << endl;
    
    dsassert((numparamvelx == numparamvely and numparamvelx == numparamvelz and numparamvelx == numparampres),
            "for now, we enrich velocity and pressure together");
    
    const vector<int> velxdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Velx);
    const vector<int> velydof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Vely);
    const vector<int> velzdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Velz);
    const vector<int> presdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Pres);
    const vector<int> tauxxdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxx);
    const vector<int> tauyydof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauyy);
    const vector<int> tauzzdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauzz);
    const vector<int> tauxydof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxy);
    const vector<int> tauxzdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxz);
    const vector<int> tauyzdof = dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauyz);

    // split velocity and pressure
    blitz::Array<double, 1> eprenp(numparampres);
    blitz::Array<double, 2> evelnp(3,numparamvelx,blitz::ColumnMajorArray<2>());
    blitz::Array<double, 2> etau(6,numparamtauxx,blitz::ColumnMajorArray<2>());
    
    for (int iparam=0; iparam<numparamvelx; ++iparam)   evelnp(0,iparam) = locval[velxdof[iparam]];
    for (int iparam=0; iparam<numparamvely; ++iparam)   evelnp(1,iparam) = locval[velydof[iparam]];
    for (int iparam=0; iparam<numparamvelz; ++iparam)   evelnp(2,iparam) = locval[velzdof[iparam]];
    for (int iparam=0; iparam<numparampres; ++iparam)   eprenp(iparam) = locval[presdof[iparam]];
    for (int iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = locval[tauxxdof[iparam]];
    for (int iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = locval[tauyydof[iparam]];
    for (int iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = locval[tauzzdof[iparam]];
    for (int iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = locval[tauxydof[iparam]];
    for (int iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = locval[tauxzdof[iparam]];
    for (int iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = locval[tauyzdof[iparam]];
    
    
    
    // stabilization parameter
    // This has to be done before anything else is calculated because
    // we use the same arrays internally.
    const double hk = XFLUID::HK<DISTYPE>(evelnp,visc,xyze_);
    const double mk = XFLUID::MK<DISTYPE>();
    
    // loop over integration cells
    for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell < domainIntCells.end(); ++cell)
    {

    // gaussian points
    const DRT::Utils::IntegrationPoints3D intpoints(gaussrule);
    
    // if cell node is on the interface, the value is not defined for a jump.
    // however, we approach the interface from one particular side and therefore,
    // -> we use the center of the cell to determine, where we come from
    //const blitz::Array<double,1> cellcenterpos(cell->GetPhysicalCenterPosition(*ele));

    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in cell coordinates \eta
      const double cell_e0 = intpoints.qxg[iquad][0];
      const double cell_e1 = intpoints.qxg[iquad][1];
      const double cell_e2 = intpoints.qxg[iquad][2];

      //const vector<double> e = cell->modifyGaussRule3D(standard_integration, cell_e0,cell_e1,cell_e2);
      const vector<double> e = modifyGaussRule3D<ASSTYPE>(*cell, cell_e0,cell_e1,cell_e2);

      // coordinates of the current integration point in element coordinates \xi
      const double e1 = e[0];
      const double e2 = e[1];
      const double e3 = e[2];
      const double detcell = e[3];

      // shape functions and their first derivatives
      BlitzVec funct_(numnode);
      BlitzMat deriv_(3,numnode,blitz::ColumnMajorArray<2>());
      DRT::Utils::shape_function_3D(funct_,e1,e2,e3,DISTYPE);
      DRT::Utils::shape_function_3D_deriv1(deriv_,e1,e2,e3,DISTYPE);
      
      const DRT::Element::DiscretizationType stressdistype = XFLUID::getStressInterpolationType3D(XFLUID::getElementFormulation(DISTYPE));
      BlitzVec funct_stress(DRT::Utils::getNumberOfElementNodes(stressdistype));
      BlitzMat deriv_stress(3, DRT::Utils::getNumberOfElementNodes(stressdistype));
      if (stressdistype != DRT::Element::dis_none)
      {
          //funct_stress = XFLUID::elemental_shape_function_3D(e1,e2,e3,stressdistype);
          //deriv_stress = XFLUID::elemental_shape_function_3D_deriv(e1,e2,e3,stressdistype);
      
          DRT::Utils::shape_function_3D(funct_stress,e1,e2,e3,stressdistype);
          DRT::Utils::shape_function_3D_deriv1(deriv_stress,e1,e2,e3,stressdistype);
      }
      //cout << "funct_stress" << funct_stress << endl;
      //cout << "deriv_stress" << deriv_stress << endl;

      // position of the gausspoint in physical coordinates
      const BlitzVec gauss_pos(blitz::sum(funct_(j)*xyze_(i,j),j));
      
      // get Jacobian matrix and determinant
      // actually compute its transpose....
      /*
        +-            -+ T      +-            -+
        | dx   dx   dx |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dr   dr   dr |
        |              |        |              |
        | dy   dy   dy |        | dx   dy   dz |
        | --   --   -- |   =    | --   --   -- |
        | dr   ds   dt |        | ds   ds   ds |
        |              |        |              |
        | dz   dz   dz |        | dx   dy   dz |
        | --   --   -- |        | --   --   -- |
        | dr   ds   dt |        | dt   dt   dt |
        +-            -+        +-            -+
      */
      // transposed of the jacobian matrix d x / d \xi
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
      BlitzMat xji_(3,3);
      xji_(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
      xji_(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
      xji_(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
      xji_(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
      xji_(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
      xji_(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
      xji_(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
      xji_(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
      xji_(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;

      // compute global derivates
      const BlitzMat derxy_(blitz::sum(xji_(i,k)*deriv_(k,j),k));
      const BlitzMat derxy_stress(blitz::sum(xji_(i,k)*deriv_stress(k,j),k));

      // compute second global derivative
      BlitzMat derxy2_(6,numnode,blitz::ColumnMajorArray<2>());
      if (higher_order_ele)
      {
        BlitzMat deriv2_(6,numnode,blitz::ColumnMajorArray<2>());
        DRT::Utils::shape_function_3D_deriv2(deriv2_,e1,e2,e3,DISTYPE);
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
                  XFEM::PHYSICS::Velx,
                  gauss_pos,
                  funct_,
                  derxy_,
                  derxy2_, 
                  enr_funct_,
                  enr_derxy_,
                  enr_derxy2_);
          
          // shape functions for element dofs
          XFEM::ComputeEnrichedStressShapefunction(
                  *ele,
                  ih,
                  dofman,
                  XFEM::PHYSICS::Tauxx,
                  gauss_pos,
                  funct_stress,
                  derxy_stress,
                  enr_funct_stress,
                  enr_derxy_stress);
      }
      else
      {
          enr_funct_  = funct_;
          enr_derxy_  = derxy_;
          enr_derxy2_ = derxy2_;
          
          enr_funct_stress  = funct_stress;
          enr_derxy_stress  = derxy_stress;
      }
      
//      cout << "enr_funct_stress" << enr_funct_stress << endl;
//      cout << "enr_derxy_stress" << enr_derxy_stress << endl;
      
      // get velocities (n+g,i) at integration point
      const BlitzVec velint_(blitz::sum(enr_funct_(j)*evelnp(i,j),j));

      // get velocity (np,i) derivatives at integration point
      const BlitzMat vderxy_(blitz::sum(enr_derxy_(j,k)*evelnp(i,k),k));
      
      //BlitzMat eps(vderxy_(i,j) + vderxy_(j,i));
      BlitzMat eps(3,3);
      for (int m = 0; m < 3; ++m)
      {
          for (int n = 0; n < 3; ++n)
          {
              eps(m,n) = vderxy_(m,n) + vderxy_(n,m);
          }
      }

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

      // get pressurrrrre
      const double press = blitz::sum(enr_funct_*eprenp);

      // get viscous stress unknowns
      BlitzMat tau(3,3,blitz::ColumnMajorArray<2>());
      tau(0,0) = blitz::sum(enr_funct_stress(_)*etau(0,_)); tau(0,1) = blitz::sum(enr_funct_stress(_)*etau(3,_)); tau(0,2) = blitz::sum(enr_funct_stress(_)*etau(4,_));
      tau(1,0) = tau(0,1);                                  tau(1,1) = blitz::sum(enr_funct_stress(_)*etau(1,_)); tau(1,2) = blitz::sum(enr_funct_stress(_)*etau(5,_));
      tau(2,0) = tau(0,2);                                  tau(2,1) = tau(1,2);                                  tau(2,2) = blitz::sum(enr_funct_stress(_)*etau(2,_));
      
      BlitzVec nabla_dot_tau(3);
      nabla_dot_tau(0) = blitz::sum(enr_derxy_stress(0,_)*etau(0,_)) + blitz::sum(enr_derxy_stress(1,_)*etau(3,_)) + blitz::sum(enr_derxy_stress(2,_)*etau(4,_));
      nabla_dot_tau(1) = blitz::sum(enr_derxy_stress(3,_)*etau(0,_)) + blitz::sum(enr_derxy_stress(1,_)*etau(1,_)) + blitz::sum(enr_derxy_stress(2,_)*etau(5,_));
      nabla_dot_tau(2) = blitz::sum(enr_derxy_stress(0,_)*etau(4,_)) + blitz::sum(enr_derxy_stress(1,_)*etau(5,_)) + blitz::sum(enr_derxy_stress(2,_)*etau(2,_));
      
//      const BlitzMat epstau(tau*reciproke_visc);
//      
//      if (ele->Id() == 0)
//      {
//          cout << endl;
//          cout << "eps^tau: " << epstau << endl;
//          cout << "eps^u:   " << eps << endl;
//          cout << "pressure = " << press << endl;
//      }


      // get bodyforce in gausspoint
      BlitzVec bodyforce_(3);
      bodyforce_ = 0.0;
  //    if ( not (XFEM::PositionWithinCondition(cellcenterpos, 1, ih->cutterdis())))
  //    {
  //        bodyforce_(1) = 1.0;
  //    }
      //////////////////////////////////////////BlitzVec bodyforce_(blitz::sum(enr_edeadng_(i,j)*enr_funct_(j),j));

      // perform integration for entire matrix and rhs

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

      // calculate tau
      // stabilization parameters for stationary case

      // compute tau_Mu
      BlitzVec tau_stab(3);
      const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
      const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
      tau_stab(0) = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);

      // compute tau_Mp
      const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
      const double xi_tau_mp = DMAX(re_tau_mp,1.0);
      tau_stab(1) = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);

      // compute tau_C
      const double xi_tau_c = DMIN(re_tau_mp, 1.0);
      tau_stab(2) = 0.5*vel_norm*hk*xi_tau_c;
      
      // stabilisation parameter
      const double tau_M  = tau_stab(0)*fac;
      const double tau_Mp = tau_stab(1)*fac;
      const double tau_C  = tau_stab(2)*fac;

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

      /*------------------------- evaluate rhs vector at integration point ---*/
      //   rhsint_ = histvec_(i) + bodyforce_(i);
      // histvec is always zero in stationary case (!):
      const BlitzVec rhsint_(bodyforce_(i));

      /*----------------- get numerical representation of single operators ---*/

      /* Convective term  u_old * grad u_old: */
      const BlitzVec conv_old_(blitz::sum(vderxy_(i, j)*velint_(j), j));
      
      BlitzVec res_old(3);
      if (numparamtauxx == 0)
      {
          /* Viscous term  div epsilon(u_old) */
          BlitzVec visc_old_(3);
          visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,3) + vderxy2_(0,2) + vderxy2_(2,4));
          visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,3) + vderxy2_(1,2) + vderxy2_(2,5));
          visc_old_(2) = vderxy2_(2,2) + 0.5 * (vderxy2_(2,0) + vderxy2_(0,4) + vderxy2_(2,1) + vderxy2_(1,5));
          
          // evaluate residual once for all stabilisation right hand sides
          res_old = -rhsint_+(conv_old_+gradp_-2.0*visc*visc_old_);
      }
      else
      {
          /* Viscous term  div epsilon(u_old) */
          res_old = -rhsint_+(conv_old_+gradp_-nabla_dot_tau);
      }
      
      /* Reactive term  u:  funct */
      /* linearise convective term */

      /*--- convective part u_old * grad (funct) --------------------------*/
      /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
         with  N .. form function matrix                                   */
      const BlitzVec enr_conv_c_(blitz::sum(enr_derxy_(j,i)*velint_(j), j));

  //    /*--- convective grid part u_G * grad (funct) -----------------------*/
  //    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
  //    enr_conv_g_ = 0.0;


      /*--- reactive part funct * grad (u_old) ----------------------------*/
      /* /                                     \
         |  u_old_x,x   u_old_x,y   u_old x,z  |
         |                                     |
         |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
         |                                     |
         |  u_old_z,x   u_old_z,y   u_old_z,z  |
         \                                     /
         with  N .. form function matrix                                   */
      const blitz::Array<double,3> enr_conv_r_(vderxy_(i, j)*shp(k));

      /*--- viscous term  - grad * epsilon(u): ----------------------------*/
      /*   /                                                \
           |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
         1 |                                                |
       - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
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

      /* pressure gradient term derxy, funct without or with integration   *
       * by parts, respectively                                            */

      /*--------------------------------- now build single stiffness terms ---*/
      {

        //----------------------------------------------------------------------
        //                            GALERKIN PART

        /* convection, convective part */
        /*
                     /                       \
                    |  / n+1       \          |
                    | | u   o nabla | Du , v  |
                    |  \ (i)       /          |
                     \                       /
        */
        assembler.template Matrix<Velx,Velx>(shp, fac, enr_conv_c_);
        assembler.template Matrix<Vely,Vely>(shp, fac, enr_conv_c_);
        assembler.template Matrix<Velz,Velz>(shp, fac, enr_conv_c_);

        /* convection */
        assembler.template Vector<Velx>(enr_conv_r_(0, 0, _), -fac*velint_(0));
        assembler.template Vector<Velx>(enr_conv_r_(0, 1, _), -fac*velint_(1));
        assembler.template Vector<Velx>(enr_conv_r_(0, 2, _), -fac*velint_(2));

        assembler.template Vector<Vely>(enr_conv_r_(1, 0, _), -fac*velint_(0));
        assembler.template Vector<Vely>(enr_conv_r_(1, 1, _), -fac*velint_(1));
        assembler.template Vector<Vely>(enr_conv_r_(1, 2, _), -fac*velint_(2));

        assembler.template Vector<Velz>(enr_conv_r_(2, 0, _), -fac*velint_(0));
        assembler.template Vector<Velz>(enr_conv_r_(2, 1, _), -fac*velint_(1));
        assembler.template Vector<Velz>(enr_conv_r_(2, 2, _), -fac*velint_(2));

        /* Viskositaetsterm */
        if (numparamtauxx == 0)
        {
            /*
                          /                        \
                         |       /  \         / \   |
                         |  tau | Du | , eps | v |  |
                         |       \  /         \ /   |
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
        }
        else
        {
            /*
                          /                  \
                         |              / \   |
                         |  Dtau , eps | v |  |
                         |              \ /   |
                          \                  /
            */
            assembler.template Matrix<Velx,Tauxx>(shp_dx, 2.0*fac*0.5, shp_tau);
            assembler.template Matrix<Velx,Tauxy>(shp_dy,     fac*0.5, shp_tau);
            assembler.template Matrix<Velx,Tauxz>(shp_dz,     fac*0.5, shp_tau);
            assembler.template Matrix<Velx,Tauyx>(shp_dy,     fac*0.5, shp_tau);
            assembler.template Matrix<Velx,Tauzx>(shp_dz,     fac*0.5, shp_tau);

            assembler.template Matrix<Vely,Tauxy>(shp_dx,     fac*0.5, shp_tau);
            assembler.template Matrix<Vely,Tauyx>(shp_dx,     fac*0.5, shp_tau);
            assembler.template Matrix<Vely,Tauyy>(shp_dy, 2.0*fac*0.5, shp_tau);
            assembler.template Matrix<Vely,Tauyz>(shp_dz,     fac*0.5, shp_tau);
            assembler.template Matrix<Vely,Tauzy>(shp_dz,     fac*0.5, shp_tau);

            assembler.template Matrix<Velz,Tauxz>(shp_dx,     fac*0.5, shp_tau);
            assembler.template Matrix<Velz,Tauyz>(shp_dy,     fac*0.5, shp_tau);
            assembler.template Matrix<Velz,Tauzx>(shp_dx,     fac*0.5, shp_tau);
            assembler.template Matrix<Velz,Tauzy>(shp_dy,     fac*0.5, shp_tau);
            assembler.template Matrix<Velz,Tauzz>(shp_dz, 2.0*fac*0.5, shp_tau);
            
            assembler.template Vector<Velx>(shp_dx, -2.0*fac*0.5*tau(0, 0));
            assembler.template Vector<Velx>(shp_dy,     -fac*0.5*tau(0, 1));
            assembler.template Vector<Velx>(shp_dy,     -fac*0.5*tau(1, 0));
            assembler.template Vector<Velx>(shp_dz,     -fac*0.5*tau(0, 2));
            assembler.template Vector<Velx>(shp_dz,     -fac*0.5*tau(2, 0));

            assembler.template Vector<Vely>(shp_dx,     -fac*0.5*tau(1, 0));
            assembler.template Vector<Vely>(shp_dx,     -fac*0.5*tau(0, 1));
            assembler.template Vector<Vely>(shp_dy, -2.0*fac*0.5*tau(1, 1));
            assembler.template Vector<Vely>(shp_dz,     -fac*0.5*tau(1, 2));
            assembler.template Vector<Vely>(shp_dz,     -fac*0.5*tau(2, 1));

            assembler.template Vector<Velz>(shp_dx,     -fac*0.5*tau(2, 0));
            assembler.template Vector<Velz>(shp_dx,     -fac*0.5*tau(0, 2));
            assembler.template Vector<Velz>(shp_dy,     -fac*0.5*tau(2, 1));
            assembler.template Vector<Velz>(shp_dy,     -fac*0.5*tau(1, 2));
            assembler.template Vector<Velz>(shp_dz, -2.0*fac*0.5*tau(2, 2));
        }


        /* Druckterm */
        /*
                        /                \
                       |                  |
                       |  Dp , nabla o v  |
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
                       /                \
                      |                  |
                      | nabla o Du  , q  |
                      |                  |
                       \                /
        */
        assembler.template Matrix<Pres,Velx>(shp, fac, shp_dx);
        assembler.template Matrix<Pres,Vely>(shp, fac, shp_dy);
        assembler.template Matrix<Pres,Velz>(shp, fac, shp_dz);
        
        assembler.template Vector<Pres>(enr_conv_r_(0, 0, _), -fac);
        assembler.template Vector<Pres>(enr_conv_r_(1, 1, _), -fac);
        assembler.template Vector<Pres>(enr_conv_r_(2, 2, _), -fac);
        
//        if (numparamtauxx != 0)
//        {
//            /* Divergenzfreiheit - continuity equation*/
//            /*
//                           /                       
//                          |                         |
//                          | trace (eps(Dtau))  , q  |
//                          |                         |
//                           \                       /
//            */
//            assembler.template Matrix<Pres,Tauxx>(shp, 0.5*reciproke_visc*fac, shp_tau);
//            assembler.template Matrix<Pres,Tauyy>(shp, 0.5*reciproke_visc*fac, shp_tau);
//            assembler.template Matrix<Pres,Tauzz>(shp, 0.5*reciproke_visc*fac, shp_tau);
//            
//            assembler.template Vector<Pres>(shp, -0.5*reciproke_visc*fac*tau(0,0));
//            assembler.template Vector<Pres>(shp, -0.5*reciproke_visc*fac*tau(1,1));
//            assembler.template Vector<Pres>(shp, -0.5*reciproke_visc*fac*tau(2,2));
//        }
        

        

        // source term of the right hand side
        assembler.template Vector<Velx>(shp, fac*rhsint_(0));
        assembler.template Vector<Vely>(shp, fac*rhsint_(1));
        assembler.template Vector<Velz>(shp, fac*rhsint_(2));


        // Hellinger-Reissner terms
        if (numparamtauxx != 0)
        {
            /* quadratic part */
            /*
                           /                     \
                          |                       |
                          | eps(Dtau)  , virt tau |
                          |                       |
                           \                     /
            */
            //cout << "HR: quadratic part" << endl;
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
                           /                   \
                          |                     |
                          | eps(Du)  , virt tau |
                          |                     |
                           \                   /
            */
            //cout << "HR: rectangular part" << endl;
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
            assembler.template Vector<Tauxy>(shp_tau,     -fac*0.5*vderxy_(0, 1));
            assembler.template Vector<Tauxy>(shp_tau,     -fac*0.5*vderxy_(1, 0));
            assembler.template Vector<Tauxz>(shp_tau,     -fac*0.5*vderxy_(0, 2));
            assembler.template Vector<Tauxz>(shp_tau,     -fac*0.5*vderxy_(2, 0));

            assembler.template Vector<Tauyx>(shp_tau,     -fac*0.5*vderxy_(1, 0));
            assembler.template Vector<Tauyx>(shp_tau,     -fac*0.5*vderxy_(0, 1));
            assembler.template Vector<Tauyy>(shp_tau, -2.0*fac*0.5*vderxy_(1, 1));
            assembler.template Vector<Tauyz>(shp_tau,     -fac*0.5*vderxy_(1, 2));
            assembler.template Vector<Tauyz>(shp_tau,     -fac*0.5*vderxy_(2, 1));

            assembler.template Vector<Tauzx>(shp_tau,     -fac*0.5*vderxy_(2, 0));
            assembler.template Vector<Tauzx>(shp_tau,     -fac*0.5*vderxy_(0, 2));
            assembler.template Vector<Tauzy>(shp_tau,     -fac*0.5*vderxy_(2, 1));
            assembler.template Vector<Tauzy>(shp_tau,     -fac*0.5*vderxy_(1, 2));
            assembler.template Vector<Tauzz>(shp_tau, -2.0*fac*0.5*vderxy_(2, 2));
        }

        //----------------------------------------------------------------------
        //                 PRESSURE STABILISATION PART

        if(pstab)
        {


          /* pressure stabilisation: convection, convective part */
          /*
                    /                            \
                   |  / n+1       \               |
                   | | u   o nabla | Du , nabla q |
                   |  \ (i)       /               |
                    \                            /
          */
          assembler.template Matrix<Pres,Velx>(shp_dx, tau_Mp, enr_conv_c_);
          assembler.template Matrix<Pres,Vely>(shp_dy, tau_Mp, enr_conv_c_);
          assembler.template Matrix<Pres,Velz>(shp_dz, tau_Mp, enr_conv_c_);

          /* pressure stabilisation: viscosity (-L_visc_u) */
          if (numparamtauxx == 0)
          {
              /*
                       /                              \
                      |               /  \             |
                      |  nabla o tau | Du | , nabla q  |
                      |               \  /             |
                       \                              /
              */
              assembler.template Matrix<Pres,Velx>(shp_dx, -2.0*visc*tau_Mp, enr_viscs2_(0, 0, _));
              assembler.template Matrix<Pres,Velx>(shp_dy, -2.0*visc*tau_Mp, enr_viscs2_(0, 1, _));
              assembler.template Matrix<Pres,Velx>(shp_dz, -2.0*visc*tau_Mp, enr_viscs2_(0, 2, _));
    
              assembler.template Matrix<Pres,Vely>(shp_dx, -2.0*visc*tau_Mp, enr_viscs2_(0, 1, _));
              assembler.template Matrix<Pres,Vely>(shp_dy, -2.0*visc*tau_Mp, enr_viscs2_(1, 1, _));
              assembler.template Matrix<Pres,Vely>(shp_dz, -2.0*visc*tau_Mp, enr_viscs2_(1, 2, _));
    
              assembler.template Matrix<Pres,Velz>(shp_dx, -2.0*visc*tau_Mp, enr_viscs2_(0, 2, _));
              assembler.template Matrix<Pres,Velz>(shp_dy, -2.0*visc*tau_Mp, enr_viscs2_(1, 2, _));
              assembler.template Matrix<Pres,Velz>(shp_dz, -2.0*visc*tau_Mp, enr_viscs2_(2, 2, _));
          }
          else
          {
              /*
                       /                        \
                      |                          |
                      |  nabla o Dtau , nabla q  |
                      |                          |
                       \                        /
              */              
              assembler.template Matrix<Pres,Tauxx>(shp_dx, -tau_Mp, shp_tau_dx);
              assembler.template Matrix<Pres,Tauyx>(shp_dy, -tau_Mp, shp_tau_dx);
              assembler.template Matrix<Pres,Tauzx>(shp_dz, -tau_Mp, shp_tau_dx);
    
              assembler.template Matrix<Pres,Tauyx>(shp_dx, -tau_Mp, shp_tau_dy);
              assembler.template Matrix<Pres,Tauyy>(shp_dy, -tau_Mp, shp_tau_dy);
              assembler.template Matrix<Pres,Tauyz>(shp_dz, -tau_Mp, shp_tau_dy);
    
              assembler.template Matrix<Pres,Tauzx>(shp_dx, -tau_Mp, shp_tau_dz);
              assembler.template Matrix<Pres,Tauzy>(shp_dy, -tau_Mp, shp_tau_dz);
              assembler.template Matrix<Pres,Tauzz>(shp_dz, -tau_Mp, shp_tau_dz);
          }
                  
          /* pressure stabilisation: pressure( L_pres_p) */
          /*
                      /                    \
                     |                      |
                     |  nabla Dp , nabla q  |
                     |                      |
                      \                    /
          */
          assembler.template Matrix<Pres,Pres>(shp_dx, tau_Mp, shp_dx);
          assembler.template Matrix<Pres,Pres>(shp_dy, tau_Mp, shp_dy);
          assembler.template Matrix<Pres,Pres>(shp_dz, tau_Mp, shp_dz);

          // pressure stabilisation
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
                 /                                           \
                |    / n+1        \        / n+1        \     |
                |   | u    o nabla | Du , | u    o nabla | v  |
                |    \ (i)        /        \ (i)        /     |
                 \                                           /
            */
            assembler.template Matrix<Velx,Velx>(enr_conv_c_, tau_M, enr_conv_c_);
            assembler.template Matrix<Vely,Vely>(enr_conv_c_, tau_M, enr_conv_c_);
            assembler.template Matrix<Velz,Velz>(enr_conv_c_, tau_M, enr_conv_c_);
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
                      /                              \
                     |              / n+1       \     |
                     |  nabla Dp , | u   o nabla | v  |
                     |              \ (i)       /     |
                      \                              /
            */
            assembler.template Matrix<Velx,Pres>(enr_conv_c_, tau_M, shp_dx);
            assembler.template Matrix<Vely,Pres>(enr_conv_c_, tau_M, shp_dy);
            assembler.template Matrix<Velz,Pres>(enr_conv_c_, tau_M, shp_dz);

            /* supg stabilisation: viscous part  (-L_visc_u) */
            if (numparamtauxx == 0)
            {
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
            }
            else
            {
                /*
                      /                                  \
                     |                 / n+1        \     |
                     |  nabla o Dtau, | u    o nabla | v  |
                     |                 \ (i)        /     |
                      \                                  /
                */
                assembler.template Matrix<Velx,Tauxx>(enr_conv_c_, -tau_M, shp_tau_dx);
                assembler.template Matrix<Velx,Tauxy>(enr_conv_c_, -tau_M, shp_tau_dy);
                assembler.template Matrix<Velx,Tauxz>(enr_conv_c_, -tau_M, shp_tau_dz);
    
                assembler.template Matrix<Vely,Tauxx>(enr_conv_c_, -tau_M, shp_tau_dx);
                assembler.template Matrix<Vely,Tauxy>(enr_conv_c_, -tau_M, shp_tau_dy);
                assembler.template Matrix<Vely,Tauxz>(enr_conv_c_, -tau_M, shp_tau_dz);
    
                assembler.template Matrix<Velz,Tauxx>(enr_conv_c_, -tau_M, shp_tau_dx);
                assembler.template Matrix<Velz,Tauxy>(enr_conv_c_, -tau_M, shp_tau_dy);
                assembler.template Matrix<Velz,Tauxz>(enr_conv_c_, -tau_M, shp_tau_dz);
            }

            if (newton)
            {
                dserror("not translated to XFEM assembler yet!");
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
          
            const double tau_C_divunp=tau_C*(vderxy_(0, 0)+vderxy_(1, 1)+vderxy_(2, 2));

            assembler.template Vector<Velx>(shp_dx, -tau_C_divunp);
            assembler.template Vector<Vely>(shp_dy, -tau_C_divunp);
            assembler.template Vector<Velz>(shp_dz, -tau_C_divunp);

        } // endif cstab
      } // end block stiffness terms
    } // end loop over gauss points
    } // end loop over integration cells
    return;
  }
  
  void callSysmat(
          XFEM::AssemblyType assembly_type,
          DRT::Elements::XFluid3* ele,
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
          {
              Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
                      ele, ih, eleDofManager, locval, estif, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex20:
          {
              Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
                      ele, ih, eleDofManager, locval, estif, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex27:
          {
              Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
                      ele, ih, eleDofManager, locval, estif, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          default:
              dserror("not templated yet");
          };
      }
      else
      {
          switch (ele->Shape())
          {
          case DRT::Element::hex8:
          {
              Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
                      ele, ih, eleDofManager, locval, estif, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex20:
          {
              Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
                      ele, ih, eleDofManager, locval, estif, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex27:
          {
              Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
                      ele, ih, eleDofManager, locval, estif, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          default:
              dserror("not templated yet");
          };
      }
  }


#endif

#endif
#endif
