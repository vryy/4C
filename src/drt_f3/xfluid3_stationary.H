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
  std::vector<DRT::Condition*> myneumcond;
  const int iel_ = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
  BlitzMat edeadng_(3,iel_,blitz::ColumnMajorArray<2>());
  dserror("not adapted to xfem (, yet)!!!");
  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<iel_;inode++)
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

  if (nodecount == iel_)
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
    for (int jnode=0; jnode<iel_; jnode++)
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
    const int iel_ = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
    Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,iel_);

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
          int iel_)
  {
      return dofman.NumDofPerField(field);
  }
  template<> int getNumParam<XFEM::standard_assembly>(
          const XFEM::ElementDofManager& dofman,
          const XFEM::PHYSICS::Field     field,
          int iel_)
  {
      return iel_;
  }
  
  /// calculating element matrix and rhs
  /*!
    Calculate matrix and rhs for stationary problem formulation
    Here the whole thing is hidden.

    \param ele            (i) the element those matrix is calculated
    \param evelnp         (i) nodal velocities at n+1
    \param eprenp         (i) nodal pressure
    \param estif          (o) element matrix to calculate
    \param esv            (o) fine-scale subgrid-viscosity element matrix
    \param eforce         (o) element rhs to calculate
    \param material       (i) fluid material
    \param pseudotime     (i) current simulation pseudotime
    \param newton         (i) full Newton or fixed-point-like
    \param pstab          (i) boolean flag for stabilisation 
    \param supg           (i) boolean flag for stabilisation
    \param vstab          (i) boolean flag for stabilisation
    \param cstab          (i) boolean flag for stabilisation
   */
  template <DRT::Element::DiscretizationType DISTYPE,
            XFEM::AssemblyType ASSTYPE>
  void Sysmat(
          DRT::Elements::XFluid3* ele,
          const RCP<XFEM::InterfaceHandle>  ih,
          const XFEM::ElementDofManager&    dofman,
          const blitz::Array<double,2>&     evelnp,
          const blitz::Array<double,1>&     eprenp,
          blitz::Array<double,2>&           estif,
          blitz::Array<double,2>&           esv,
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
    //! number of nodes for element
    const int iel_ = DRT::Utils::getNumberOfElementNodes<DISTYPE>();
      
    //! information about domain integration cells
    const XFEM::DomainIntCells   domainIntCells   = ih->GetDomainIntCells(ele->Id(),DISTYPE);
    //! information about boundary integration cells
    const XFEM::BoundaryIntCells boundaryIntCells = ih->GetBoundaryIntCells(ele->Id());

    // get node coordinates
    DRT::Node** const nodes = ele->Nodes();
    BlitzMat xyze_(3,iel_,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<iel_; inode++)
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

    // stabilization parameter
    // This has to be done before anything else is calculated because
    // we use the same arrays internally.
    const double hk = DRT::Elements::FLUID::HK<DISTYPE>(evelnp,visc,xyze_);
    const double mk = DRT::Elements::FLUID::MK<DISTYPE>();

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
    
    LocalAssembler<DISTYPE, ASSTYPE> localAssembler;
    
    const blitz::Range _  = blitz::Range::all();
    
    // We define the variables i,j,k to be indices to blitz arrays.
    // These are used for array expressions, that is matrix-vector
    // products in the following.

    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, iel_);

    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamtauxx = 0;//getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, iel_);

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
      BlitzVec funct_(iel_);
      BlitzMat deriv_(3,iel_,blitz::ColumnMajorArray<2>());
      DRT::Utils::shape_function_3D(funct_,e1,e2,e3,DISTYPE);
      DRT::Utils::shape_function_3D_deriv1(deriv_,e1,e2,e3,DISTYPE);

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
      /// transposed of the jacobian matrix d x / d \xi
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

      // compute second global derivative
      BlitzMat derxy2_(6,iel_,blitz::ColumnMajorArray<2>());
      if (higher_order_ele)
      {
        BlitzMat deriv2_(6,iel_,blitz::ColumnMajorArray<2>());
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
      if (ASSTYPE == XFEM::xfem_assembly)
      {
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
      }
      else
      {
          enr_funct_  = funct_;
          enr_derxy_  = derxy_;
          enr_derxy2_ = derxy2_;
      }
      
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

      // get pressurrrrre
      const double press = blitz::sum(enr_funct_*eprenp);

      // get viscous stress unknowns per interface
      vector<BlitzMat> taus;
      const int numinterface = 1;
      for (int iiface = 0; iiface < numinterface; ++iiface)
      {
          BlitzMat tau(3,3,blitz::ColumnMajorArray<2>());

          // fill tau

          // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          taus.push_back(tau);

      }

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

      const BlitzVec shp(enr_funct_);
      const BlitzVec shp_dx(enr_derxy_(0,_));
      const BlitzVec shp_dy(enr_derxy_(1,_));
      const BlitzVec shp_dz(enr_derxy_(2,_));
//      BlitzVec shp_dx(enr_derxy2_(0,_));
//      BlitzVec shp_dy(enr_derxy2_(1,_));
//      BlitzVec shp_dz(enr_derxy2_(2,_));
      
      const BlitzVec shp_tau;
      const BlitzVec shp_tau_dx;
      const BlitzVec shp_tau_dy;
      const BlitzVec shp_tau_dz;

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
      const BlitzVec res_old(-rhsint_+(conv_old_+gradp_-2.0*visc*visc_old_));

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
      enr_viscs2_(0,0,_) = - 0.5 * (2.0 * enr_derxy2_(0,_) + enr_derxy2_(1,_) + enr_derxy2_(2,_));
      enr_viscs2_(0,1,_) = - 0.5 *  enr_derxy2_(3,_);
      enr_viscs2_(0,2,_) = - 0.5 *  enr_derxy2_(4,_);
      enr_viscs2_(1,0,_) = - 0.5 *  enr_derxy2_(3,_);
      enr_viscs2_(1,1,_) = - 0.5 * (enr_derxy2_(0,_) + 2.0 * enr_derxy2_(1,_) + enr_derxy2_(2,_));
      enr_viscs2_(1,2,_) = - 0.5 *  enr_derxy2_(5,_);
      enr_viscs2_(2,0,_) = - 0.5 *  enr_derxy2_(4,_);
      enr_viscs2_(2,1,_) = - 0.5 *  enr_derxy2_(5,_);
      enr_viscs2_(2,2,_) = - 0.5 * (enr_derxy2_(0,_) + enr_derxy2_(1,_) + 2.0 * enr_derxy2_(2,_));

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
        localAssembler.template Matrix<Velx,Velx>(dofman, estif, shp, fac, enr_conv_c_);
        localAssembler.template Matrix<Vely,Vely>(dofman, estif, shp, fac, enr_conv_c_);
        localAssembler.template Matrix<Velz,Velz>(dofman, estif, shp, fac, enr_conv_c_);

        /* convection */
        localAssembler.template Vector<Velx>(dofman, eforce, enr_conv_r_(0, 0, _), -fac*velint_(0));
        localAssembler.template Vector<Velx>(dofman, eforce, enr_conv_r_(0, 1, _), -fac*velint_(1));
        localAssembler.template Vector<Velx>(dofman, eforce, enr_conv_r_(0, 2, _), -fac*velint_(2));

        localAssembler.template Vector<Vely>(dofman, eforce, enr_conv_r_(1, 0, _), -fac*velint_(0));
        localAssembler.template Vector<Vely>(dofman, eforce, enr_conv_r_(1, 1, _), -fac*velint_(1));
        localAssembler.template Vector<Vely>(dofman, eforce, enr_conv_r_(1, 2, _), -fac*velint_(2));

        localAssembler.template Vector<Velz>(dofman, eforce, enr_conv_r_(2, 0, _), -fac*velint_(0));
        localAssembler.template Vector<Velz>(dofman, eforce, enr_conv_r_(2, 1, _), -fac*velint_(1));
        localAssembler.template Vector<Velz>(dofman, eforce, enr_conv_r_(2, 2, _), -fac*velint_(2));

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
            localAssembler.template Matrix<Velx,Velx>(dofman, estif, shp_dx, 2.0*visc*fac, shp_dx);
            localAssembler.template Matrix<Velx,Velx>(dofman, estif, shp_dy,     visc*fac, shp_dy);
            localAssembler.template Matrix<Velx,Velx>(dofman, estif, shp_dz,     visc*fac, shp_dz);
            localAssembler.template Matrix<Velx,Vely>(dofman, estif, shp_dy,     visc*fac, shp_dx);
            localAssembler.template Matrix<Velx,Velz>(dofman, estif, shp_dz,     visc*fac, shp_dx);

            localAssembler.template Matrix<Vely,Velx>(dofman, estif, shp_dx,     visc*fac, shp_dy);
            localAssembler.template Matrix<Vely,Vely>(dofman, estif, shp_dx,     visc*fac, shp_dx);
            localAssembler.template Matrix<Vely,Vely>(dofman, estif, shp_dy, 2.0*visc*fac, shp_dy);
            localAssembler.template Matrix<Vely,Vely>(dofman, estif, shp_dz,     visc*fac, shp_dz);
            localAssembler.template Matrix<Vely,Velz>(dofman, estif, shp_dz,     visc*fac, shp_dy);

            localAssembler.template Matrix<Velz,Velx>(dofman, estif, shp_dx,     visc*fac, shp_dz);
            localAssembler.template Matrix<Velz,Vely>(dofman, estif, shp_dy,     visc*fac, shp_dz);
            localAssembler.template Matrix<Velz,Velz>(dofman, estif, shp_dx,     visc*fac, shp_dx);
            localAssembler.template Matrix<Velz,Velz>(dofman, estif, shp_dy,     visc*fac, shp_dy);
            localAssembler.template Matrix<Velz,Velz>(dofman, estif, shp_dz, 2.0*visc*fac, shp_dz);
            
            localAssembler.template Vector<Velx>(dofman, eforce, shp_dx, -2.0*visc*fac*vderxy_(0, 0));
            localAssembler.template Vector<Velx>(dofman, eforce, shp_dy,     -visc*fac*vderxy_(0, 1));
            localAssembler.template Vector<Velx>(dofman, eforce, shp_dy,     -visc*fac*vderxy_(1, 0));
            localAssembler.template Vector<Velx>(dofman, eforce, shp_dz,     -visc*fac*vderxy_(0, 2));
            localAssembler.template Vector<Velx>(dofman, eforce, shp_dz,     -visc*fac*vderxy_(2, 0));

            localAssembler.template Vector<Vely>(dofman, eforce, shp_dx,     -visc*fac*vderxy_(1, 0));
            localAssembler.template Vector<Vely>(dofman, eforce, shp_dx,     -visc*fac*vderxy_(0, 1));
            localAssembler.template Vector<Vely>(dofman, eforce, shp_dy, -2.0*visc*fac*vderxy_(1, 1));
            localAssembler.template Vector<Vely>(dofman, eforce, shp_dz,     -visc*fac*vderxy_(1, 2));
            localAssembler.template Vector<Vely>(dofman, eforce, shp_dz,     -visc*fac*vderxy_(2, 1));

            localAssembler.template Vector<Velz>(dofman, eforce, shp_dx,     -visc*fac*vderxy_(2, 0));
            localAssembler.template Vector<Velz>(dofman, eforce, shp_dx,     -visc*fac*vderxy_(0, 2));
            localAssembler.template Vector<Velz>(dofman, eforce, shp_dy,     -visc*fac*vderxy_(2, 1));
            localAssembler.template Vector<Velz>(dofman, eforce, shp_dy,     -visc*fac*vderxy_(1, 2));
            localAssembler.template Vector<Velz>(dofman, eforce, shp_dz, -2.0*visc*fac*vderxy_(2, 2));
        }
        else
        {
            for (int iiface = 0; iiface < numinterface; ++iiface)
            {
                /*
                              /                  \
                             |              / \   |
                             |  Dtau , eps | v |  |
                             |              \ /   |
                              \                  /
                */
                localAssembler.template Matrix<Velx,Tauxx>(dofman, estif, shp_dx, 2.0*fac, shp_tau);
                localAssembler.template Matrix<Velx,Tauxy>(dofman, estif, shp_dy,     fac, shp_tau);
                localAssembler.template Matrix<Velx,Tauxz>(dofman, estif, shp_dz,     fac, shp_tau);
                localAssembler.template Matrix<Velx,Tauyx>(dofman, estif, shp_dy,     fac, shp_tau);
                localAssembler.template Matrix<Velx,Tauzx>(dofman, estif, shp_dz,     fac, shp_tau);
    
                localAssembler.template Matrix<Vely,Tauxy>(dofman, estif, shp_dx,     fac, shp_tau);
                localAssembler.template Matrix<Vely,Tauyx>(dofman, estif, shp_dx,     fac, shp_tau);
                localAssembler.template Matrix<Vely,Tauyy>(dofman, estif, shp_dy, 2.0*fac, shp_tau);
                localAssembler.template Matrix<Vely,Tauyz>(dofman, estif, shp_dz,     fac, shp_tau);
                localAssembler.template Matrix<Vely,Tauzy>(dofman, estif, shp_dz,     fac, shp_tau);
    
                localAssembler.template Matrix<Velz,Tauxz>(dofman, estif, shp_dx,     fac, shp_tau);
                localAssembler.template Matrix<Velz,Tauyz>(dofman, estif, shp_dy,     fac, shp_tau);
                localAssembler.template Matrix<Velz,Tauzx>(dofman, estif, shp_dx,     fac, shp_tau);
                localAssembler.template Matrix<Velz,Tauzy>(dofman, estif, shp_dy,     fac, shp_tau);
                localAssembler.template Matrix<Velz,Tauzz>(dofman, estif, shp_dz, 2.0*fac, shp_tau);
                
                BlitzMat tau = taus[iiface];
    
                localAssembler.template Vector<Velx>(dofman, eforce, shp_dx, -2.0*fac*tau(0, 0));
                localAssembler.template Vector<Velx>(dofman, eforce, shp_dy,     -fac*tau(0, 1));
                localAssembler.template Vector<Velx>(dofman, eforce, shp_dy,     -fac*tau(1, 0));
                localAssembler.template Vector<Velx>(dofman, eforce, shp_dz,     -fac*tau(0, 2));
                localAssembler.template Vector<Velx>(dofman, eforce, shp_dz,     -fac*tau(2, 0));
    
                localAssembler.template Vector<Vely>(dofman, eforce, shp_dx,     -fac*tau(1, 0));
                localAssembler.template Vector<Vely>(dofman, eforce, shp_dx,     -fac*tau(0, 1));
                localAssembler.template Vector<Vely>(dofman, eforce, shp_dy, -2.0*fac*tau(1, 1));
                localAssembler.template Vector<Vely>(dofman, eforce, shp_dz,     -fac*tau(1, 2));
                localAssembler.template Vector<Vely>(dofman, eforce, shp_dz,     -fac*tau(2, 1));
    
                localAssembler.template Vector<Velz>(dofman, eforce, shp_dx,     -fac*tau(2, 0));
                localAssembler.template Vector<Velz>(dofman, eforce, shp_dx,     -fac*tau(0, 2));
                localAssembler.template Vector<Velz>(dofman, eforce, shp_dy,     -fac*tau(2, 1));
                localAssembler.template Vector<Velz>(dofman, eforce, shp_dy,     -fac*tau(1, 2));
                localAssembler.template Vector<Velz>(dofman, eforce, shp_dz, -2.0*fac*tau(2, 2));
            }
        }


        /* Druckterm */
        /*
                        /                \
                       |                  |
                       |  Dp , nabla o v  |
                       |                  |
                        \                /
        */
        localAssembler.template Matrix<Velx,Pres>(dofman, estif, shp_dx, -fac, shp);
        localAssembler.template Matrix<Vely,Pres>(dofman, estif, shp_dy, -fac, shp);
        localAssembler.template Matrix<Velz,Pres>(dofman, estif, shp_dz, -fac, shp);
        
        localAssembler.template Vector<Velx>(dofman, eforce, shp_dx, press*fac);
        localAssembler.template Vector<Vely>(dofman, eforce, shp_dy, press*fac);
        localAssembler.template Vector<Velz>(dofman, eforce, shp_dz, press*fac);


        /* Divergenzfreiheit - continuity equation*/
        /*
                       /                \
                      |                  |
                      | nabla o Du  , q  |
                      |                  |
                       \                /
        */
        localAssembler.template Matrix<Pres,Velx>(dofman, estif, shp, fac, shp_dx);
        localAssembler.template Matrix<Pres,Vely>(dofman, estif, shp, fac, shp_dy);
        localAssembler.template Matrix<Pres,Velz>(dofman, estif, shp, fac, shp_dz);
        
        localAssembler.template Vector<Pres>(dofman, eforce, enr_conv_r_(0, 0, _), -fac);
        localAssembler.template Vector<Pres>(dofman, eforce, enr_conv_r_(1, 1, _), -fac);
        localAssembler.template Vector<Pres>(dofman, eforce, enr_conv_r_(2, 2, _), -fac);

        

        // source term of the right hand side
        localAssembler.template Vector<Velx>(dofman, eforce, shp, fac*rhsint_(0));
        localAssembler.template Vector<Vely>(dofman, eforce, shp, fac*rhsint_(1));
        localAssembler.template Vector<Velz>(dofman, eforce, shp, fac*rhsint_(2));


        // Hellinger-Reissner terms
        if (numparamtauxx > 0)
        {
            for (int iiface = 0; iiface < numinterface; ++iiface)
            {
                BlitzMat tau = taus[iiface];
                /* quadratic part */
                /*
                               /                     \
                              |                       |
                              | eps(Dtau)  , virt tau |
                              |                       |
                               \                     /
                */
                localAssembler.template Matrix<Tauxx,Tauxx>(dofman, estif, shp_tau, -reciproke_visc*fac, shp_tau);
                localAssembler.template Matrix<Tauyy,Tauyy>(dofman, estif, shp_tau, -reciproke_visc*fac, shp_tau);
                localAssembler.template Matrix<Tauzz,Tauzz>(dofman, estif, shp_tau, -reciproke_visc*fac, shp_tau);
                localAssembler.template Matrix<Tauxy,Tauxy>(dofman, estif, shp_tau, -reciproke_visc*fac, shp_tau);
                localAssembler.template Matrix<Tauxz,Tauxz>(dofman, estif, shp_tau, -reciproke_visc*fac, shp_tau);
                localAssembler.template Matrix<Tauyz,Tauyz>(dofman, estif, shp_tau, -reciproke_visc*fac, shp_tau);

                localAssembler.template Vector<Tauxx>(dofman, eforce, shp_tau, reciproke_visc*fac*tau(0,0));
                localAssembler.template Vector<Tauyy>(dofman, eforce, shp_tau, reciproke_visc*fac*tau(1,1));
                localAssembler.template Vector<Tauzz>(dofman, eforce, shp_tau, reciproke_visc*fac*tau(2,2));
                localAssembler.template Vector<Tauxy>(dofman, eforce, shp_tau, reciproke_visc*fac*tau(0,1));
                localAssembler.template Vector<Tauxz>(dofman, eforce, shp_tau, reciproke_visc*fac*tau(0,2));
                localAssembler.template Vector<Tauyz>(dofman, eforce, shp_tau, reciproke_visc*fac*tau(1,2));

                /* rectangular part */
                /*
                               /                   \
                              |                     |
                              | eps(Du)  , virt tau |
                              |                     |
                               \                   /
                */
                localAssembler.template Matrix<Tauxx,Velx>(dofman, estif, shp_tau, 2.0*fac, shp_dx);
                localAssembler.template Matrix<Tauxy,Velx>(dofman, estif, shp_tau,     fac, shp_dy);
                localAssembler.template Matrix<Tauxz,Velx>(dofman, estif, shp_tau,     fac, shp_dz);
                localAssembler.template Matrix<Tauyx,Velx>(dofman, estif, shp_tau,     fac, shp_dy);
                localAssembler.template Matrix<Tauzx,Velx>(dofman, estif, shp_tau,     fac, shp_dz);

                localAssembler.template Matrix<Tauxy,Vely>(dofman, estif, shp_tau,     fac, shp_dx);
                localAssembler.template Matrix<Tauyx,Vely>(dofman, estif, shp_tau,     fac, shp_dx);
                localAssembler.template Matrix<Tauyy,Vely>(dofman, estif, shp_tau, 2.0*fac, shp_dy);
                localAssembler.template Matrix<Tauyz,Vely>(dofman, estif, shp_tau,     fac, shp_dz);
                localAssembler.template Matrix<Tauzy,Vely>(dofman, estif, shp_tau,     fac, shp_dz);

                localAssembler.template Matrix<Tauxz,Velz>(dofman, estif, shp_tau,     fac, shp_dx);
                localAssembler.template Matrix<Tauyz,Velz>(dofman, estif, shp_tau,     fac, shp_dy);
                localAssembler.template Matrix<Tauzx,Velz>(dofman, estif, shp_tau,     fac, shp_dx);
                localAssembler.template Matrix<Tauzy,Velz>(dofman, estif, shp_tau,     fac, shp_dy);
                localAssembler.template Matrix<Tauzz,Velz>(dofman, estif, shp_tau, 2.0*fac, shp_dz);


                localAssembler.template Vector<Tauxx>(dofman, eforce, shp_tau, 2.0*fac*vderxy_(0, 0));
                localAssembler.template Vector<Tauxy>(dofman, eforce, shp_tau,     fac*vderxy_(0, 1));
                localAssembler.template Vector<Tauxz>(dofman, eforce, shp_tau,     fac*vderxy_(0, 2));
                localAssembler.template Vector<Tauyx>(dofman, eforce, shp_tau,     fac*vderxy_(0, 1));
                localAssembler.template Vector<Tauzx>(dofman, eforce, shp_tau,     fac*vderxy_(0, 2));

                localAssembler.template Vector<Tauxy>(dofman, eforce, shp_tau,     fac*vderxy_(1, 0));
                localAssembler.template Vector<Tauyx>(dofman, eforce, shp_tau,     fac*vderxy_(1, 0));
                localAssembler.template Vector<Tauyy>(dofman, eforce, shp_tau, 2.0*fac*vderxy_(1, 1));
                localAssembler.template Vector<Tauyz>(dofman, eforce, shp_tau,     fac*vderxy_(1, 2));
                localAssembler.template Vector<Tauzy>(dofman, eforce, shp_tau,     fac*vderxy_(1, 2));

                localAssembler.template Vector<Tauxz>(dofman, eforce, shp_tau,     fac*vderxy_(2, 0));
                localAssembler.template Vector<Tauyz>(dofman, eforce, shp_tau,     fac*vderxy_(2, 1));
                localAssembler.template Vector<Tauzx>(dofman, eforce, shp_tau,     fac*vderxy_(2, 0));
                localAssembler.template Vector<Tauzy>(dofman, eforce, shp_tau,     fac*vderxy_(2, 1));
                localAssembler.template Vector<Tauzz>(dofman, eforce, shp_tau, 2.0*fac*vderxy_(2, 2));
            }
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
          localAssembler.template Matrix<Pres,Velx>(dofman, estif, shp_dx, tau_Mp, enr_conv_c_);
          localAssembler.template Matrix<Pres,Vely>(dofman, estif, shp_dy, tau_Mp, enr_conv_c_);
          localAssembler.template Matrix<Pres,Velz>(dofman, estif, shp_dz, tau_Mp, enr_conv_c_);

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
              localAssembler.template Matrix<Pres,Velx>(dofman, estif, shp_dx, 2.0*visc*tau_Mp, enr_viscs2_(0, 0, _));
              localAssembler.template Matrix<Pres,Velx>(dofman, estif, shp_dy, 2.0*visc*tau_Mp, enr_viscs2_(0, 1, _));
              localAssembler.template Matrix<Pres,Velx>(dofman, estif, shp_dz, 2.0*visc*tau_Mp, enr_viscs2_(0, 2, _));
    
              localAssembler.template Matrix<Pres,Vely>(dofman, estif, shp_dx, 2.0*visc*tau_Mp, enr_viscs2_(0, 1, _));
              localAssembler.template Matrix<Pres,Vely>(dofman, estif, shp_dy, 2.0*visc*tau_Mp, enr_viscs2_(1, 1, _));
              localAssembler.template Matrix<Pres,Vely>(dofman, estif, shp_dz, 2.0*visc*tau_Mp, enr_viscs2_(1, 2, _));
    
              localAssembler.template Matrix<Pres,Velz>(dofman, estif, shp_dx, 2.0*visc*tau_Mp, enr_viscs2_(0, 2, _));
              localAssembler.template Matrix<Pres,Velz>(dofman, estif, shp_dy, 2.0*visc*tau_Mp, enr_viscs2_(1, 2, _));
              localAssembler.template Matrix<Pres,Velz>(dofman, estif, shp_dz, 2.0*visc*tau_Mp, enr_viscs2_(2, 2, _));
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
              for (int iiface = 0; iiface < numinterface; ++iiface)
              {
                  BlitzMat tau = taus[iiface];
                  localAssembler.template Matrix<Pres,Tauxx>(dofman, estif, shp_dx, tau_Mp, shp_tau_dx);
                  localAssembler.template Matrix<Pres,Tauyx>(dofman, estif, shp_dy, tau_Mp, shp_tau_dx);
                  localAssembler.template Matrix<Pres,Tauzx>(dofman, estif, shp_dz, tau_Mp, shp_tau_dx);
        
                  localAssembler.template Matrix<Pres,Tauxy>(dofman, estif, shp_dx, tau_Mp, shp_tau_dy);
                  localAssembler.template Matrix<Pres,Tauyy>(dofman, estif, shp_dy, tau_Mp, shp_tau_dy);
                  localAssembler.template Matrix<Pres,Tauzy>(dofman, estif, shp_dz, tau_Mp, shp_tau_dy);
        
                  localAssembler.template Matrix<Pres,Tauxz>(dofman, estif, shp_dx, tau_Mp, shp_tau_dz);
                  localAssembler.template Matrix<Pres,Tauyz>(dofman, estif, shp_dy, tau_Mp, shp_tau_dz);
                  localAssembler.template Matrix<Pres,Tauzz>(dofman, estif, shp_dz, tau_Mp, shp_tau_dz);
              }
          }
                  
          /* pressure stabilisation: pressure( L_pres_p) */
          /*
                      /                    \
                     |                      |
                     |  nabla Dp , nabla q  |
                     |                      |
                      \                    /
          */
          localAssembler.template Matrix<Pres,Pres>(dofman, estif, shp_dx, tau_Mp, shp_dx);
          localAssembler.template Matrix<Pres,Pres>(dofman, estif, shp_dy, tau_Mp, shp_dy);
          localAssembler.template Matrix<Pres,Pres>(dofman, estif, shp_dz, tau_Mp, shp_dz);

          // pressure stabilisation
          localAssembler.template Vector<Pres>(dofman, eforce, shp_dx, -tau_Mp*res_old(0));
          localAssembler.template Vector<Pres>(dofman, eforce, shp_dy, -tau_Mp*res_old(1));
          localAssembler.template Vector<Pres>(dofman, eforce, shp_dz, -tau_Mp*res_old(2));
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
            localAssembler.template Matrix<Velx,Velx>(dofman, estif, enr_conv_c_, tau_M, enr_conv_c_);
            localAssembler.template Matrix<Vely,Vely>(dofman, estif, enr_conv_c_, tau_M, enr_conv_c_);
            localAssembler.template Matrix<Velz,Velz>(dofman, estif, enr_conv_c_, tau_M, enr_conv_c_);
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
                      /                              \
                     |              / n+1       \     |
                     |  nabla Dp , | u   o nabla | v  |
                     |              \ (i)       /     |
                      \                              /
            */
            localAssembler.template Matrix<Velx,Pres>(dofman, estif, enr_conv_c_, tau_M, shp_dx);
            localAssembler.template Matrix<Vely,Pres>(dofman, estif, enr_conv_c_, tau_M, shp_dy);
            localAssembler.template Matrix<Velz,Pres>(dofman, estif, enr_conv_c_, tau_M, shp_dz);

            /* supg stabilisation: viscous part  (-L_visc_u) */
            /*
                  /                                        \
                 |               /  \    / n+1        \     |
                 |  nabla o eps | Du |, | u    o nabla | v  |
                 |               \  /    \ (i)        /     |
                  \                                        /
            */
            localAssembler.template Matrix<Velx,Velx>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(0, 0, _));
            localAssembler.template Matrix<Velx,Vely>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(0, 1, _));
            localAssembler.template Matrix<Velx,Velz>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(0, 2, _));

            localAssembler.template Matrix<Vely,Velx>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(0, 1, _));
            localAssembler.template Matrix<Vely,Vely>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(1, 1, _));
            localAssembler.template Matrix<Vely,Velz>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(1, 2, _));

            localAssembler.template Matrix<Velz,Velx>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(0, 2, _));
            localAssembler.template Matrix<Velz,Vely>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(1, 2, _));
            localAssembler.template Matrix<Velz,Velz>(dofman, estif, enr_conv_c_, 2.0*visc*tau_M, enr_viscs2_(2, 2, _));

          if (newton)
          {
              dserror("not translated to XFEM localAssembler yet!");
              } // if newton

          // supg stabilisation
          localAssembler.template Vector<Velx>(dofman, eforce, enr_conv_c_, -tau_M*res_old(0));
          localAssembler.template Vector<Vely>(dofman, eforce, enr_conv_c_, -tau_M*res_old(1));
          localAssembler.template Vector<Velz>(dofman, eforce, enr_conv_c_, -tau_M*res_old(2));
        }


        //----------------------------------------------------------------------
        //                       STABILISATION, VISCOUS PART

        if(vstab)
        {
            //!!!!
            //  viscous part of stabilisation is switched off!
            //  vstab is set to false within fluid3_evaluate.cpp
            //!!!!
            dserror("not translated to XFEM localAssembler yet!");
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
            localAssembler.template Matrix<Velx,Velx>(dofman, estif, shp_dx, tau_C, shp_dx);
            localAssembler.template Matrix<Velx,Vely>(dofman, estif, shp_dx, tau_C, shp_dy);
            localAssembler.template Matrix<Velx,Velz>(dofman, estif, shp_dx, tau_C, shp_dz);
          
            localAssembler.template Matrix<Vely,Velx>(dofman, estif, shp_dy, tau_C, shp_dx);
            localAssembler.template Matrix<Vely,Vely>(dofman, estif, shp_dy, tau_C, shp_dy);
            localAssembler.template Matrix<Vely,Velz>(dofman, estif, shp_dy, tau_C, shp_dz);
          
            localAssembler.template Matrix<Velz,Velx>(dofman, estif, shp_dz, tau_C, shp_dx);
            localAssembler.template Matrix<Velz,Vely>(dofman, estif, shp_dz, tau_C, shp_dy);
            localAssembler.template Matrix<Velz,Velz>(dofman, estif, shp_dz, tau_C, shp_dz);
          
            const double tau_C_divunp=tau_C*(vderxy_(0, 0)+vderxy_(1, 1)+vderxy_(2, 2));

            localAssembler.template Vector<Velx>(dofman, eforce, shp_dx, -tau_C_divunp);
            localAssembler.template Vector<Vely>(dofman, eforce, shp_dy, -tau_C_divunp);
            localAssembler.template Vector<Velz>(dofman, eforce, shp_dz, -tau_C_divunp);

        } // endif cstab
      } // end block stiffness terms
    } // end loop over gauss points
    } // end loop over integration cells
    return;
  }
  
  void callSysmat(DRT::Element::DiscretizationType distype,
          XFEM::AssemblyType assembly_type,
          DRT::Elements::XFluid3* ele,
          const RCP<XFEM::InterfaceHandle>  ih,
          const XFEM::ElementDofManager&    eleDofManager,
          const blitz::Array<double,2>&     evelnp,
          const blitz::Array<double,1>&     eprenp,
          blitz::Array<double,2>&           estif,
          blitz::Array<double,2>&           esv,
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
          switch (distype)
          {
          case DRT::Element::hex8:
          {
              Sysmat<DRT::Element::hex8,XFEM::standard_assembly>(
                      ele, ih, eleDofManager, evelnp, eprenp, estif, esv, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex20:
          {
              Sysmat<DRT::Element::hex20,XFEM::standard_assembly>(
                      ele, ih, eleDofManager, evelnp, eprenp, estif, esv, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex27:
          {
              Sysmat<DRT::Element::hex27,XFEM::standard_assembly>(
                      ele, ih, eleDofManager, evelnp, eprenp, estif, esv, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          default:
              dserror("not templated yet");
          };
      }
      else
      {
          switch (distype)
          {
          case DRT::Element::hex8:
          {
              Sysmat<DRT::Element::hex8,XFEM::xfem_assembly>(
                      ele, ih, eleDofManager, evelnp, eprenp, estif, esv, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex20:
          {
              Sysmat<DRT::Element::hex20,XFEM::xfem_assembly>(
                      ele, ih, eleDofManager, evelnp, eprenp, estif, esv, eforce,
                      material, pseudotime, newton, pstab, supg, vstab, cstab);
              break;
          }
          case DRT::Element::hex27:
          {
              Sysmat<DRT::Element::hex27,XFEM::xfem_assembly>(
                      ele, ih, eleDofManager, evelnp, eprenp, estif, esv, eforce,
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
