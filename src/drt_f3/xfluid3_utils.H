/*!
\file xfluid3_utils.H

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFLUID3_UTILS_H
#define XFLUID3_UTILS_H

#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/dof_management_element.H"
#include "../drt_xfem/xfem_enums.H"
#include <Epetra_SerialDenseSolver.h>


namespace XFEM
{
  
/*! get number of unknowns for this field
 * use the dofmanager in the general case,
 * otherwise use a given number numnode
 */
template<XFEM::AssemblyType>
inline int getNumParam(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        const int)
{
    return dofman.NumDofPerField(field);
}
template<> 
inline int getNumParam<XFEM::standard_assembly>(
        const XFEM::ElementDofManager&,
        const XFEM::PHYSICS::Field,
        const int numnode)
{
    return numnode;
}

//! determine number of parameter
template<int NUMNODE, XFEM::AssemblyType ASSTYPE>
struct NumParam
{
  static int get(
     const XFEM::ElementDofManager& dofman,
     const XFEM::PHYSICS::Field     field
     )
  {
    return dofman.NumDofPerField(field);
  }
};

/// determine number of parameter for standard_assembly
template<int NUMNODE> struct NumParam<NUMNODE,XFEM::standard_assembly>
{
  static int get(
      const XFEM::ElementDofManager& dofman,
      const XFEM::PHYSICS::Field     field
      )
  {
    return NUMNODE;
  }
};


template<class M1, class M2, class M3>
static void fill_tau(const int numparamtau, const M1&  shp_tau, const M2& etau, M3& tau)
{
    tau = 0.0;
    for (int i = 0; i < numparamtau; ++i)
    {
        const double fac = shp_tau(i);
        tau(0,0) += fac*etau(0,i);
        tau(1,1) += fac*etau(1,i);
        tau(2,2) += fac*etau(2,i);
        tau(0,1) += fac*etau(3,i);
        tau(0,2) += fac*etau(4,i);
        tau(1,2) += fac*etau(5,i);
    }
    tau(1,0) = tau(0,1);
    tau(2,0) = tau(0,2);
    tau(2,1) = tau(1,2);
}

/*!
 * \brief calculate the body force from the nodes condition
 * 
 *  get the body force in the nodes of the element (private) gammi 04/07
 *  the Neumann condition associated with the nodes is stored in the
 *  array edeadng only if all nodes have a VolumeNeumann condition
 *
 */

template <DRT::Element::DiscretizationType DISTYPE>
LINALG::SerialDenseMatrix BodyForce(DRT::Node** nodes, const double pseudotime)
{
  dserror("not adapted to xfem (, yet)!!!");
    
  std::vector<DRT::Condition*> myneumcond;
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  // check whether all nodes have a unique VolumeNeumann condition
  int nodecount = 0;
  for (int inode=0;inode<numnode;++inode)
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

  LINALG::SerialDenseMatrix edeadng_(3,numnode);
  if (nodecount == numnode)
  {
    // find out whether we will use a (pseudo-)time curve
    const std::vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
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
    for (int jnode=0; jnode<numnode; ++jnode)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const std::vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const std::vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

      for(int isd=0;isd<3;++isd)
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
template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2, class M3, class M4, class M5>
void gder2(
        const M1& xjm,
        const M2& derxy,
        const M3& deriv2,
        const M4& xyze,
        M5& derxy2)
{
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // initialize and zero out everything
    static LINALG::FixedSizeSerialDenseMatrix<6,6> bm;

    // calculate elements of jacobian_bar matrix
    bm(0,0) = xjm(0,0)*xjm(0,0);
    bm(1,0) = xjm(1,0)*xjm(1,0);
    bm(2,0) = xjm(2,0)*xjm(2,0);
    bm(3,0) = xjm(0,0)*xjm(1,0);
    bm(4,0) = xjm(0,0)*xjm(2,0);
    bm(5,0) = xjm(2,0)*xjm(1,0);

    bm(0,1) = xjm(0,1)*xjm(0,1);
    bm(1,1) = xjm(1,1)*xjm(1,1);
    bm(2,1) = xjm(2,1)*xjm(2,1);
    bm(3,1) = xjm(0,1)*xjm(1,1);
    bm(4,1) = xjm(0,1)*xjm(2,1);
    bm(5,1) = xjm(2,1)*xjm(1,1);

    bm(0,2) = xjm(0,2)*xjm(0,2);
    bm(1,2) = xjm(1,2)*xjm(1,2);
    bm(2,2) = xjm(2,2)*xjm(2,2);
    bm(3,2) = xjm(0,2)*xjm(1,2);
    bm(4,2) = xjm(0,2)*xjm(2,2);
    bm(5,2) = xjm(2,2)*xjm(1,2);

    bm(0,3) = 2.*xjm(0,0)*xjm(0,1);
    bm(1,3) = 2.*xjm(1,0)*xjm(1,1);
    bm(2,3) = 2.*xjm(2,0)*xjm(2,1);
    bm(3,3) = xjm(0,0)*xjm(1,1)+xjm(1,0)*xjm(0,1);
    bm(4,3) = xjm(0,0)*xjm(2,1)+xjm(2,0)*xjm(0,1);
    bm(5,3) = xjm(1,0)*xjm(2,1)+xjm(2,0)*xjm(1,1);

    bm(0,4) = 2.*xjm(0,0)*xjm(0,2);
    bm(1,4) = 2.*xjm(1,0)*xjm(1,2);
    bm(2,4) = 2.*xjm(2,0)*xjm(2,2);
    bm(3,4) = xjm(0,0)*xjm(1,2)+xjm(1,0)*xjm(0,2);
    bm(4,4) = xjm(0,0)*xjm(2,2)+xjm(2,0)*xjm(0,2);
    bm(5,4) = xjm(1,0)*xjm(2,2)+xjm(2,0)*xjm(1,2);

    bm(0,5) = 2.*xjm(0,1)*xjm(0,2);
    bm(1,5) = 2.*xjm(1,1)*xjm(1,2);
    bm(2,5) = 2.*xjm(2,1)*xjm(2,2);
    bm(3,5) = xjm(0,1)*xjm(1,2)+xjm(1,1)*xjm(0,2);
    bm(4,5) = xjm(0,1)*xjm(2,2)+xjm(2,1)*xjm(0,2);
    bm(5,5) = xjm(1,1)*xjm(2,2)+xjm(2,1)*xjm(1,2);

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
    |                                            +-+-+-+
    |
    |        xder2               deriv2          xyze^T
    |
    |
    |                                     +-                  -+
    |                                     | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                                     | dr^2   dr^2   dr^2 |
    |                                     |                    |
    |                                     | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                                     | ds^2   ds^2   ds^2 |
    |                                     |                    |
    |                                     | d^2x   d^2y   d^2z |
    |                                     | ----   ----   ---- |
    |                                     | dt^2   dt^2   dt^2 |
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
    |                                     +-                  -+
    |
    |
    */

    //xder2_ = deriv2_(i,k)*xyze_(j,k);
    static LINALG::FixedSizeSerialDenseMatrix<6,3> xder2;
    xder2.MultiplyNT(deriv2,xyze);

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

    //derxy2_ = - xder2_(i,k)*derxy_(k,j);
    derxy2.Multiply(-1.0,xder2,derxy);

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
    derxy2.Update(1.0,deriv2,1.0);

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
    |                                            |
    |                                            chainrulerhs[i]
    |
    |       yields
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
    |                      +-+-+-+-+
    */

    LINALG::FixedSizeSerialDenseSolver<6,6,numnode> solver;
    solver.SetMatrix(bm);

    // No need for a separate rhs. We assemble the rhs to the solution
    // vector. The solver will destroy the rhs and return the solution.
    solver.SetVectors(derxy2,derxy2);
    solver.Solve();
}


/*!
 \brief get optimal gaussrule for given discretization type
        template version
 \return enum with optimal gausrule for this element
 */
static DRT::UTILS::GaussRule3D getOptimalVelocityGaussrule(DRT::Element::DiscretizationType  distype)
{
    DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    switch (distype)
    {
    case DRT::Element::hex8:
        rule = DRT::UTILS::intrule_hex_8point;
        break;
    case DRT::Element::hex20: case DRT::Element::hex27:
        rule = DRT::UTILS::intrule_hex_27point;
        break;
    case DRT::Element::tet4:
        rule = DRT::UTILS::intrule_tet_4point;
        break;
    case DRT::Element::tet10:
        rule = DRT::UTILS::intrule_tet_11point;
        break;
    default:
        dserror("unknown distype for gaussrule initialization");
  }
  return rule;
}

/*!
 \brief get optimal gaussrule for given discretization type

 \return enum with optimal gausrule for this element
 */
inline DRT::UTILS::GaussRule3D getXFEMGaussrule(
    const bool                              intersected,
    const DRT::Element::DiscretizationType  celldistype,
    const DRT::Element::DiscretizationType  xeledistype
    )
{
  GEO::EleGeoType eleGeoType = GEO::CARTESIAN;
  
  DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
  if (not intersected)
  {
    rule = getOptimalVelocityGaussrule(celldistype);
  }
  else
  {
    if (celldistype == DRT::Element::tet4 or celldistype == DRT::Element::tet10)
    {
      switch (xeledistype)
      {
        case DRT::Element::hex8:
          if (eleGeoType == GEO::CARTESIAN) // constant jacobian
            rule = DRT::UTILS::intrule_tet_24point;
          else
            rule = DRT::UTILS::intrule_tet_45point;
          break;
        case DRT::Element::hex20:
          if (eleGeoType == GEO::CARTESIAN) // constant jacobian
            rule = DRT::UTILS::intrule_tet_45point;
          else
            rule = DRT::UTILS::intrule_tet_125point_peano;
          break;
        case DRT::Element::hex27:
          rule = DRT::UTILS::intrule_tet_343point_peano;
          break;
        default:
          dserror("no rule defined for this xele");
      }
    }
    else
    {
      dserror("gaussrule only for tetrahedral integration cells defined. Add your stuff here!");
    }
  }

  return rule;
}

/*!
  tell, whether second derivatives have to be calculated
 */
template <DRT::Element::DiscretizationType  DISTYPE>
bool secondDerivativesAvailable()
{
  bool hoel = true;
  switch (DISTYPE)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
  case DRT::Element::pyramid5:
    hoel = true;
    break;
  case DRT::Element::tet4:
    hoel = false;
    break;
  default:
    dserror("DISTYPE unknown!");
  }
  return hoel;
}

//! decide, whether second derivatives can occur
template <DRT::Element::DiscretizationType  DISTYPE>
bool isHigherOrderElement()
{
  bool hoel = true;
  switch (DISTYPE)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::wedge15:
    hoel = true;
    break;
  case DRT::Element::tet4:
  case DRT::Element::wedge6:
  case DRT::Element::pyramid5:
    hoel = false;
    break;
  default:
    dserror("DISTYPE unknown!");
  }
  return hoel;
}


} // namespace XFEM

namespace XFLUID
{

//! interpolate from nodal vector array to integration point vector using the shape function
template <class M, class VS>
static LINALG::FixedSizeSerialDenseMatrix<3,1> interpolateVectorFieldToIntPoint(
    const M&  eleVectorField,       ///< array with nodal vector values
    const VS& shp,                  ///< array with nodal shape function
    const int numparam              ///< number of parameters
    )
{
  LINALG::FixedSizeSerialDenseMatrix<3,1> v;
  const int nsd = 3;
  for (int isd = 0; isd < nsd; ++isd)
  {
      v(isd) = 0.0;
      for (int iparam = 0; iparam < numparam; ++iparam)
          v(isd) += eleVectorField(isd,iparam)*shp(iparam);
  }
  return v;
}

} // namespace XFLUID

#endif
#endif  // #ifdef CCADISCRET
