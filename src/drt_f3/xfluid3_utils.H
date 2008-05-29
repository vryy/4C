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

#include "../drt_xfem/xfem.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/xdofmapcreation.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_xfem/dofkey.H"
#include "../drt_xfem/blitz_tiny_operation.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_timecurve.H"
#include <Epetra_SerialDenseSolver.h>
#include "../drt_xfem/coordinate_transformation.H"


namespace XFEM
{

using namespace XFEM::PHYSICS;
  
/*! get number of unknowns for this field
 * use the dofmanager in the general case,
 * otherwise use a given number numnode
 */
template<XFEM::AssemblyType>
int getNumParam(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        const int)
{
    return dofman.NumDofPerField(field);
}
template<> 
int getNumParam<XFEM::standard_assembly>(
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field     field,
        const int numnode)
{
    return numnode;
}

template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void fillElementUnknownsArrays2(
        const XFEM::ElementDofManager& dofman,
        const vector<double>& locval,
        const vector<double>& locval_hist,
        BlitzMat& evelnp,
        BlitzMat& evelnp_hist,
        BlitzVec& eprenp,
        BlitzMat& etau,
        BlitzVec& ediscprenp
        )
{
    
    const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
    
    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
    const int numparamvely = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Vely, numnode);
    const int numparamvelz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velz, numnode);
    const int numparampres = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Pres, numnode);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const int numparamtauxx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 1);
    const int numparamdiscpres = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 1);
    
    const vector<int>& velxdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Velx));
    const vector<int>& velydof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Vely));
    const vector<int>& velzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Velz));
    const vector<int>& presdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Pres));
    
    for (int iparam=0; iparam<numparamvelx; ++iparam)
    {
        evelnp(     0,iparam) = locval[     velxdof[iparam]];
        evelnp_hist(0,iparam) = locval_hist[velxdof[iparam]];
    }
    for (int iparam=0; iparam<numparamvely; ++iparam)
    {
        evelnp(     1,iparam) = locval[     velydof[iparam]];
        evelnp_hist(1,iparam) = locval_hist[velydof[iparam]];
    }
    for (int iparam=0; iparam<numparamvelz; ++iparam)
    {
        evelnp(     2,iparam) = locval[     velzdof[iparam]];
        evelnp_hist(2,iparam) = locval_hist[velzdof[iparam]];
    }
    for (int iparam=0; iparam<numparampres; ++iparam)
        eprenp(iparam) = locval[presdof[iparam]];
    const bool tauele_unknowns_present = (getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 0) > 0);
    if (tauele_unknowns_present)
    {
        const int numparamtauyy = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyy, 1);
        const int numparamtauzz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauzz, 1);
        const int numparamtauxy = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxy, 1);
        const int numparamtauxz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxz, 1);
        const int numparamtauyz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyz, 1);
        const vector<int>& tauxxdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxx));
        const vector<int>& tauyydof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauyy));
        const vector<int>& tauzzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauzz));
        const vector<int>& tauxydof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxy));
        const vector<int>& tauxzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxz));
        const vector<int>& tauyzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauyz));
        for (int iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = locval[tauxxdof[iparam]];
        for (int iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = locval[tauyydof[iparam]];
        for (int iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = locval[tauzzdof[iparam]];
        for (int iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = locval[tauxydof[iparam]];
        for (int iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = locval[tauxzdof[iparam]];
        for (int iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = locval[tauyzdof[iparam]];
    }
    const bool discpres_unknowns_present = (getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::DiscPres, 0) > 0);
    if (discpres_unknowns_present)
    {
        const vector<int>& discpresdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::DiscPres));
        for (int iparam=0; iparam<numparamdiscpres; ++iparam)   ediscprenp(iparam) = locval[discpresdof[iparam]];
    }
}

template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
void fillElementUnknownsArrays3(
        const XFEM::ElementDofManager& dofman,
        const vector<double>& locval,
        const vector<double>& locval_hist,
        BlitzMat& evelnp,
        BlitzMat& evelnp_hist,
        BlitzVec& eprenp,
        BlitzMat& etau
        )
{
    
    const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
    
    // number of parameters for each field (assumed to be equal for each velocity component and the pressure)
    const int numparamvelx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velx, numnode);
    const int numparamvely = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Vely, numnode);
    const int numparamvelz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Velz, numnode);
    const int numparampres = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Pres, numnode);
    // put one here to create arrays of size 1, since they are not needed anyway
    // in the xfem assembly, the numparam is determined by the dofmanager
    const int numparamtauxx = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 1);
    
    const vector<int>& velxdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Velx));
    const vector<int>& velydof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Vely));
    const vector<int>& velzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Velz));
    const vector<int>& presdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Pres));
    
    for (int iparam=0; iparam<numparamvelx; ++iparam)
    {
        evelnp(     0,iparam) = locval[     velxdof[iparam]];
        evelnp_hist(0,iparam) = locval_hist[velxdof[iparam]];
    }
    for (int iparam=0; iparam<numparamvely; ++iparam)
    {
        evelnp(     1,iparam) = locval[     velydof[iparam]];
        evelnp_hist(1,iparam) = locval_hist[velydof[iparam]];
    }
    for (int iparam=0; iparam<numparamvelz; ++iparam)
    {
        evelnp(     2,iparam) = locval[     velzdof[iparam]];
        evelnp_hist(2,iparam) = locval_hist[velzdof[iparam]];
    }
    for (int iparam=0; iparam<numparampres; ++iparam)
        eprenp(iparam) = locval[presdof[iparam]];
    const bool tauele_unknowns_present = (getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxx, 0) > 0);
    if (tauele_unknowns_present)
    {
        const int numparamtauyy = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyy, 1);
        const int numparamtauzz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauzz, 1);
        const int numparamtauxy = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxy, 1);
        const int numparamtauxz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauxz, 1);
        const int numparamtauyz = getNumParam<ASSTYPE>(dofman, XFEM::PHYSICS::Tauyz, 1);
        const vector<int>& tauxxdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxx));
        const vector<int>& tauyydof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauyy));
        const vector<int>& tauzzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauzz));
        const vector<int>& tauxydof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxy));
        const vector<int>& tauxzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauxz));
        const vector<int>& tauyzdof(dofman.LocalDofPosPerField(XFEM::PHYSICS::Tauyz));
        for (int iparam=0; iparam<numparamtauxx; ++iparam)   etau(0,iparam) = locval[tauxxdof[iparam]];
        for (int iparam=0; iparam<numparamtauyy; ++iparam)   etau(1,iparam) = locval[tauyydof[iparam]];
        for (int iparam=0; iparam<numparamtauzz; ++iparam)   etau(2,iparam) = locval[tauzzdof[iparam]];
        for (int iparam=0; iparam<numparamtauxy; ++iparam)   etau(3,iparam) = locval[tauxydof[iparam]];
        for (int iparam=0; iparam<numparamtauxz; ++iparam)   etau(4,iparam) = locval[tauxzdof[iparam]];
        for (int iparam=0; iparam<numparamtauyz; ++iparam)   etau(5,iparam) = locval[tauyzdof[iparam]];
    }
}


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
  const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
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
    for (int jnode=0; jnode<numnode; ++jnode)
    {
      nodes[jnode]->GetCondition("VolumeNeumann",myneumcond);

      // get values and switches from the condition
      const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
      const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );

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
template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2, class M3, class M4>
void gder2(
        const M1& xjm,
        const M2& derxy,
        const M3& deriv2,
        const M4& xyze,
        BlitzMat& derxy2_)
{
    const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;

    // initialize and zero out everything
    static Epetra_SerialDenseMatrix bm(6,6);

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

    //const BlitzMat xder2(blitz::sum(deriv2(i,k)*xyze(j,k),k));
    static blitz::TinyMatrix<double,6,3> xder2;
    BLITZTINY::MMt_product<6,3,numnode>(deriv2,xyze,xder2);

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

    //derxy2_ = deriv2 - blitz::sum(xder2(i,k)*derxy(k,j),k);
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < numnode; ++j)
        {
            derxy2_(i,j) = deriv2(i,j);
            for (int k = 0; k < 3; ++k)
            {
                derxy2_(i,j) -= xder2(i,k)*derxy(k,j);
            }
        }
    }

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
    Epetra_SerialDenseMatrix ederxy2(View,derxy2_.data(),6,6,numnode);

    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(bm);

    // No need for a separate rhs. We assemble the rhs to the solution
    // vector. The solver will destroy the rhs and return the solution.
    solver.SetVectors(ederxy2,ederxy2);
    solver.Solve();

    return;
}


/*!
 \brief get optimal gaussrule for given discretization type
        template version
 \return enum with optimal gausrule for this element
 */
template <DRT::Element::DiscretizationType  distype>
inline DRT::UTILS::GaussRule3D getOptimalVelocityGaussrule()
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
        rule = DRT::UTILS::intrule_tet_5point;
        break;
    default:
        dserror("unknown distype for gaussrule initialization");
  }
  return rule;
}

/*!
 \brief get optimal gaussrule for given discretization type
        template version
 \return enum with optimal gausrule for this element
 */
template <DRT::Element::DiscretizationType  distype, XFEM::AssemblyType asstype>
inline DRT::UTILS::GaussRule3D getXFEMGaussrule(const bool intersected)
{
    DRT::UTILS::GaussRule3D rule = getOptimalVelocityGaussrule<distype>();
    if (asstype == XFEM::xfem_assembly)
    {
        if (intersected)
        {
          // set integration rule for tetrahedral integrtion cells
          // attention: use "open" rules -> no integration points on boundary of element
          // don't use gausspoint rules that have gausspoints on the element boundary 
          // (so called closed rules)
          // -> enrichment values might be undefined
            if (distype == DRT::Element::hex8)
            {
                rule = DRT::UTILS::intrule_tet_24point;
            }
            else
            {
              dserror("we need higher order gauss rules here, but they still need to be implemented :-(");
                rule = DRT::UTILS::intrule_tet_24point;
            }
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

  /*!
    Calculate inverse of a matrix

   */
template<class M1, class M2>
void Inverse3x3(const M1& xjm_, const double det, M2& inverse)
{
    inverse(0,0) = (  xjm_(1,1)*xjm_(2,2) - xjm_(2,1)*xjm_(1,2))/det;
    inverse(1,0) = (- xjm_(1,0)*xjm_(2,2) + xjm_(2,0)*xjm_(1,2))/det;
    inverse(2,0) = (  xjm_(1,0)*xjm_(2,1) - xjm_(2,0)*xjm_(1,1))/det;
    inverse(0,1) = (- xjm_(0,1)*xjm_(2,2) + xjm_(2,1)*xjm_(0,2))/det;
    inverse(1,1) = (  xjm_(0,0)*xjm_(2,2) - xjm_(2,0)*xjm_(0,2))/det;
    inverse(2,1) = (- xjm_(0,0)*xjm_(2,1) + xjm_(2,0)*xjm_(0,1))/det;
    inverse(0,2) = (  xjm_(0,1)*xjm_(1,2) - xjm_(1,1)*xjm_(0,2))/det;
    inverse(1,2) = (- xjm_(0,0)*xjm_(1,2) + xjm_(1,0)*xjm_(0,2))/det;
    inverse(2,2) = (  xjm_(0,0)*xjm_(1,1) - xjm_(1,0)*xjm_(0,1))/det;
    return;
}


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


/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE,
          XFEM::AssemblyType ASSTYPE>
double AreaRatio(
        const DRT::Element*               ele,           ///< the element those matrix is calculated
        const RCP<XFEM::InterfaceHandle>  ih             ///< connection to the interface handler
        )
{
    
    // number of nodes for element
    const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
    
    // dimension for 3d fluid element
    const int nsd = 3;
    
    // get node coordinates of the current element
    static blitz::TinyMatrix<double,nsd,numnode> xyze;
    DRT::UTILS::fillPositionArray<DISTYPE>(ele, xyze);
    
    const DRT::UTILS::GaussRule3D gaussrule = getXFEMGaussrule<DISTYPE,ASSTYPE>(ih->ElementIntersected(ele->Id()));
    
  
    //double 
    double area_ele  = 0.0;
    double area_fict = 0.0;
    
    // information about domain integration cells
    const XFEM::DomainIntCells&  domainIntCells(ih->GetDomainIntCells(ele->Id(),DISTYPE));
    // loop over integration cells
    for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {

        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            static XFEM::PosEtaDomain pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static XFEM::PosXiDomain posXiDomain;
            XFEM::mapEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = XFEM::detEtaToXi3D<ASSTYPE>(*cell, pos_eta_domain);
            
            // shape functions and their first derivatives
            static BlitzVec funct(numnode);
            static BlitzMat deriv(3, numnode, blitz::ColumnMajorArray<2>());
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
            DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
            // position of the gausspoint in physical coordinates
            static BlitzVec3 gauss_pos_xyz;
            BLITZTINY::MV_product<3,numnode>(xyze,funct,gauss_pos_xyz);
      
            // get transposed of the jacobian matrix d x / d \xi
            //BlitzMat xjm(3,3);
            //xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);
            static BlitzMat3x3 xjm;
            BLITZTINY::MMt_product<3,3,numnode>(deriv,xyze,xjm);

            const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                               xjm(0,1)*xjm(1,2)*xjm(2,0)+
                               xjm(0,2)*xjm(1,0)*xjm(2,1)-
                               xjm(0,2)*xjm(1,1)*xjm(2,0)-
                               xjm(0,0)*xjm(1,2)*xjm(2,1)-
                               xjm(0,1)*xjm(1,0)*xjm(2,2);
            const double fac = intpoints.qwgt[iquad]*det*detcell;

            if (det < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
            }

            // inverse of jacobian
            static BlitzMat xji(3,3);
            Inverse3x3(xjm, det, xji);

            // compute global derivates
            BlitzMat derxy(3, numnode);
            //BlitzMat derxy_stress(3, DRT::UTILS::getNumberOfElementNodes(stressdistype),blitz::ColumnMajorArray<2>());
            //BlitzMat derxy_discpres(3, DRT::UTILS::getNumberOfElementNodes(discpresdistype),blitz::ColumnMajorArray<2>());
            //derxy          = blitz::sum(xji(i,k)*deriv(k,j),k);
            for (int isd = 0; isd < nsd; ++isd)
            {
              for (int inode = 0; inode < numnode; ++inode)
              {
                derxy(isd,inode) = 0.0;
                for (int jsd = 0; jsd < nsd; ++jsd)
                {
                   derxy(isd,inode) += xji(isd,jsd)*deriv(jsd,inode);
                }
              }
            }
            
            area_ele += fac;
            
            const BlitzVec3 cellcenter(cell->GetPhysicalCenterPosition(*ele));
                        
            map<int,bool> posInCondition;
            PositionWithinCondition(cellcenter, *ih, posInCondition);
            bool in_solid = false;
            for (map<int,bool>::const_iterator p = posInCondition.begin(); p != posInCondition.end(); ++p)
            {
              if (p->second == true)
              {
                in_solid = true;
                break;
              }
            }
            
            if (in_solid)
            {
              area_fict += fac;
            }
            
        } // end loop over gauss points
    } // end loop over integration cells

    return pow(area_fict,3) / (area_ele);
}

} // namespace XFEM

#endif
#endif  // #ifdef CCADISCRET
