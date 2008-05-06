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
#include "../drt_lib/drt_discret.H"


namespace XFEM
{

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
void fillElementUnknownsArrays(
        const XFEM::ElementDofManager& dofman,
        const int numnode,
        const vector<double>& locval,
        const vector<double>& locval_hist,
        BlitzMat& evelnp,
        BlitzMat& evelnp_hist,
        BlitzVec& eprenp,
        BlitzMat& etau,
        BlitzVec& ediscprenp
        )
{
    
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
} // namespace XFEM

#endif  // #ifndef DOF_MANAGEMENT_H
#endif  // #ifdef CCADISCRET
