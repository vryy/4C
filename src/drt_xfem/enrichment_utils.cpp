/*!
\file enrichment_utils.cpp

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET

#include <string>
#include <sstream>
#include <blitz/array.h>
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "intersection_service.H"
#include "xfem.H"
#include "physics.H"
#include "enrichment_utils.H"
#include "dof_management.H"
#include "interface.H"




std::map<XFEM::Enrichment, double> XFEM::computeEnrvalMap(
        const RCP<XFEM::InterfaceHandle>      ih,
        const std::set<XFEM::Enrichment>&     enrset,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection
        )
{
    std::map<XFEM::Enrichment, double> enrvals;
    //TODO: achtung: bei mehreren enrichments ist approach from plus nicht mehr so einfach
    for (std::set<XFEM::Enrichment>::const_iterator enriter =
        enrset.begin(); enriter != enrset.end(); ++enriter)
    {
        const double enrval = enriter->EnrValue(actpos, *ih, approachdirection);
        enrvals.insert(make_pair((*enriter), enrval));
    }
    return enrvals;
}



//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedNodalShapefunction(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>      ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection,
        const BlitzVec&                       funct,
        BlitzVec&                             enr_funct
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection));
    
    const DRT::Node*const* nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedNodalShapefunction(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>      ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection,
        const BlitzVec&                       funct,
        const BlitzMat&                       derxy,
        BlitzVec&                             enr_funct,
        BlitzMat&                             enr_derxy
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection));
    
    blitz::Range _  = blitz::Range::all();
    
    const DRT::Node*const* nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedNodalShapefunction(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>      ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection,
        const BlitzVec&                       funct,
        const BlitzMat&                       derxy,
        const BlitzMat&                       derxy2,
        BlitzVec&                             enr_funct,
        BlitzMat&                             enr_derxy,
        BlitzMat&                             enr_derxy2
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection));
    
    blitz::Range _  = blitz::Range::all();
    
    const DRT::Node*const* nodes = ele.Nodes();
    
    int dofcounter = 0;
    for (int inode=0; inode<ele.NumNode(); ++inode)
    {
        const int gid = nodes[inode]->Id();
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                enr_derxy2(_,dofcounter) = derxy2(_,inode) * enrval;
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}


//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedElementShapefunction(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>      ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection,
        const BlitzVec&                       funct,
        BlitzVec&                             enr_funct
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection));
    
    int dofcounter = 0;
    for (int inode = 0; inode < dofman.NumVirtualNodes(); ++inode)
    {
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerVirtualElementNode(inode);
        dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                dofcounter += 1;                
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedElementShapefunction(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>      ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection,
        const BlitzVec&                       funct,
        const BlitzMat&                       derxy,
        BlitzVec&                             enr_funct,
        BlitzMat&                             enr_derxy
        )
{
    
    // compute enrichment values for all available enrichemnts in this dofmap (saves lots of time)
    std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
            ih,
            dofman.getUniqueEnrichments(),
            actpos,
            approachdirection));
    
    blitz::Range _  = blitz::Range::all();
    
    int dofcounter = 0;
    for (int inode = 0; inode < dofman.NumVirtualNodes(); ++inode)
    {
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerVirtualElementNode(inode);
        dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                enr_derxy(_,dofcounter) = derxy(_,inode) * enrval;
                dofcounter += 1;
            }
        }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

void XFEM::computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzVec& elementvalues,
  BlitzVec&      cellvalues
  )
{
  const BlitzMat* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const BlitzVec3 cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  cellvalues = 0.0;
  for (int inen = 0; inen < cell.NumNode(); ++inen)
  {
    BlitzVec funct(ele.NumNode());
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,inen),
      (*nodalPosXiDomain)(1,inen),
      (*nodalPosXiDomain)(2,inen),
      ele.Shape());

    const int numparam  = dofman.NumDofPerField(field);
    BlitzVec enr_funct(numparam);
    XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, cellcenterpos, XFEM::Enrichment::approachUnknown, funct, enr_funct);
    // interpolate value
    for (int iparam = 0; iparam < numparam; ++iparam)
    {
      cellvalues(inen) += elementvalues(iparam) * enr_funct(iparam);
    }
  }
  return;
}


void XFEM::computeScalarCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzVec& elementvalues,
  BlitzVec&      cellvalues
  )
{
  const BlitzMat* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const BlitzVec3 cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  for (int incn = 0; incn < cell.NumNode(); ++incn)
  {
    BlitzVec funct(ele.NumNode());
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,incn),
      (*nodalPosXiDomain)(1,incn),
      (*nodalPosXiDomain)(2,incn),
      ele.Shape());

    const int numparam  = dofman.NumDofPerField(field);
    if (numparam == 0)
      continue;
    BlitzVec enr_funct(numparam);
    XFEM::ComputeEnrichedElementShapefunction(ele, ih, dofman, field, cellcenterpos, XFEM::Enrichment::approachUnknown, funct, enr_funct);
    // interpolate value
    cellvalues(incn) = 0.0;
    for (int iparam = 0; iparam < numparam; ++iparam)
    {
      cellvalues(incn) += elementvalues(iparam) * enr_funct(iparam);
    }
  }
  return;
}


void XFEM::computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const blitz::Array<double,2>& elementvalues,
  blitz::Array<double,2>&      cellvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);
  const int nsd = 3;

  const blitz::Array<double,2>* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const blitz::TinyVector<double,3> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  // cell corner nodes
  //const blitz::Array<double,2> cellnodeposvectors = cell.NodalPosXYZ(ele);
  blitz::Array<double,1> enr_funct(numparam);
  //blitz::Array<double,1> funct(DRT::UTILS::getNumberOfElementNodes(ele.Shape()));
  static blitz::Array<double,1> funct(27);
  cellvalues = 0.0;
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,inen),
      (*nodalPosXiDomain)(1,inen),
      (*nodalPosXiDomain)(2,inen),
      ele.Shape());

    XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, cellcenterpos, XFEM::Enrichment::approachUnknown, funct, enr_funct);
    // interpolate value
    for (int iparam = 0; iparam < numparam; ++iparam)
    {
      for (int isd = 0; isd < nsd; ++isd)
      {
        cellvalues(isd,inen) += elementvalues(isd,iparam) * enr_funct(iparam);        
      }
    }
  }
  return;
}

void XFEM::computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::BoundaryIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const blitz::Array<double,2>& elementvalues,
  blitz::Array<double,2>&      cellvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);
  const int nsd = 3;

  const blitz::Array<double,2>* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const blitz::TinyVector<double,3> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  // cell corner nodes
  //const blitz::Array<double,2> cellnodeposvectors = cell.NodalPosXYZ(ele);
  blitz::Array<double,1> enr_funct(numparam);
  //blitz::Array<double,1> funct(DRT::UTILS::getNumberOfElementNodes(ele.Shape()));
  static blitz::Array<double,1> funct(27);
  cellvalues = 0.0;
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,inen),
      (*nodalPosXiDomain)(1,inen),
      (*nodalPosXiDomain)(2,inen),
      ele.Shape());

    XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, cellcenterpos, XFEM::Enrichment::approachFromPlus, funct, enr_funct);
    // interpolate value
    for (int iparam = 0; iparam < numparam; ++iparam)
    {
      for (int isd = 0; isd < nsd; ++isd)
      {
        cellvalues(isd,inen) += elementvalues(isd,iparam) * enr_funct(iparam);        
      }
    }
  }
  return;
}

#endif  // #ifdef CCADISCRET
