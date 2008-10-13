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

#include "enrichment_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_geometry/intersection_service.H"
#include "interfacexfsi.H"
#include "../drt_geometry/blitz_tiny_operation.H"
#include "../drt_geometry/coordinate_transformation.H"


class DRT::Discretization;


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
    const DRT::Element&                        ele,
    const RCP<XFEM::InterfaceHandle>       ih,
    const XFEM::ElementDofManager&             dofman,
    const XFEM::PHYSICS::Field                 field,
    const std::map<XFEM::Enrichment, double>&  enrvals,
    const BlitzVec&                            funct,
    BlitzVec&                                  enr_funct
    )
{
    const int* nodeids = ele.NodeIds();
    
    int dofcounter = 0;
    for (int inode=0; inode < ele.NumNode(); ++inode)
    {
        const int gid = nodeids[inode];
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
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec&                       funct,
        const BlitzMat&                       derxy,
        BlitzVec&                             enr_funct,
        BlitzMat&                             enr_derxy
        )
{
    const int* nodeids = ele.NodeIds();
    
    int dofcounter = 0;
    for (int inode=0; inode < ele.NumNode(); ++inode)
    {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                for (int isd = 0; isd < 3; ++isd)
                {
                  enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
                }
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
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec&                       funct,
        const BlitzMat&                       derxy,
        const BlitzMat&                       derxy2,
        BlitzVec&                             enr_funct,
        BlitzMat&                             enr_derxy,
        BlitzMat&                             enr_derxy2
        )
{
    const int* nodeids = ele.NodeIds();
    
    int dofcounter = 0;
    for (int inode=0; inode < ele.NumNode(); ++inode)
    {
        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);
        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
                enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
            if (enrfield->getField() == field)
            {
                const double enrval = enrvals.find(enrfield->getEnrichment())->second;
                enr_funct(dofcounter) = funct(inode) * enrval;
                for (int isd = 0; isd < 3; ++isd)
                {
                  enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
                }
                for (int isd = 0; isd < 6; ++isd)
                {
                  enr_derxy2(isd,dofcounter) = derxy2(isd,inode) * enrval;
                }
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
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec&                       funct,
        BlitzVec&                             enr_funct
        )
{
    int dofcounter = 0;
    
    const std::set<XFEM::FieldEnr>& enrfieldset = dofman.getEnrichedFieldsPerEleField(field);
    const DRT::Element::DiscretizationType distype = dofman.getDisTypePerField(field);
    const int numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
    dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      for (int inode = 0; inode < numvirtualnode; ++inode)
      {
        const double enrval = enrvals.find(enrfield->getEnrichment())->second;
        enr_funct(dofcounter) = funct(inode) * enrval;
        dofcounter += 1;
      }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

//
// For a given situation compute the enriched shape functions
// 
void XFEM::ComputeEnrichedElementShapefunction(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&        dofman,
        const XFEM::PHYSICS::Field            field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec&                       funct,
        const BlitzMat&                       derxy,
        BlitzVec&                             enr_funct,
        BlitzMat&                             enr_derxy
        )
{
    int dofcounter = 0;

    const std::set<XFEM::FieldEnr>& enrfieldset = dofman.getEnrichedFieldsPerEleField(field);
    const DRT::Element::DiscretizationType distype = dofman.getDisTypePerField(field);
    const int numvirtualnode = DRT::UTILS::getNumberOfElementNodes(distype);
    dsassert(enrfieldset.size() > 0, "empty enrfieldset not allowed at this point!");
    for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
    {
      for (int inode = 0; inode < numvirtualnode; ++inode)
      {
        const double enrval = enrvals.find(enrfield->getEnrichment())->second;
        enr_funct(dofcounter) = funct(inode) * enrval;
        for (int isd = 0; isd < 3; ++isd)
        {
          enr_derxy(isd,dofcounter) = derxy(isd,inode) * enrval;
        }
        dofcounter += 1;
      }
    }
    dsassert(dofcounter == dofman.NumDofPerField(field), "mismatch in information from eledofmanager!");
}

void XFEM::computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
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

  std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
        ih,
        dofman.getUniqueEnrichments(),
        cellcenterpos,
        XFEM::Enrichment::approachUnknown));
  
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
    XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, enrvals, funct, enr_funct);
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
  const GEO::DomainIntCell& cell,
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
  
  std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
        ih,
        dofman.getUniqueEnrichments(),
        cellcenterpos,
        XFEM::Enrichment::approachUnknown));

  cellvalues = 0.0;
  for (int incn = 0; incn < cell.NumNode(); ++incn)
  {
    const int numparam  = dofman.NumDofPerField(field);
    if (numparam == 0)
    {
      continue;
    }
    
    const DRT::Element::DiscretizationType eleval_distype = dofman.getDisTypePerField(field);
    const int numvirtnode = DRT::UTILS::getNumberOfElementNodes(eleval_distype);
    if (numvirtnode != numparam) dserror("bug");
    
    BlitzVec funct(numparam);
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,incn),
      (*nodalPosXiDomain)(1,incn),
      (*nodalPosXiDomain)(2,incn),
      eleval_distype);

    BlitzVec enr_funct(numparam);
    XFEM::ComputeEnrichedElementShapefunction(ele, ih, dofman, field, enrvals, funct, enr_funct);
    // interpolate value
    for (int iparam = 0; iparam < numparam; ++iparam)
    {
      cellvalues(incn) += elementvalues(iparam) * enr_funct(iparam);
    }
  }
  return;
}


void XFEM::computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&      cellvalues
  )
{
  const BlitzMat* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const BlitzVec3 cellcenterpos(cell.GetPhysicalCenterPosition(ele));
  flush(cout);
  std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
        ih,
        dofman.getUniqueEnrichments(),
        cellcenterpos,
        XFEM::Enrichment::approachUnknown));
  
  cellvalues = 0.0;
  for (int incn = 0; incn < cell.NumNode(); ++incn)
  {
    const int numparam  = dofman.NumDofPerField(field);
    if (numparam == 0)
    {
      continue;
    }
    
    const DRT::Element::DiscretizationType eleval_distype = dofman.getDisTypePerField(field);
    const int numvirtnode = DRT::UTILS::getNumberOfElementNodes(eleval_distype);
    if (numvirtnode != numparam) dserror("bug");
    
    BlitzVec funct(numparam);
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,incn),
      (*nodalPosXiDomain)(1,incn),
      (*nodalPosXiDomain)(2,incn),
      eleval_distype);

    BlitzVec enr_funct(numparam);
    XFEM::ComputeEnrichedElementShapefunction(ele, ih, dofman, field, enrvals, funct, enr_funct);
    // interpolate value
    for (int iparam = 0; iparam < numparam; ++iparam)
    {
      for (int ientry = 0; ientry < 9; ++ientry)
      {
        cellvalues(ientry,incn) += elementvalues(ientry,iparam) * enr_funct(iparam);        
      }
    }
  }
  return;
}


void XFEM::computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&      cellvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);
  const int nsd = 3;

  const BlitzMat* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());
  
  BlitzMat xyz_cell(3,nen_cell);
  cell.NodalPosXYZ(ele, xyz_cell);

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const blitz::TinyVector<double,3> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
        ih,
        dofman.getUniqueEnrichments(),
        cellcenterpos,
        XFEM::Enrichment::approachUnknown));
  
  // cell corner nodes
  //const BlitzMat cellnodeposvectors = cell.NodalPosXYZ(ele);
  BlitzVec enr_funct(numparam);
  //BlitzVec funct(DRT::UTILS::getNumberOfElementNodes(ele.Shape()));
  static BlitzVec funct(27);
  cellvalues = 0.0;
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,inen),
      (*nodalPosXiDomain)(1,inen),
      (*nodalPosXiDomain)(2,inen),
      ele.Shape());
    if (cell.Shape() == DRT::Element::tet4 or cell.Shape() == DRT::Element::tet10)
    {
      XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, enrvals, funct, enr_funct);
    }
    else
    {
      static BlitzVec3 actpos;
      actpos(0) = xyz_cell(0,inen);
      actpos(1) = xyz_cell(1,inen);
      actpos(2) = xyz_cell(2,inen);
      XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, enrvals, funct, enr_funct);
    }
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
  const GEO::BoundaryIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&      cellvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);
  const int nsd = 3;

  const BlitzMat* nodalPosXiDomain(cell.NodalPosXiDomainBlitz());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const blitz::TinyVector<double,3> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  std::map<XFEM::Enrichment, double> enrvals(computeEnrvalMap(
        ih,
        dofman.getUniqueEnrichments(),
        cellcenterpos,
        XFEM::Enrichment::approachFromPlus));
  
  // cell corner nodes
  //const BlitzMat cellnodeposvectors = cell.NodalPosXYZ(ele);
  BlitzVec enr_funct(numparam);
  //BlitzVec funct(DRT::UTILS::getNumberOfElementNodes(ele.Shape()));
  static BlitzVec funct(27);
  cellvalues = 0.0;
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      (*nodalPosXiDomain)(0,inen),
      (*nodalPosXiDomain)(1,inen),
      (*nodalPosXiDomain)(2,inen),
      ele.Shape());

    XFEM::ComputeEnrichedNodalShapefunction(ele, ih, dofman, field, enrvals, funct, enr_funct);
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


/*!
  Calculate ratio between fictitious element size and normal size
  */
template <DRT::Element::DiscretizationType DISTYPE>
double DomainCoverageRatioT(
        const DRT::Element&               ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
    
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    
    // dimension for 3d fluid element
    const int nsd = 3;
    
    // get node coordinates of the current element
    static blitz::TinyMatrix<double,nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(&ele, xyze);
    
    //double 
    double area_ele  = 0.0;
    double area_fict = 0.0;
    
    // information about domain integration cells
    const GEO::DomainIntCells&  domainIntCells(ih.GetDomainIntCells(ele.Id(),DISTYPE));
    // loop over integration cells
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {

      DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
      switch (cell->Shape())
      {
        case DRT::Element::hex8:
        case DRT::Element::hex20:
        case DRT::Element::hex27:
        {
          gaussrule = DRT::UTILS::intrule_hex_8point;
          break;
        }
        case DRT::Element::tet4:
        case DRT::Element::tet10:
        {
          gaussrule = DRT::UTILS::intrule_tet_4point;
          break;
        }
        default:
          dserror("add your element type here...");
      }
      
        // gaussian points
        const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

        // integration loop
        for (int iquad=0; iquad<intpoints.nquad; ++iquad)
        {
            // coordinates of the current integration point in cell coordinates \eta
            static GEO::PosEtaDomain pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static GEO::PosXiDomain posXiDomain;
            GEO::mapEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain);
            
            // shape functions and their first derivatives
            static blitz::TinyVector<double,numnode> funct;
            static blitz::TinyMatrix<double,nsd,numnode> deriv;
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
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele.Id(), det);
            }

            // inverse of jacobian
            static BlitzMat3x3 xji;
            GEO::Inverse3x3(xjm, det, xji);

            // compute global derivates
            static blitz::TinyMatrix<double,nsd,numnode> derxy;
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
            
            const BlitzVec3 cellcenter(cell->GetPhysicalCenterPosition(ele));
                        
            const int label = ih.PositionWithinConditionNP(cellcenter);
            
            if (label != 0)
            {
              area_fict += fac;
            }
            
        } // end loop over gauss points
    } // end loop over integration cells

    return area_fict / area_ele;
}

double XFEM::DomainCoverageRatio(
        const DRT::Element&               ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return DomainCoverageRatioT<DRT::Element::hex8>(ele,ih);
    case DRT::Element::hex20:
      return DomainCoverageRatioT<DRT::Element::hex20>(ele,ih);
    case DRT::Element::hex27:
      return DomainCoverageRatioT<DRT::Element::hex27>(ele,ih);
    default:
      dserror("add you distype here...");
      exit(1);
  }
}


/*!
  Calculate ratio between fictitious element size and normal size
  */
template <DRT::Element::DiscretizationType DISTYPE>
double BoundaryCoverageRatioT(
        const DRT::Element&               ele,           ///< the element whose boundary ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
  
  double area_fict = 0.0;
  
  // information about boundary integration cells
  const GEO::BoundaryIntCells& boundaryIntCells = ih.GetBoundaryIntCells(ele.Id());
  
  double base_area = 0.0;
  if (DISTYPE == DRT::Element::tet10 or DISTYPE == DRT::Element::tet4)
  {
    base_area = 0.5;
  }
  else if (DISTYPE == DRT::Element::hex8 or DISTYPE == DRT::Element::hex20 or DISTYPE == DRT::Element::hex27)
  {
    base_area = 4.0;
  }
  else
  {
    dserror("think about it. factor at the end of this function needs another values");
  }
  
  // loop over boundary integration cells
  for (GEO::BoundaryIntCells::const_iterator cell = boundaryIntCells.begin(); cell != boundaryIntCells.end(); ++cell)
  {
    
    DRT::UTILS::GaussRule2D gaussrule = DRT::UTILS::intrule2D_undefined;
    switch (cell->Shape())
    {
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      gaussrule = DRT::UTILS::intrule_tri_1point;
      break;
    }
    default:
      dserror("add your element type here...");
    }
    
    // gaussian points
    const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);
    
    // get the right boundary element
//    const DRT::Element* boundaryele = ih.GetBoundaryEle(cell->GetSurfaceEleGid());
    //cout << (*boundaryele) << endl;
//    const int numnode_boundary = boundaryele->NumNode();
    //        cout << "numnode_boundary: " << numnode_boundary << endl;
    
    // get current node coordinates
//    const std::map<int,blitz::TinyVector<double,3> >* positions = ih.cutterposnp();
//    const BlitzMat xyze_boundary(GEO::getCurrentNodalPositions(boundaryele, *positions));
    
    const BlitzMat* nodalpos_xi_domain = cell->NodalPosXiDomainBlitz();
    const int numnode_cell = cell->NumNode();
    
    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in cell coordinates \eta^\boundary
      GEO::PosEtaBoundary pos_eta_boundary;
      pos_eta_boundary(0) = intpoints.qxg[iquad][0];
      pos_eta_boundary(1) = intpoints.qxg[iquad][1];
      
      //            cout << pos_eta_boundary << endl;
      GEO::PosXiDomain posXiDomain;
      mapEtaBToXiD(*cell, pos_eta_boundary, posXiDomain);
      //            cout << cell->toString() << endl;
      //            cout << posXiDomain << endl;
      
      // shape functions and their first derivatives
      BlitzVec funct_boundary(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
      DRT::UTILS::shape_function_2D(funct_boundary, pos_eta_boundary(0),pos_eta_boundary(1),cell->Shape());
      BlitzMat deriv_boundary(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, pos_eta_boundary(0),pos_eta_boundary(1),cell->Shape());
      
      // shape functions and their first derivatives
      static BlitzVec funct(DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement);
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
      // get jacobian matrix d x / d \xi  (3x2)
      static BlitzMat3x2 dxyzdrs;
      //dxyzdrs = blitz::sum(xyze_boundary(i,k)*deriv_boundary(j,k),k);
      for (int isd = 0; isd < 3; ++isd)
      {
        for (int j = 0; j < 2; ++j)
        {
          dxyzdrs(isd,j) = 0.0;
          for (int k = 0; k < numnode_cell; ++k)
          {
            dxyzdrs(isd,j) += (*nodalpos_xi_domain)(isd,k)*deriv_boundary(j,k);
          }
        }
      }
      
      // compute covariant metric tensor G for surface element (2x2)
      static BlitzMat2x2 metric;
      //metric = blitz::sum(dxyzdrs(k,i)*dxyzdrs(k,j),k);
      BLITZTINY::MtM_product<2,2,3>(dxyzdrs,dxyzdrs,metric);
      //const BlitzMat metric = computeMetricTensor(xyze_boundary,deriv_boundary);
      
      const double detmetric = sqrt(metric(0,0)*metric(1,1) - metric(0,1)*metric(1,0));
      if (detmetric <= 0.0)
      {
        dserror("negative detmetric! should be a bug!");
      }
      
      const double fac = intpoints.qwgt[iquad]*detmetric;//*detcell;
      if (fac <= 0.0)
      {
        dserror("negative fac! should be a bug!");
      }
      
      area_fict += fac;
      
    } // end loop over gauss points
  } // end loop over integration cells
  
  // scale result by area of one surface of the volume element
  return area_fict / base_area;
}

double XFEM::BoundaryCoverageRatio(
        const DRT::Element&               ele,           ///< the element whose boundary ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return BoundaryCoverageRatioT<DRT::Element::hex8>(ele,ih);
    case DRT::Element::hex20:
      return BoundaryCoverageRatioT<DRT::Element::hex20>(ele,ih);
    case DRT::Element::hex27:
      return BoundaryCoverageRatioT<DRT::Element::hex27>(ele,ih);
    default:
      dserror("add you distype here...");
      exit(1);
  }
}

#endif  // #ifdef CCADISCRET
