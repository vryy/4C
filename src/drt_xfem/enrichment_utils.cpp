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
#include "coordinate_transformation.H"


XFEM::ElementEnrichmentValues::ElementEnrichmentValues(
        const DRT::Element&                   ele,
        const RCP<XFEM::InterfaceHandle>&     ih,               ///< interface information
        const XFEM::ElementDofManager&        dofman,
        const LINALG::Matrix<3,1>&            actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection
        ) :
          ele_(ele),
          dofman_(dofman) 
{
    enrvals_.clear();
    const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
    //TODO: achtung: bei mehreren enrichments ist approach from plus nicht mehr so einfach
    for (std::set<XFEM::Enrichment>::const_iterator enriter =
        enrset.begin(); enriter != enrset.end(); ++enriter)
    {
        const double enrval = enriter->EnrValue(actpos, *ih, approachdirection);
//        enrvals_.insert(make_pair((*enriter), enrval));
        enrvals_[*enriter] = enrval;
    }
    return;
}

void XFEM::computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const LINALG::SerialDenseVector& elementvalues,
  LINALG::SerialDenseVector&      cellvalues
  )
{
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.NodalPosXiDomain());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const LINALG::Matrix<3,1> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  const XFEM::ElementEnrichmentValues enrvals(
        ele,
        ih,
        dofman,
        cellcenterpos,
        XFEM::Enrichment::approachUnknown);
  
  cellvalues.Zero();
  for (int inen = 0; inen < cell.NumNode(); ++inen)
  {
    LINALG::SerialDenseVector funct(ele.NumNode());
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      nodalPosXiDomain(0,inen),
      nodalPosXiDomain(1,inen),
      nodalPosXiDomain(2,inen),
      ele.Shape());

    const int numparam  = dofman.NumDofPerField(field);
    LINALG::SerialDenseVector enr_funct(numparam);
    enrvals.ComputeEnrichedNodalShapefunction(field, funct, enr_funct);
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
  const LINALG::SerialDenseVector& elementvalues,
  LINALG::SerialDenseVector&      cellvalues
  )
{
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.NodalPosXiDomain());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const LINALG::Matrix<3,1> cellcenterpos(cell.GetPhysicalCenterPosition(ele));
  
  const XFEM::ElementEnrichmentValues enrvals(
        ele,
        ih,
        dofman,
        cellcenterpos,
        XFEM::Enrichment::approachUnknown);

  cellvalues.Zero();
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
    
    LINALG::SerialDenseVector funct(numparam);
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      nodalPosXiDomain(0,incn),
      nodalPosXiDomain(1,incn),
      nodalPosXiDomain(2,incn),
      eleval_distype);

    LINALG::SerialDenseVector enr_funct(numparam);
    enrvals.ComputeEnrichedElementShapefunction(field, funct, enr_funct);
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
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&       cellvalues
  )
{
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.NodalPosXiDomain());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const LINALG::Matrix<3,1> cellcenterpos(cell.GetPhysicalCenterPosition(ele));
  flush(cout);
  const XFEM::ElementEnrichmentValues enrvals(
        ele,
        ih,
        dofman,
        cellcenterpos,
        XFEM::Enrichment::approachUnknown);
  
  cellvalues.Zero();
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
    
    LINALG::SerialDenseVector funct(numparam);
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      nodalPosXiDomain(0,incn),
      nodalPosXiDomain(1,incn),
      nodalPosXiDomain(2,incn),
      eleval_distype);

    LINALG::SerialDenseVector enr_funct(numparam);
    enrvals.ComputeEnrichedElementShapefunction(field, funct, enr_funct);
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
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&       cellvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);
  const int nsd = 3;

  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.NodalPosXiDomain());
  
  LINALG::SerialDenseMatrix xyz_cell(3,nen_cell);
  cell.NodalPosXYZ(ele, xyz_cell);

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const LINALG::Matrix<3,1> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

  const XFEM::ElementEnrichmentValues enrvals(
        ele,
        ih,
        dofman,
        cellcenterpos,
        XFEM::Enrichment::approachUnknown);
  
  // cell corner nodes
  //const LINALG::SerialDenseMatrix cellnodeposvectors = cell.NodalPosXYZ(ele);
  LINALG::SerialDenseVector enr_funct(numparam);
  //LINALG::SerialDenseVector funct(DRT::UTILS::getNumberOfElementNodes(ele.Shape()));
  static LINALG::SerialDenseVector funct(27);
  cellvalues.Zero();
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      nodalPosXiDomain(0,inen),
      nodalPosXiDomain(1,inen),
      nodalPosXiDomain(2,inen),
      ele.Shape());
    enrvals.ComputeEnrichedNodalShapefunction(field, funct, enr_funct);
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
  const LINALG::SerialDenseMatrix& elementvalues,
  LINALG::SerialDenseMatrix&      cellvalues
  )
{
  const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
  const int numparam  = dofman.NumDofPerField(field);
  const int nsd = 3;

  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.NodalPosXiDomain());

  // if cell node is on the interface, the value is not defined for a jump.
  // however, we approach the interface from one particular side and therefore,
  // -> we use the center of the cell to determine, where we come from
  const LINALG::Matrix<3,1> cellcenterpos(true);

  const XFEM::ElementEnrichmentValues enrvals(
        ele,
        ih,
        dofman,
        cellcenterpos,
        XFEM::Enrichment::approachFromPlus);
  
  // cell corner nodes
  //const LINALG::SerialDenseMatrix cellnodeposvectors = cell.NodalPosXYZ(ele);
  LINALG::SerialDenseVector enr_funct(numparam);
  //LINALG::SerialDenseVector funct(DRT::UTILS::getNumberOfElementNodes(ele.Shape()));
  static LINALG::SerialDenseVector funct(27);
  cellvalues.Zero();
  for (int inen = 0; inen < nen_cell; ++inen)
  {
    // fill shape functions
    DRT::UTILS::shape_function_3D(funct,
      nodalPosXiDomain(0,inen),
      nodalPosXiDomain(1,inen),
      nodalPosXiDomain(2,inen),
      ele.Shape());

    enrvals.ComputeEnrichedNodalShapefunction(field, funct, enr_funct);
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
        const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
    
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    
    // dimension for 3d fluid element
    const int nsd = 3;
    
    // get node coordinates of the current element
    static LINALG::Matrix<nsd,numnode> xyze;
    GEO::fillInitialPositionArray<DISTYPE>(&ele, xyze);
    
    //double 
    double area_ele  = 0.0;
    double area_fict = 0.0;
    
    // information about domain integration cells
    const GEO::DomainIntCells&  domainIntCells(ih.GetDomainIntCells(ele.Id(),DISTYPE));
    // loop over integration cells
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {

      const LINALG::Matrix<3,1> cellcenter(cell->GetPhysicalCenterPosition(ele));
                  
      const int label = ih.PositionWithinConditionNP(cellcenter);
      
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
            static LINALG::Matrix<3,1> pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static LINALG::Matrix<3,1> posXiDomain;
            GEO::mapEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain);
            
            // shape functions and their first derivatives
            static LINALG::Matrix<numnode,1> funct;
            static LINALG::Matrix<nsd,numnode> deriv;
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
            DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
            // get transposed of the jacobian matrix d x / d \xi
            //xjm = deriv(i,k)*xyze(j,k);
            static LINALG::Matrix<3,3> xjm;
            xjm.MultiplyNT(deriv,xyze);

            const double det = xjm.Determinant();
            const double fac = intpoints.qwgt[iquad]*det*detcell;

            if (det < 0.0)
            {
                dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele.Id(), det);
            }

            area_ele += fac;
            
            if (label != 0)
            {
              area_fict += fac;
            }
            
        } // end loop over gauss points
    } // end loop over integration cells

    return area_fict / area_ele;
}

double XFEM::DomainCoverageRatio(
        const DRT::Element&           ele,
        const XFEM::InterfaceHandle&  ih
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
vector<double> DomainCoverageRatioPerNodeT(
        const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
    
    // number of nodes for element
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
    
    // dimension for 3d fluid element
//    const int nsd = 3;
    
    // get node coordinates of the current element
//    static LINALG::Matrix<nsd,numnode> xyze;
//    GEO::fillInitialPositionArray<DISTYPE>(&ele, xyze);
    
    //double 
    double area_ele  = 0.0;
//    double area_fict = 0.0;
    vector<double> portions(numnode,0.0);
    
    // information about domain integration cells
    const GEO::DomainIntCells&  domainIntCells(ih.GetDomainIntCells(ele.Id(),DISTYPE));
    // loop over integration cells
    for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
    {

      const LINALG::Matrix<3,1> cellcenter(cell->GetPhysicalCenterPosition(ele));
                  
      const int label = ih.PositionWithinConditionNP(cellcenter);
      
      DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
      switch (cell->Shape())
      {
        case DRT::Element::hex8:
        {
          gaussrule = DRT::UTILS::intrule_hex_8point;
          break;
        }
        case DRT::Element::hex20:
        case DRT::Element::hex27:
        {
          gaussrule = DRT::UTILS::intrule_hex_27point;
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
            static LINALG::Matrix<3,1> pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static LINALG::Matrix<3,1> posXiDomain;
            GEO::mapEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain);
            
            // shape functions and their first derivatives
            static LINALG::Matrix<numnode,1> funct;
//            static LINALG::Matrix<nsd,numnode> deriv;
            DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
//            DRT::UTILS::shape_function_3D_deriv1(deriv,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
            const double fac = intpoints.qwgt[iquad]*detcell;

            area_ele += fac;
            
            if (label == 0)
            {
              for (int inode = 0;inode < numnode;++inode)
              {
                portions[inode] += funct(inode) * fac;
              }
            }
            
        } // end loop over gauss points
    } // end loop over integration cells

    for (int inode = 0;inode < numnode;++inode)
    {
      portions[inode] /= area_ele;
    }    
    
    return portions;
}

vector<double> XFEM::DomainCoverageRatioPerNode(
        const DRT::Element&           ele,
        const XFEM::InterfaceHandle&  ih
        )
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return DomainCoverageRatioPerNodeT<DRT::Element::hex8>(ele,ih);
    case DRT::Element::hex20:
      return DomainCoverageRatioPerNodeT<DRT::Element::hex20>(ele,ih);
    case DRT::Element::hex27:
      return DomainCoverageRatioPerNodeT<DRT::Element::hex27>(ele,ih);
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
  const Epetra_BLAS blas;
  
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
    
    const LINALG::SerialDenseMatrix& nodalpos_xi_domain(cell->NodalPosXiDomain());
    const int numnode_cell = cell->NumNode();
    
    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
      // coordinates of the current integration point in cell coordinates \eta^\boundary
      LINALG::Matrix<2,1> pos_eta_boundary;
      pos_eta_boundary(0) = intpoints.qxg[iquad][0];
      pos_eta_boundary(1) = intpoints.qxg[iquad][1];
      
      LINALG::Matrix<3,1> posXiDomain;
      mapEtaBToXiD(*cell, pos_eta_boundary, posXiDomain);
      
      // shape functions and their first derivatives
      LINALG::SerialDenseVector funct_boundary(DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
      DRT::UTILS::shape_function_2D(funct_boundary, pos_eta_boundary(0),pos_eta_boundary(1),cell->Shape());
      LINALG::SerialDenseMatrix deriv_boundary(3, DRT::UTILS::getNumberOfElementNodes(cell->Shape()));
      DRT::UTILS::shape_function_2D_deriv1(deriv_boundary, pos_eta_boundary(0),pos_eta_boundary(1),cell->Shape());
      
      // shape functions and their first derivatives
      static LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement,1> funct;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
      // get jacobian matrix d x / d \xi  (3x2)
      // dxyzdrs = xyze_boundary(i,k)*deriv_boundary(j,k);
      static LINALG::Matrix<3,2> dxyzdrs;
      blas.GEMM('N','T',3,2,numnode_cell,1.0,nodalpos_xi_domain.A(),nodalpos_xi_domain.LDA(),deriv_boundary.A(),deriv_boundary.LDA(),0.0,dxyzdrs.A(),dxyzdrs.M());
      
      // compute covariant metric tensor G for surface element (2x2)
      // metric = dxyzdrs(k,i)*dxyzdrs(k,j);
      static LINALG::Matrix<2,2> metric;
      metric.MultiplyTN(dxyzdrs,dxyzdrs);
      
      const double detmetric = sqrt(metric.Determinant());
      
      const double fac = intpoints.qwgt[iquad]*detmetric;//*detcell;
      if (fac < 0.0)
      {
        cout << endl;
        cout << "detmetric = " << detmetric << endl;
        cout << "fac       = " << fac << endl;
        dserror("negative fac! should be a bug!");
      }
      
      area_fict += fac;
      
    } // end loop over gauss points
  } // end loop over integration cells
  
  // scale result by area of one surface of the volume element
  return area_fict / base_area;
}

double XFEM::BoundaryCoverageRatio(
        const DRT::Element&           ele,
        const XFEM::InterfaceHandle&  ih
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

/*!
  Calculate ratio between fictitious element size and normal size
  */
template <DRT::Element::DiscretizationType DISTYPE>
vector<double> DomainIntCellCoverageRatioT(
        const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
  // number of nodes for element
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
 
  // information about domain integration cells
  const GEO::DomainIntCells&  domainIntCells(ih.GetDomainIntCells(ele.Id(),DISTYPE));
  
  double area_ele  = 0.0;
  
  std::vector<double> portions(domainIntCells.size(),0.0);
  
  // loop over integration cells
  int cellcount = 0;
  for (GEO::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
  {

    const LINALG::Matrix<3,1> cellcenter(cell->GetPhysicalCenterPosition(ele));
    
    DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
    switch (cell->Shape())
    {
      case DRT::Element::hex8:
      {
        gaussrule = DRT::UTILS::intrule_hex_8point;
        break;
      }
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      {
        gaussrule = DRT::UTILS::intrule_hex_27point;
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
      static LINALG::Matrix<3,1> pos_eta_domain;
      pos_eta_domain(0) = intpoints.qxg[iquad][0];
      pos_eta_domain(1) = intpoints.qxg[iquad][1];
      pos_eta_domain(2) = intpoints.qxg[iquad][2];

      // coordinates of the current integration point in element coordinates \xi
      static LINALG::Matrix<3,1> posXiDomain;
      GEO::mapEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain, posXiDomain);
      const double detcell = GEO::detEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain);
      
      // shape functions and their first derivatives
      static LINALG::Matrix<numnode,1> funct;
      DRT::UTILS::shape_function_3D(funct,posXiDomain(0),posXiDomain(1),posXiDomain(2),DISTYPE);
      
      const double fac = intpoints.qwgt[iquad]*detcell;

      area_ele += fac;
      portions[cellcount] += fac;
        
    } // end loop over gauss points
    cellcount++;
  } // end loop over integration cells

  for (unsigned icell = 0;icell < domainIntCells.size();++icell)
  {
    portions[icell] /= area_ele;
  }    
  
  return portions;
}

std::vector<double> XFEM::DomainIntCellCoverageRatio(
        const DRT::Element&           ele,
        const XFEM::InterfaceHandle&  ih
        )
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return DomainIntCellCoverageRatioT<DRT::Element::hex8>(ele,ih);
    case DRT::Element::hex20:
      return DomainIntCellCoverageRatioT<DRT::Element::hex20>(ele,ih);
    case DRT::Element::hex27:
      return DomainIntCellCoverageRatioT<DRT::Element::hex27>(ele,ih);
    default:
      dserror("add you distype here...");
      exit(1);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::AssemblyType XFEM::CheckForStandardEnrichmentsOnly(
    const ElementDofManager&   eleDofManager,
    const int                  numnode,
    const int*                 nodeids
)
{
  // find out whether we can use standard assembly or need xfem assembly
  XFEM::AssemblyType assembly_type = XFEM::standard_assembly;
  for (int inode = 0; inode < numnode; ++inode)
  {
    if (assembly_type == XFEM::xfem_assembly)
    {
      break;
    }
    const int gid = nodeids[inode];
    const std::set<XFEM::FieldEnr>& fields = eleDofManager.FieldEnrSetPerNode(gid);
    if (fields.size() != 4)
    {
      assembly_type = XFEM::xfem_assembly;
      break;
    };
    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fields.begin(); fieldenr != fields.end(); ++fieldenr)
    {
      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
        assembly_type = XFEM::xfem_assembly;
        break;
      };
    };
  };
  const int eledof = eleDofManager.NumElemDof();
  if (eledof != 0)
  {
    assembly_type = XFEM::xfem_assembly;
  }
  
  return assembly_type;
}

#endif  // #ifdef CCADISCRET
