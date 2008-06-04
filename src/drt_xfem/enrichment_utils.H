/*!
\file enrichment_utils.H

\brief describes the enrichment types and classes

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef ENRICHMENT_UTILS_H
#define ENRICHMENT_UTILS_H

#include <blitz/array.h>
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "xfem.H"
#include "physics.H"
#include "dof_management.H"
#include "interface.H"
#include "coordinate_transformation.H"


namespace XFEM
{


std::map<XFEM::Enrichment, double> computeEnrvalMap(
        const RCP<XFEM::InterfaceHandle>      ih,
        const std::set<XFEM::Enrichment>&     enrset,
        const BlitzVec3&                      actpos,
        const XFEM::Enrichment::ApproachFrom  approachdirection
        );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedNodalShapefunction(
        const DRT::Element&                     ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager&    dofman,
        const XFEM::PHYSICS::Field        field,
        const BlitzVec3&                  actpos,      ///< Position in physical coordinates
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const BlitzVec&                   funct,
        BlitzVec&                         enr_funct
        );

void ComputeEnrichedNodalShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec3& actpos,                  ///< Position in physical coordinates
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy
        );

void ComputeEnrichedNodalShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec3& actpos,                   ///< Position in physical coordinates
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        const BlitzMat& derxy2,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy,
        BlitzMat& enr_derxy2
        );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedElementShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec3& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const BlitzVec& funct,
        BlitzVec& enr_funct
        );

void ComputeEnrichedElementShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec3& actpos,
        const XFEM::Enrichment::ApproachFrom approachdirection,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy
        );


void computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzVec& elementvalues,
  BlitzVec&     cellvalues
  );

void computeScalarCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzVec& elementvalues,
  BlitzVec&     cellvalues
  );

void computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&     cellvalues
  );

void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&     cellvalues
  );

void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::BoundaryIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&     cellvalues
  );


/*!
  Calculate matrix and rhs for stationary problem formulation
  */
template <DRT::Element::DiscretizationType DISTYPE>
double AreaRatioT(
        const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
    
    // number of nodes for element
    const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
    
    // dimension for 3d fluid element
    const int nsd = 3;
    
    // get node coordinates of the current element
    static blitz::TinyMatrix<double,nsd,numnode> xyze;
    DRT::UTILS::fillPositionArray<DISTYPE>(&ele, xyze);
    
    //double 
    double area_ele  = 0.0;
    double area_fict = 0.0;
    
    // information about domain integration cells
    const XFEM::DomainIntCells&  domainIntCells(ih.GetDomainIntCells(ele.Id(),DISTYPE));
    // loop over integration cells
    for (XFEM::DomainIntCells::const_iterator cell = domainIntCells.begin(); cell != domainIntCells.end(); ++cell)
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
            static XFEM::PosEtaDomain pos_eta_domain;
            pos_eta_domain(0) = intpoints.qxg[iquad][0];
            pos_eta_domain(1) = intpoints.qxg[iquad][1];
            pos_eta_domain(2) = intpoints.qxg[iquad][2];

            // coordinates of the current integration point in element coordinates \xi
            static XFEM::PosXiDomain posXiDomain;
            XFEM::mapEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain, posXiDomain);
            const double detcell = XFEM::detEtaToXi3D<XFEM::xfem_assembly>(*cell, pos_eta_domain);
            
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
            XFEM::Inverse3x3(xjm, det, xji);

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
                        
            map<int,bool> posInCondition;
            PositionWithinCondition(cellcenter, ih, posInCondition);
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

    return area_fict / area_ele;
}

inline double AreaRatio(
        const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
        const XFEM::InterfaceHandle&  ih             ///< connection to the interface handler
        )
{
  switch (ele.Shape())
  {
    case DRT::Element::hex8:
      return XFEM::AreaRatioT<DRT::Element::hex8>(ele,ih);
    case DRT::Element::hex27:
      return XFEM::AreaRatioT<DRT::Element::hex27>(ele,ih);
    default:
      dserror("add you distype here...");
      exit(1);
  }
}

}


#endif  // #ifndef ENRICHMENT_UTILS_H
#endif  // #ifdef CCADISCRET
