/*!
\file enrichment_utils.cpp

\brief describes the enrichment types and classes

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

\warning this combustion module related file will be deleted within the next time!!!
*/


#include "enrichment_utils.H"
#include "../drt_combust/combust_defines.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../linalg/linalg_serialdensevector.H"


/*----------------------------------------------------------------------*
 | constructor called for xfsi problems (void enrichment)            ag |
 *----------------------------------------------------------------------*/
XFEM::ElementEnrichmentValues::ElementEnrichmentValues(
        const DRT::Element&                    ele,
        COMBUST::InterfaceHandleCombust*       ih,                ///< interface information
        const XFEM::ElementDofManager&         dofman,
        const LINALG::Matrix<3,1>&             actpos,
        const bool                             boundary_integral,
        const int                              boundary_label
        ) :
          ele_(ele),
          dofman_(dofman)
{
    enrvals_.clear();
    enrvalnodes_.clear();  // not used for void enrichment
    enrvalderxy_.clear();  // not used for void enrichment
    enrvalderxy2_.clear(); // not used for void enrichment

    const std::set<XFEM::Enrichment>& enrset(dofman.getUniqueEnrichments());
    for (std::set<XFEM::Enrichment>::const_iterator enriter = enrset.begin(); enriter != enrset.end(); ++enriter)
    {
        XFEM::Enrichment::ApproachFrom   approachdirection = XFEM::Enrichment::approachUnknown;
        if (boundary_integral and enriter->XFEMConditionLabel() == boundary_label)
          approachdirection = XFEM::Enrichment::approachFromPlus;
        else
          approachdirection = XFEM::Enrichment::approachUnknown;

        const double enrval = enriter->EnrValue(actpos, *ih, approachdirection);
        enrvals_[*enriter] = enrval;
    }
    return;
}


/*----------------------------------------------------------------------*
 | interpolate field from element node values to cell node values based |
 | on a level-set field                                     henke 10/09 |
 | remark: function used for modified jump enrichment strategy          |
 *----------------------------------------------------------------------*/
void XFEM::InterpolateCellValuesFromElementValuesLevelSet(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues)
{
  if (ele.Shape() != DRT::Element::hex8)
    dserror("OutputToGmsh() only available for hex8 elements! However, this is easy to extend.");
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  const size_t numparam  = dofman.NumDofPerField(field);

  // copy element phi vector from std::vector (phivalues) to LINALG::Matrix (phi)
  LINALG::Matrix<numnode,1> phi;
  for (size_t inode=0; inode<numnode; ++inode)
    phi(inode) = phivalues[inode];

  // get coordinates of cell vertices
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.CellNodalPosXiDomain());

  // compute enrichment values based on a level set field 'phi'
#ifndef COMBUST_SXFEM
  const XFEM::ElementEnrichmentValues enrvals(ele, dofman, cell, phi);
  LINALG::SerialDenseVector enr_funct(numparam);
  LINALG::SerialDenseVector funct(numnode);
#else
  LINALG::SerialDenseVector enr_funct(numparam);
  LINALG::Matrix<8,1> funct(true);
  LINALG::Matrix<3,8> derxy(true);
  LINALG::Matrix<6,8> derxy2(true);
#endif
  cellvalues.Zero();
  for (int inode = 0; inode < cell.NumNode(); ++inode)
  {
    // evaluate shape functions
    DRT::UTILS::shape_function_3D(funct,nodalPosXiDomain(0,inode),nodalPosXiDomain(1,inode),nodalPosXiDomain(2,inode),DRT::Element::hex8);

#ifdef COMBUST_SXFEM
    const XFEM::ElementEnrichmentValues enrvals(ele,dofman,phi,cell,funct,derxy,derxy2);
#endif
    // evaluate enriched shape functions
    enrvals.ComputeModifiedEnrichedNodalShapefunction(field, funct, enr_funct);

    switch (field)
    {
    // scalar fields
    case XFEM::PHYSICS::Pres:
    {
    // interpolate value
    for (size_t iparam = 0; iparam < numparam; ++iparam)
      cellvalues(0,inode) += elementvalues(0,iparam) * enr_funct(iparam);
    break;
    }
    // vector fields
    case XFEM::PHYSICS::Velx:
    {
      for (std::size_t iparam = 0; iparam < numparam; ++iparam)
        for (std::size_t isd = 0; isd < 3; ++isd)
          cellvalues(isd,inode) += elementvalues(isd,iparam) * enr_funct(iparam);
      break;
    }
    default:
      dserror("interpolation to cells not available for this field");
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | interpolate field from element node values to cell node values based |
 | on a level-set field                                     henke 01/11 |
 | remark: function used for modified jump normal enrichment strategy   |
 *----------------------------------------------------------------------*/
void XFEM::InterpolateCellValuesFromElementValuesLevelSetNormal(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const LINALG::Matrix<3,8>&            gradphi,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues)
{
  if (ele.Shape() != DRT::Element::hex8)
    dserror("OutputToGmsh() only available for hex8 elements! However, this is easy to extend.");
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  const size_t numparam  = dofman.NumDofPerField(field);

  // copy element phi vector from std::vector (phivalues) to LINALG::Matrix (phi)
  LINALG::Matrix<numnode,1> phi;
  for (size_t inode=0; inode<numnode; ++inode)
    phi(inode) = phivalues[inode];

  // get coordinates of cell vertices
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.CellNodalPosXiDomain());

  // compute enrichment values based on a level set field 'phi'
  const XFEM::ElementEnrichmentValues enrvals(ele, dofman, cell, phi);

  LINALG::SerialDenseVector enr_funct(numparam,true);
  LINALG::SerialDenseVector funct(numnode,true);

  cellvalues.Zero();
  for (int ivertex = 0; ivertex < cell.NumNode(); ++ivertex)
  {
    // evaluate shape functions
    DRT::UTILS::shape_function_3D(funct,nodalPosXiDomain(0,ivertex),nodalPosXiDomain(1,ivertex),nodalPosXiDomain(2,ivertex),DRT::Element::hex8);

    XFEM::ApproxFuncNormalVector<0,8> shp(true);
    // fill approximation functions for XFEM
    for (size_t iparam = 0; iparam != numparam; ++iparam)
    {
      shp.velx.d0.s(iparam) = funct(iparam);
      shp.vely.d0.s(iparam) = funct(iparam);
      shp.velz.d0.s(iparam) = funct(iparam);
    }

    LINALG::Matrix<3,1> normal(true);
#ifdef COLLAPSE_FLAME
    // get coordinates of cell vertices
    const LINALG::SerialDenseMatrix& nodalPosXYZ(cell.CellNodalPosXYZ());
    normal(0) = nodalPosXYZ(0,ivertex);
    normal(1) = nodalPosXYZ(1,ivertex);
    normal(2) = 0.0;
    const double norm = normal.Norm2(); // sqrt(normal(0)*normal(0) + normal(1)*normal(1) + normal(2)*normal(2))
    if (norm == 0.0) dserror("norm of normal vector is zero!");
    normal.Scale(1.0/norm);
#endif

    // shape functions and derivatives for nodal parameters (dofs)
    enrvals.ComputeNormalShapeFunction(funct,gradphi,
//#ifdef COLLAPSE_FLAME
//        normal,
//#endif
        shp);

    switch (field)
    {
    // vector fields
    case XFEM::PHYSICS::Velx:
    {
      const int* nodeids = ele.NodeIds();

      std::size_t velncounter = 0;
      for (std::size_t inode=0; inode<numnode; ++inode)
      {
        // standard shape functions are identical for all vector components
        // shp.velx.d0.s == shp.vely.d0.s == shp.velz.d0.s
        cellvalues(0,ivertex) += elementvalues(0,inode)*shp.velx.d0.s(inode);
        cellvalues(1,ivertex) += elementvalues(1,inode)*shp.vely.d0.s(inode);
        cellvalues(2,ivertex) += elementvalues(2,inode)*shp.velz.d0.s(inode);

        const int gid = nodeids[inode];
        const std::set<XFEM::FieldEnr>& enrfieldset = dofman.FieldEnrSetPerNode(gid);

        for (std::set<XFEM::FieldEnr>::const_iterator enrfield =
            enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
        {
          if (enrfield->getField() == XFEM::PHYSICS::Veln)
          {
            cellvalues(0,ivertex) += elementvalues(3,velncounter)*shp.velx.d0.n(velncounter);
            cellvalues(1,ivertex) += elementvalues(3,velncounter)*shp.vely.d0.n(velncounter);
            cellvalues(2,ivertex) += elementvalues(3,velncounter)*shp.velz.d0.n(velncounter);

            velncounter += 1;
          }
        }
      }
      if (velncounter != dofman.NumDofPerField(XFEM::PHYSICS::Veln)) dserror("Alles falsch, du Depp!");
      dsassert(velncounter == dofman.NumDofPerField(XFEM::PHYSICS::Veln), "mismatch in information from eledofmanager!");

      break;
    }
    default:
      dserror("interpolation to cells not available for this field");
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | interpolate field from element node values to cell node values based |
 | on a level-set field                                     henke 10/09 |
 | remark: function used for modified kink enrichment strategy          |
 *----------------------------------------------------------------------*/
void XFEM::InterpolateCellValuesFromElementValuesLevelSetKink(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues)
{
  if (ele.Shape() != DRT::Element::hex8)
    dserror("OutputToGmsh() only available for hex8 elements! However, this is easy to extend.");
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  const size_t numparam  = dofman.NumDofPerField(field);

  // copy element phi vector from std::vector (phivalues) to LINALG::Matrix (phi)
//  LINALG::Matrix<numnode,1> phi;
  LINALG::SerialDenseVector phi(numnode);
  for (size_t inode=0; inode<numnode; ++inode)
    phi(inode) = phivalues[inode];

  // get coordinates of cell vertices
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.CellNodalPosXiDomain());

  LINALG::SerialDenseVector enr_funct(numparam);
  enr_funct.Zero();
  LINALG::SerialDenseVector funct(numnode);

  cellvalues.Zero();
  for (int inode = 0; inode < cell.NumNode(); ++inode)
  {
    // evaluate shape functions
    DRT::UTILS::shape_function_3D(funct,nodalPosXiDomain(0,inode),nodalPosXiDomain(1,inode),nodalPosXiDomain(2,inode),DRT::Element::hex8);
    // first and second derivatives are dummy matrices needed to call XFEM::ElementEnrichmentValues for kink enrichment
    const LINALG::Matrix<3,numnode> derxy(true);
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::hex8>::numderiv2, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement> derxy2(true);

    // compute enrichment values based on a level set field 'phi'
    const XFEM::ElementEnrichmentValues enrvals(ele, dofman, phi, funct, derxy, derxy2);

    // evaluate enriched shape functions
    // remark: since we do not compute derivatives of enriched shape functions (chain rule for kink
    //         enrichment!), we can use the same function for all types of enrichments.
    enrvals.ComputeModifiedEnrichedNodalShapefunction(field, funct, enr_funct);

    switch (field)
    {
    // scalar fields
    case XFEM::PHYSICS::Pres:
    {
    // interpolate value
    for (size_t iparam = 0; iparam < numparam; ++iparam)
      cellvalues(0,inode) += elementvalues(0,iparam) * enr_funct(iparam);
    break;
    }
    // vector fields
    case XFEM::PHYSICS::Velx:
    {
      for (std::size_t iparam = 0; iparam < numparam; ++iparam)
        for (std::size_t isd = 0; isd < 3; ++isd)
          cellvalues(isd,inode) += elementvalues(isd,iparam) * enr_funct(iparam);
      break;
    }
    default:
      dserror("interpolation to cells not available for this field");
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | interpolate field from element node values to cell node values based |
 | on a level-set field                                 schott 05/17/10 |
 | remark: function used for 2-phase-flow                               |
 | with jump enrichments in pressure and kink enrichment in velocity    |
 *----------------------------------------------------------------------*/
void XFEM::InterpolateCellValuesFromElementValuesLevelSetKinkJump(
  const DRT::Element&                   ele,
  const XFEM::ElementDofManager&        dofman,
  const GEO::DomainIntCell&             cell,
  const std::vector<double>&            phivalues,
  const XFEM::PHYSICS::Field            field,
  const LINALG::SerialDenseMatrix&      elementvalues,
  LINALG::SerialDenseMatrix&            cellvalues)
{
  if (ele.Shape() != DRT::Element::hex8)
    dserror("OutputToGmsh() only available for hex8 elements! However, this is easy to extend.");
  const size_t numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
  const size_t numparam  = dofman.NumDofPerField(field);

  // copy element phi vector from std::vector (phivalues) to LINALG::Matrix (phi)
//  LINALG::Matrix<numnode,1> phi;
  LINALG::SerialDenseVector phi(numnode);
  for (size_t inode=0; inode<numnode; ++inode)
    phi(inode) = phivalues[inode];

  // get coordinates of cell vertices
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.CellNodalPosXiDomain());

  LINALG::SerialDenseVector funct(numnode);

  cellvalues.Zero();

  for (int inode = 0; inode < cell.NumNode(); ++inode)
  {
    // evaluate shape functions
    DRT::UTILS::shape_function_3D(funct,nodalPosXiDomain(0,inode),nodalPosXiDomain(1,inode),nodalPosXiDomain(2,inode),DRT::Element::hex8);
    // first and second derivatives are dummy matrices needed to call XFEM::ElementEnrichmentValues for kink enrichment
    const LINALG::Matrix<3,numnode> derxy(true);
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumDeriv2<DRT::Element::hex8>::numderiv2, DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement> derxy2(true);

    const XFEM::ElementEnrichmentValues enrvals(ele, dofman, phi, cell, funct, derxy, derxy2);
    LINALG::SerialDenseVector enr_funct(numparam);
    enr_funct.Zero();

    enrvals.ComputeModifiedEnrichedNodalShapefunction(field, funct, enr_funct);

    switch (field)
    {
    // scalar fields
    case XFEM::PHYSICS::Pres:
    {
    // interpolate value
    for (size_t iparam = 0; iparam < numparam; ++iparam)
      cellvalues(0,inode) += elementvalues(0,iparam) * enr_funct(iparam);
    break;
    }
    // vector fields
    case XFEM::PHYSICS::Velx:
    {
      for (std::size_t iparam = 0; iparam < numparam; ++iparam)
        for (std::size_t isd = 0; isd < 3; ++isd)
            cellvalues(isd,inode) += elementvalues(isd,iparam) * enr_funct(iparam);
      break;
    }
    default:
      dserror("interpolation to cells not available for this field");
    }
  }
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFEM::computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&                 ele,
  COMBUST::InterfaceHandleCombust*    ih,
  const XFEM::ElementDofManager&      dofman,
  const GEO::DomainIntCell&           cell,
  const XFEM::PHYSICS::Field          field,
  const LINALG::SerialDenseMatrix&    elementvalues,
  LINALG::SerialDenseMatrix&          cellvalues)
{
  const LINALG::SerialDenseMatrix& nodalPosXiDomain(cell.CellNodalPosXiDomain());

  const XFEM::ElementEnrichmentValues enrvals(
        ele,
        ih,
        dofman,
        cell.GetPhysicalCenterPosition(),
        false, -1);

  cellvalues.Zero();
  for (int incn = 0; incn < cell.NumNode(); ++incn)
  {
    const int numparam  = dofman.NumDofPerField(field);
    if (numparam == 0)
      continue;

    const DRT::Element::DiscretizationType eleval_distype = dofman.getDisTypePerField(field);
    const int numvirtnode = DRT::UTILS::getNumberOfElementNodes(eleval_distype);
//    if (numvirtnode != numparam) dserror("bug");

    LINALG::SerialDenseVector funct(numvirtnode);
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
      for (int ientry = 0; ientry < 9; ++ientry)
        cellvalues(ientry,incn) += elementvalues(ientry,iparam) * enr_funct(iparam);
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::AssemblyType XFEM::ComputeAssemblyType(
    const ElementDofManager&   eleDofManager,
    const std::size_t          numnode,
    const int*                 nodeids)
{
  // find out whether we can use standard assembly or need xfem assembly
  XFEM::AssemblyType assembly_type = XFEM::standard_assembly;
  for (std::size_t inode = 0; inode < numnode; ++inode)
  {
    if (assembly_type == XFEM::xfem_assembly)
      break;

    const int gid = nodeids[inode];
    const std::set<XFEM::FieldEnr>& fields = eleDofManager.FieldEnrSetPerNode(gid);
    if (fields.size() != 4)
    {
      assembly_type = XFEM::xfem_assembly;
      break;
    }

    for (std::set<XFEM::FieldEnr>::const_iterator fieldenr = fields.begin(); fieldenr != fields.end(); ++fieldenr)
      if (fieldenr->getEnrichment().Type() != XFEM::Enrichment::typeStandard)
      {
        assembly_type = XFEM::xfem_assembly;
        break;
      };
  };

  if (eleDofManager.NumElemDof() != 0)
    assembly_type = XFEM::xfem_assembly;


  return assembly_type;
}
