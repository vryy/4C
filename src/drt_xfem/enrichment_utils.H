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

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_geometry/vector_definitions.H"
#include "physics.H"
#include "dof_management.H"
#include "../drt_geometry/coordinate_transformation.H"


class DRT::Element;
class DRT::Discretization;

namespace XFEM
{

  class InterfaceHandle;

  /*!
   * \brief for each given Enrichment in enrset, compute the enrichment value at actpos
   * 
   * if one already knows, an which side of the discontinuity one computes something,
   * the information can be given via approachdirection
   */
std::map<XFEM::Enrichment, double> computeEnrvalMap(
        const RCP<XFEM::InterfaceHandle>      ih,               ///< interface onformation
        const std::set<XFEM::Enrichment>&     enrset,           ///< enrichments for which one seeks enrichment values
        const BlitzVec3&                      actpos,           ///< a 3d position
        const XFEM::Enrichment::ApproachFrom  approachdirection ///< is location with respect to interface is known or not
        );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedNodalShapefunction(
    const DRT::Element&                        ele,
    const RCP<XFEM::InterfaceHandle>           ih,
    const XFEM::ElementDofManager&             dofman,
    const XFEM::PHYSICS::Field                 field,
    const std::map<XFEM::Enrichment, double>&  enrvals,
    const BlitzVec&                            funct,
    BlitzVec&                                  enr_funct
    );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedNodalShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy
        );

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedNodalShapefunction(
        const DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const std::map<XFEM::Enrichment, double>&  enrvals,
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
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec& funct,
        BlitzVec& enr_funct
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
        const std::map<XFEM::Enrichment, double>&  enrvals,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy
        );

//! for output, compute values on integration cells
void computeScalarCellNodeValuesFromNodalUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzVec& elementvalues,
  BlitzVec&     cellvalues
  );

//! for output, compute values on the integration cells
void computeScalarCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzVec& elementvalues,
  BlitzVec&     cellvalues
  );

//! for output, compute values on the integration cells
void computeTensorCellNodeValuesFromElementUnknowns(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&     cellvalues
  );

//! for output, compute values on the integration cells
void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&     cellvalues
  );

//! for output, compute values on the integration cells
void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const GEO::BoundaryIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const BlitzMat& elementvalues,
  BlitzMat&     cellvalues
  );

/*! 
 * \brief calculates the fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 */
double DomainCoverageRatio(
    const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
    const XFEM::InterfaceHandle&  ih             ///< interface information
    );

/*! 
 * \brief calculates the area fraction of an element that is covered by fictitious domains.
 * \warning should be done separately for each enrichment, but currently it's not done
 */
double BoundaryCoverageRatio(
    const DRT::Element&           ele,           ///< the element whose area ratio we want to compute
    const XFEM::InterfaceHandle&  ih             ///< interface information
    );

}


#endif  // #ifndef ENRICHMENT_UTILS_H
#endif  // #ifdef CCADISCRET
