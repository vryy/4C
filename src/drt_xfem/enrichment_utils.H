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
#include "xfem.H"
#include "physics.H"
#include "dof_management.H"
#include "interface.H"

using namespace Teuchos;


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


void computeScalarCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const blitz::Array<double,1>& elementvalues,
  blitz::Array<double,1>&     cellvalues
  );

void computeVectorCellNodeValues(
  const DRT::Element&  ele,
  const RCP<XFEM::InterfaceHandle>&  ih,
  const XFEM::ElementDofManager& dofman,
  const XFEM::DomainIntCell& cell,
  const XFEM::PHYSICS::Field field,
  const blitz::Array<double,2>& elementvalues,
  blitz::Array<double,2>&     cellvalues
  );

}


#endif  // #ifndef ENRICHMENT_UTILS_H
#endif  // #ifdef CCADISCRET
