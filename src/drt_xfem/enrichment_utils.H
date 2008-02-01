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

/*!
 * \brief For a given situation compute the enriched shape functions
 * 
 * simplest case: jump or void enrichment
 * -> no chain rule, since enrichment function derivative is zero
 * TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
 */
void ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const BlitzVec& funct,
        BlitzVec& enr_funct
        );

void ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        const BlitzMat& derxy2,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy,
        BlitzMat& enr_derxy2
        );

void ComputeEnrichedStressShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const BlitzVec& actpos,
        const BlitzVec& funct,
        const BlitzMat& derxy,
        BlitzVec& enr_funct,
        BlitzMat& enr_derxy
        );

}


#endif  // #ifndef ENRICHMENT_UTILS_H
#endif  // #ifdef CCADISCRET
