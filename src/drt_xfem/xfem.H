/*!
\file xfem.H

\brief general xfem routines without a class

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFEM_H
#define XFEM_H

#include <blitz/array.h>
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "physics.H"
#include "dof_management.H"
#include "interface.H"

using namespace Teuchos;


namespace XFEM
{

enum Assembly
{
    xfem_assembly,
    standard_assembly
};

enum Integration
{
    sub_integration,
    standard_integration
};

inline blitz::Array<double,1> toBlitzArray(const double* x);

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
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& funct,
        blitz::Array<double,1>& enr_funct
        );

void ComputeEnrichedShapefunction(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1>& actpos,
        const blitz::Array<double,1>& funct,
        const blitz::Array<double,2>& derxy,
        const blitz::Array<double,2>& derxy2,
        blitz::Array<double,1>& enr_funct,
        blitz::Array<double,2>& enr_derxy,
        blitz::Array<double,2>& enr_derxy2
        );

}


#endif  // #ifndef XFEM_H
#endif  // #ifdef CCADISCRET
