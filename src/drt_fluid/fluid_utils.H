/*!----------------------------------------------------------------------
\file fluid_utils.H
\brief utility functions for fluid problems

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUID_UTILS_H
#define FLUID_UTILS_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "../drt_lib/drt_discret.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_lib/linalg_mapextractor.H"

/*!
 * \brief namespace for fluid related stuff
 */
namespace FLUID_UTILS
{

/// Split dof row map of discretization into velocity and pressure maps
/*!

  We assume that each node possesses ndim velocity dofs and (optionally)
  one pressure dof. The dof row map is thus split in two.

  The velocity dofs are assigned to map 0 (the other map) and the pressure
  dofs are assigned to map 1 (the condition map).

  \param dis : (in) Discretization
  \param ndim : (in) dimensions of problem (2 or 3)
  \param extractor : (out) ready made map splitter

  \author u.kue
  \date 02/08
 */
void SetupFluidSplit(const DRT::Discretization& dis,
                     int ndim,
                     LINALG::MapExtractor& extractor);

/// Split dof row map of discretization into velocity and pressure maps
/*!

  We assume that each node possesses ndim velocity dofs and (optionally)
  one pressure dof. The dof row map is thus split in two.

  The velocity dofs are assigned to map 0 (the other map) and the pressure
  dofs are assigned to map 1 (the condition map).

  \author a.ger
  \date 03/08
 */
void SetupXFluidSplit(
        const DRT::Discretization&   dis,       ///< fluid discretization
        const RCP<XFEM::DofManager>  dofman,    ///< information about the dofs per node
        LINALG::MapExtractor&        extractor  ///< ready made map splitter
        );

}

#endif  // #ifndef FLUID_UTILS_H
#endif  // #ifdef CCADISCRET
