/*----------------------------------------------------------------------*/
/*!
\file fsi_create_boundary.H

\brief Basis of all FSI algorithms

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef FSI_CREATE_BOUNDARY_H
#define FSI_CREATE_BOUNDARY_H

/*----------------------------------------------------------------------*/
// create ale discretization parallel to the fluid one
/*----------------------------------------------------------------------*/
void CreateBoundaryDiscretization(RCP<DRT::Discretization> cutterdis);

#endif
#endif
