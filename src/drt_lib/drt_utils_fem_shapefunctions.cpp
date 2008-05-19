/*!
 \file drt_utils_fem_shapefunctions.cpp

 \brief Provide a node numbering scheme together with a set of shape functions

 Provided are 1D, 2D and 3D shape functions

 The surface mapping gives the node numbers such that the 2D shape functions can be used
 Nodal mappings describe the relation between volume, surface and line node numbering.
 They should be used as the only reference for such relationships.
 The corresponding graphics and a detailed description can be found in the Baci guide in the Convention chapter.
 The numbering of lower order elements is included in the higher order element, such that
 e.g. the hex8 volume element uses only the first 8 nodes of the hex27 mapping

 !!!!
 The corresponding graphics and a detailed description can be found
 in the Baci guide in the Convention chapter.
 !!!!

 \author Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15236
 */
#ifdef CCADISCRET

//#include "drt_element.H"
//#include "drt_utils.H"
//#include "drt_dserror.H"


#endif  // #ifdef CCADISCRET
