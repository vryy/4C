/*----------------------------------------------------------------------*/
/*!
\file drt_meshfree_utils.cpp

\brief service methods for a given meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#include "drt_meshfree_utils.H"
#include "drt_meshfree_node.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/standardtypes_cpp.H"

/*--------------------------------------------------------------------------*
 | reduces node position to dimension of face                     nis Jan14 |
 *--------------------------------------------------------------------------*/
std::vector<int> DRT::MESHFREE::ReduceDimensionOfFaceNodes(
  const LINALG::SerialDenseMatrix & xyz,
  LINALG::SerialDenseMatrix       & xyz_reduced
  )
{
  // initilize vector that holds dimensions in which face lies
  std::vector<int> dims;

  // find coordinates in which face lies
  int ndims = 0;
  const int nnode = xyz.ColDim(); // number of nodes
  const int dim   = xyz.RowDim(); // dimension of input nodes
  for(int idim=0; idim<dim; ++idim)
  {
    for(int inode=0; inode<(nnode-1); ++inode)
    {
      if (!(abs(xyz(idim,inode)-xyz(idim,inode+1))<EPS12))
      {
        if (ndims<(dim-1))
        {
          dims.push_back(idim);
          ndims++;
          break;
        }
        else
          dserror("So far, faces have to be aligned with a coordinate axis!");
      }
    }
  }

  // set size of output matrix
  xyz_reduced.LightShape(ndims,nnode);

  // copy node positions for reduced dimensions
  for(int idim=0; idim<ndims; ++idim)
    for(int inode=0; inode<nnode; ++inode)
      xyz_reduced(idim,inode) = xyz(dims[idim],inode);

  return dims;
}

// this is a secret place!
