/*!
\file spacetime_boundary.cpp

\brief integration cell classes for domain and boundary integration

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "spacetime_boundary.H"
#include <string>
#include <sstream>
#include <blitz/array.h>
#include "../drt_lib/drt_utils.H"
//#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


//! little helper function
template <int dim>
static BlitzMat ConvertPosArrayToBlitz(
        const vector<vector<double> >&         pos_array,
        const DRT::Element::DiscretizationType distype
        )
{
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    BlitzMat pos_array_blitz(dim,numnode);
    for (int inode=0; inode<numnode; ++inode)
    {
        for (int isd=0; isd<dim; ++isd)
        {
            pos_array_blitz(isd,inode) = pos_array[inode][isd];
        }
    }    
    return pos_array_blitz;
}


/*!
 * \brief create array with physical coordinates based an local coordinates of a parent element
 */
template<class Cell>
static void ComputePhysicalCoordinates(
        const DRT::Element&  ele,  ///< parent element
        const Cell&          cell, ///< integration cell whose coordinates we'd like to transform
        BlitzMat&            physicalCoordinates
        )
{
    const BlitzMat eleCoord(DRT::UTILS::InitialPositionArrayBlitz(&ele));
    //DRT::UTILS::fillInitialPositionArray(&ele, eleCoord);
    const BlitzMat* nodalPosXiDomain = cell.NodalPosXiDomainBlitz();
    
    const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // return value
    //BlitzMat physicalCoordinates(3, nen_cell);
    physicalCoordinates = 0.0;
    // for each cell node, compute physical position
    const int nen_ele = ele.NumNode();
    BlitzVec funct(nen_ele);
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        DRT::UTILS::shape_function_3D(funct,
                (*nodalPosXiDomain)(0, inen),
                (*nodalPosXiDomain)(1, inen),
                (*nodalPosXiDomain)(2, inen),
                ele.Shape());

        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < nen_ele; ++j)
            {
              physicalCoordinates(i,inen) += eleCoord(i, j) * funct(j);
            }
        }
    };
    return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell(
    const DRT::Element* boundaryele,
    const BlitzMat&     posnp,
    const BlitzMat&     posn
    ) :
      boundaryele_(boundaryele),
      posnp_(posnp),
      posn_(posn)
{
    return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
XFEM::SpaceTimeBoundaryCell::SpaceTimeBoundaryCell(
    const SpaceTimeBoundaryCell& old
    ) : 
      boundaryele_(old.boundaryele_),
      posnp_(old.posnp_),
      posn_(old.posn_)
{
    return;   
}
        
 
//
// Print method
//
std::string XFEM::SpaceTimeBoundaryCell::toString() const
{
    std::stringstream s;
    s << "SpaceTimeBoundaryCell" << endl;
//    MCONST_FOREACH(vector< vector<double> >, coordinate, nodalpos_xi_domain_)
//    {
//        s << "[";
//        MPFOREACH(vector<double>, val, coordinate)
//        {
//            s << *val << " ";
//        };
//        s << "]" << endl;
//    };
    return s.str();
}




     


#endif  // #ifdef CCADISCRET
