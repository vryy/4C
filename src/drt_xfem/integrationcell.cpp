/*!
\file integrationcell.cpp

\brief integration cell classes for domain and boundary integration

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "integrationcell.H"
#include <string>
#include <sstream>
#include <blitz/array.h>
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


//! little helper function
static BlitzMat ConvertPosArrayToBlitz(
        const vector<vector<double> >&         pos_array,
        const DRT::Element::DiscretizationType distype,
        const int                              dim
        )
{
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    BlitzMat pos_array_blitz(dim,numnode,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<numnode; ++inode)
    {
        for (int isd=0; isd<dim; ++isd)
        {
            pos_array_blitz(isd,inode) = pos_array[inode][isd];
        }
    }    
    return pos_array_blitz;
}


//
//  ctor
//
XFEM::IntCell::IntCell(
        const DRT::Element::DiscretizationType distype) :
            distype_(distype)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
XFEM::IntCell::IntCell(
        const IntCell& old) : 
            distype_(old.distype_)
{
    return;   
}
 
//
// Print method
//
std::string XFEM::IntCell::toString() const
{
  return "";
}


//
//  ctor
//
XFEM::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const vector< vector<double> >& domainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_blitz_(ConvertPosArrayToBlitz(domainCoordinates, distype, 3))
{
    return;
}


//
//  ctor
//
XFEM::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const BlitzMat&                        domainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_blitz_(domainCoordinates)
{
    return;
}

        
//
//  ctor for dummy cells
//
XFEM::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype) :
            IntCell(distype),
            nodalpos_xi_domain_blitz_(GetDefaultCoordinates(distype))
{
    return;
}
        
/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
XFEM::DomainIntCell::DomainIntCell(
        const DomainIntCell& old) :
            IntCell(old),
            nodalpos_xi_domain_blitz_(old.nodalpos_xi_domain_blitz_)
{
    return;   
}
     
std::string XFEM::DomainIntCell::toString() const
{
    std::stringstream s;
    s << "DomainIntCell" << endl;
    s << nodalpos_xi_domain_blitz_ << endl;
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

// set element nodal coordinates according to given distype
BlitzMat XFEM::DomainIntCell::GetDefaultCoordinates(
        const DRT::Element::DiscretizationType distype) const
{
    const int nsd = 3;
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    BlitzMat coords(3,numnode);
    
    switch (distype)
    {
        case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        {
            for(int inode = 0; inode < numnode; ++inode)
            {
                for(int k = 0; k < nsd; ++k)
                {
                    coords(k,inode) = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[inode][k];
                }
            }
            break;
        }
        case DRT::Element::tet4: case DRT::Element::tet10:
        {
            for(int inode = 0; inode < numnode; ++inode)
            {
                for(int k = 0; k < nsd; ++k)
                {
                    coords(k,inode) = DRT::UTILS::eleNodeNumbering_tet10_nodes_reference[inode][k];
                }
            }
            break;
        }
        default:
            dserror("not supported in integrationcells. can be coded easily... ;-)");
    }
    return coords;
}

//
// return the center of the cell in physical coordinates
//
BlitzVec3 XFEM::DomainIntCell::GetPhysicalCenterPosition(const DRT::Element& ele) const
{
    // number of space dimensions
    //const int nsd = 3;
    
    // physical positions of cell nodes
    BlitzMat physcoord(3,27);
    this->NodalPosXYZ(ele, physcoord);
    
    // center in local coordinates
    BlitzVec3 localcenterpos;
    localcenterpos = DRT::UTILS::getLocalCenterPosition<3>(this->Shape());

    // shape functions
    BlitzVec funct(DRT::UTILS::getNumberOfElementNodes(this->Shape()));
    DRT::UTILS::shape_function_3D(funct,
            localcenterpos(0),
            localcenterpos(1),
            localcenterpos(2),
            this->Shape());
    
    //interpolate position to x-space
    BlitzVec3 x_interpol;
    for (int isd = 0; isd < 3; ++isd)
    {
        x_interpol(isd) = 0.0;
        for (int inode = 0; inode < DRT::UTILS::getNumberOfElementNodes(this->Shape()); ++inode)
        {
            x_interpol(isd) += funct(inode)*physcoord(isd,inode);
        }
    }
    
    return x_interpol;
}


//
//  ctor
//
XFEM::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType    distype,
        const int                                 surface_ele_gid,
        const vector< vector<double> >&           domainCoordinates,
        const vector< vector<double> >&           boundaryCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_blitz_(ConvertPosArrayToBlitz(domainCoordinates, distype, 3)),
            nodalpos_xi_boundary_blitz_(ConvertPosArrayToBlitz(boundaryCoordinates, distype, 2))
{
//    if (abs(blitz::sum(nodalpos_xi_boundary_blitz_)) < 1.0e-7 )
//    {
//      for(int ii=0; ii < boundaryCoordinates.size(); ii++)
//          for(int jj=0; jj < 3; jj++)
//                  printf("boundary = %f\n",boundaryCoordinates[ii][jj]);
//      cout << "Surface Ele Id: " << surface_ele_gid << endl;
//      cout << DRT::DistypeToString(distype) << endl;
//      dserror("something went wrong! A");
//    }
    return;
}

//
//  ctor
//
XFEM::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType    distype,
        const int                                 surface_ele_gid,
        const BlitzMat&                           domainCoordinates,
        const BlitzMat&                           boundaryCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_blitz_(domainCoordinates),
            nodalpos_xi_boundary_blitz_(boundaryCoordinates)
{
//    if (abs(blitz::sum(nodalpos_xi_boundary_blitz_)) < 1.0e-7 )
//    {
//      dserror("something went wrong! B");
//    }
    return;
}
        
        
/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
XFEM::BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            surface_ele_gid_(old.surface_ele_gid_),
            nodalpos_xi_domain_blitz_(old.nodalpos_xi_domain_blitz_),
            nodalpos_xi_boundary_blitz_(old.nodalpos_xi_boundary_blitz_)
{
    return;   
}
     
std::string XFEM::BoundaryIntCell::toString() const
{
    std::stringstream s;
    s << "BoundaryIntCell" << endl;
    s << nodalpos_xi_domain_blitz_ << endl;
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

//
// return the center of the cell in physical coordinates
//
BlitzVec3 XFEM::BoundaryIntCell::GetPhysicalCenterPosition(const DRT::Element& ele) const
{
    // number of space dimensions
    //const int nsd = 3;
    
    // physical positions of cell nodes
    BlitzMat physcoord(3,27);
    this->NodalPosXYZ(ele, physcoord);
    
    // center in local coordinates
    BlitzVec2 localcenterpos;
    localcenterpos = DRT::UTILS::getLocalCenterPosition<2>(this->Shape());

    // shape functions
    BlitzVec funct(DRT::UTILS::getNumberOfElementNodes(this->Shape()));
    DRT::UTILS::shape_function_2D(funct,
            localcenterpos(0),
            localcenterpos(1),
            this->Shape());
    
    //interpolate position to x-space
    BlitzVec3 x_interpol;
    for (int isd = 0; isd < 3; ++isd)
    {
        x_interpol(isd) = 0.0;
        for (int inode = 0; inode < DRT::UTILS::getNumberOfElementNodes(this->Shape()); ++inode)
        {
            x_interpol(isd) += funct(inode)*physcoord(isd,inode);
        }
    }
    
    return x_interpol;
}

#endif  // #ifdef CCADISCRET
