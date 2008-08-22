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
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


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


//
//  ctor
//
GEO::IntCell::IntCell(
        const DRT::Element::DiscretizationType distype) :
            distype_(distype)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
GEO::IntCell::IntCell(
        const IntCell& old) : 
            distype_(old.distype_)
{
    return;   
}
 
//
// Print method
//
std::string GEO::IntCell::toString() const
{
  return "";
}


//
//  ctor
//
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const vector< vector<double> >& domainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_blitz_(ConvertPosArrayToBlitz<3>(domainCoordinates, distype))
{
    return;
}


//
//  ctor
//
GEO::DomainIntCell::DomainIntCell(
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
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype) :
            IntCell(distype),
            nodalpos_xi_domain_blitz_(GetDefaultCoordinates(distype))
{
    return;
}
        
/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DomainIntCell& old) :
            IntCell(old),
            nodalpos_xi_domain_blitz_(old.nodalpos_xi_domain_blitz_)
{
    return;   
}
     
std::string GEO::DomainIntCell::toString() const
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void GEO::DomainIntCell::NodalPosXYZ(const DRT::Element& ele, BlitzMat& xyz_cell) const
{
    ComputePhysicalCoordinates(ele, (*this), xyz_cell);
    return;
}

// set element nodal coordinates according to given distype
BlitzMat GEO::DomainIntCell::GetDefaultCoordinates(
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
BlitzVec3 GEO::DomainIntCell::GetPhysicalCenterPosition(const DRT::Element& ele) const
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
GEO::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType    distype,
        const int                                 surface_ele_gid,
        const vector< vector<double> >&           domainCoordinates,
        const vector< vector<double> >&           boundaryCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_blitz_(ConvertPosArrayToBlitz<3>(domainCoordinates, distype)),
            nodalpos_xi_boundary_blitz_(ConvertPosArrayToBlitz<2>(boundaryCoordinates, distype))
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
GEO::BoundaryIntCell::BoundaryIntCell(
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
GEO::BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            surface_ele_gid_(old.surface_ele_gid_),
            nodalpos_xi_domain_blitz_(old.nodalpos_xi_domain_blitz_),
            nodalpos_xi_boundary_blitz_(old.nodalpos_xi_boundary_blitz_)
{
    return;   
}
     
std::string GEO::BoundaryIntCell::toString() const
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void GEO::BoundaryIntCell::NodalPosXYZ(const DRT::Element& ele, BlitzMat& xyz_cell) const
{
    ComputePhysicalCoordinates(ele, (*this), xyz_cell);
    return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BlitzVec3 GEO::BoundaryIntCell::GetPhysicalCenterPosition(const DRT::Element& ele) const
{
    // number of space dimensions
    //const int nsd = 3;
    
    // physical positions of cell nodes
    static BlitzMat physcoord(3,27);
    this->NodalPosXYZ(ele, physcoord);
    
    // center in local coordinates
    const BlitzVec2 localcenterpos(DRT::UTILS::getLocalCenterPosition<2>(this->Shape()));

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
