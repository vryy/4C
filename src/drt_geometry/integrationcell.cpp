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
#include "intersection_service.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


//! translate between std::vector and blitz array
template <int dim>
static LINALG::SerialDenseMatrix ConvertPosArrayToLINALG(
        const std::vector<std::vector<double> >&   pos_array,
        const DRT::Element::DiscretizationType     distype
        )
{
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    LINALG::SerialDenseMatrix pos_array_blitz(dim,numnode);
    for (int inode=0; inode<numnode; ++inode)
    {
        for (int isd=0; isd<dim; ++isd)
        {
            pos_array_blitz(isd,inode) = pos_array[inode][isd];
        }
    }    
    return pos_array_blitz;
}


//! create array with physical coordinates based an local coordinates of a parent element
template<class Cell>
static void ComputePhysicalCoordinates(
        const DRT::Element&        ele,  ///< parent element
        const Cell&                cell, ///< integration cell whose coordinates we'd like to transform
        LINALG::SerialDenseMatrix& physicalCoordinates
        )
{
    const BlitzMat eleCoord(GEO::InitialPositionArrayBlitz(&ele));
    //DRT::UTILS::fillInitialPositionArray(&ele, eleCoord);
    const LINALG::SerialDenseMatrix* nodalPosXiDomain = cell.NodalPosXiDomain();
    
    const int nen_cell = DRT::UTILS::getNumberOfElementNodes(cell.Shape());
    
    // return value
    //BlitzMat physicalCoordinates(3, nen_cell);
    physicalCoordinates.Zero();
    // for each cell node, compute physical position
    const int nen_ele = ele.NumNode();
    LINALG::SerialDenseVector funct(nen_ele);
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        DRT::UTILS::shape_function_3D(funct,
                (*nodalPosXiDomain)(0, inen),
                (*nodalPosXiDomain)(1, inen),
                (*nodalPosXiDomain)(2, inen),
                ele.Shape());

        for (int j = 0; j < nen_ele; ++j)
        {
            for (int i = 0; i < 3; ++i)
            {
              physicalCoordinates(i,inen) += eleCoord(i, j) * funct(j);
            }
        }
    };
    return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::IntCell::IntCell(
        const DRT::Element::DiscretizationType distype) :
            distype_(distype)
{
    return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::IntCell::IntCell(
        const IntCell& old) : 
            distype_(old.distype_)
{
    return;   
}
 
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string GEO::IntCell::toString() const
{
  return "";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const std::vector< std::vector<double> >& domainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_(ConvertPosArrayToLINALG<3>(domainCoordinates, distype))
{
    return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const LINALG::SerialDenseMatrix&       domainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_(domainCoordinates)
{
    return;
}

        
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype) :
            IntCell(distype),
            nodalpos_xi_domain_(GetDefaultCoordinates(distype))
{
    return;
}
        
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::DomainIntCell::DomainIntCell(
        const DomainIntCell& old) :
            IntCell(old),
            nodalpos_xi_domain_(old.nodalpos_xi_domain_)
{
    return;   
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::string GEO::DomainIntCell::toString() const
{
    std::stringstream s;
    s << "DomainIntCell" << endl;
    s << nodalpos_xi_domain_ << endl;
//    MCONST_FOREACH(std::vector< std::vector<double> >, coordinate, nodalpos_xi_domain_)
//    {
//        s << "[";
//        MPFOREACH(std::vector<double>, val, coordinate)
//        {
//            s << *val << " ";
//        };
//        s << "]" << endl;
//    };
    return s.str();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void GEO::DomainIntCell::NodalPosXYZ(const DRT::Element& ele, LINALG::SerialDenseMatrix& xyz_cell) const
{
    ComputePhysicalCoordinates(ele, (*this), xyz_cell);
    return;
}

// set element nodal coordinates according to given distype
LINALG::SerialDenseMatrix GEO::DomainIntCell::GetDefaultCoordinates(
        const DRT::Element::DiscretizationType distype) const
{
    const int nsd = 3;
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    LINALG::SerialDenseMatrix coords(3,numnode);
    
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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
BlitzVec3 GEO::DomainIntCell::GetPhysicalCenterPosition(const DRT::Element& ele) const
{
    // number of space dimensions
    const int nsd = 3;
    
    // physical positions of cell nodes
    static LINALG::SerialDenseMatrix physcoord(nsd,27);
    this->NodalPosXYZ(ele, physcoord);
    
    const int numnodecell = DRT::UTILS::getNumberOfElementNodes(this->Shape());
    
    // center in local coordinates
    LINALG::Matrix<3,1> localcenterpos;
    localcenterpos = DRT::UTILS::getLocalCenterPosition<3>(this->Shape());

    // shape functions
    LINALG::SerialDenseVector funct(numnodecell);
    DRT::UTILS::shape_function_3D(funct,
            localcenterpos(0),
            localcenterpos(1),
            localcenterpos(2),
            this->Shape());
    
    //interpolate position to x-space
    BlitzVec3 x_interpol;
    x_interpol = 0.0;
    for (int inode = 0; inode < numnodecell; ++inode)
    {
        for (int isd = 0; isd < nsd; ++isd)
        {
            x_interpol(isd) += funct(inode)*physcoord(isd,inode);
        }
    }
    
    return x_interpol;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType    distype,
        const int                                 surface_ele_gid,
        const std::vector< std::vector<double> >&           domainCoordinates,
        const std::vector< std::vector<double> >&           boundaryCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_(ConvertPosArrayToLINALG<3>(domainCoordinates, distype)),
            nodalpos_xi_boundary_(ConvertPosArrayToLINALG<2>(boundaryCoordinates, distype))
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


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType    distype,
        const int                                 surface_ele_gid,
        const LINALG::SerialDenseMatrix&          domainCoordinates,
        const LINALG::SerialDenseMatrix&          boundaryCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_(domainCoordinates),
            nodalpos_xi_boundary_(boundaryCoordinates)
{
    return;
}
        
        
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GEO::BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            surface_ele_gid_(old.surface_ele_gid_),
            nodalpos_xi_domain_(old.nodalpos_xi_domain_),
            nodalpos_xi_boundary_(old.nodalpos_xi_boundary_)
{
    return;   
}
     
std::string GEO::BoundaryIntCell::toString() const
{
    std::stringstream s;
    s << "BoundaryIntCell" << endl;
    s << nodalpos_xi_domain_ << endl;
//    MCONST_FOREACH(std::vector< std::vector<double> >, coordinate, nodalpos_xi_domain_)
//    {
//        s << "[";
//        MPFOREACH(std::vector<double>, val, coordinate)
//        {
//            s << *val << " ";
//        };
//        s << "]" << endl;
//    };
    return s.str();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void GEO::BoundaryIntCell::NodalPosXYZ(const DRT::Element& ele, LINALG::SerialDenseMatrix& xyz_cell) const
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
    static LINALG::SerialDenseMatrix physcoord(3,27);
    this->NodalPosXYZ(ele, physcoord);
    
    // center in local coordinates
    const LINALG::Matrix<2,1> localcenterpos(DRT::UTILS::getLocalCenterPosition<2>(this->Shape()));

    // shape functions
    LINALG::SerialDenseVector funct(DRT::UTILS::getNumberOfElementNodes(this->Shape()));
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
