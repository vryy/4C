/*!
\file integrationcell.cpp

\brief integration cell

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
#include "../drt_lib/drt_utils_integration.H"
#include "../drt_lib/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_utils_local_connectivity_matrices.H"

using namespace std;
using namespace XFEM;

#define MFOREACH(TYPE,VAL,VALS) for( TYPE::iterator VAL = VALS.begin(); VAL != VALS.end(); ++VAL )
#define MCONST_FOREACH(TYPE,VAL,VALS) for( TYPE::const_iterator VAL = VALS.begin(); VAL != VALS.end(); ++VAL )
#define MPFOREACH(TYPE,VAL,VALS) for( TYPE::const_iterator VAL = VALS->begin(); VAL != VALS->end(); ++VAL )

//
//  ctor
//
IntCell::IntCell(
        const DRT::Element::DiscretizationType distype) :
            distype_(distype)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
IntCell::IntCell(
        const IntCell& old) : 
            distype_(old.distype_)
{
    return;   
}
 
/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
IntCell::~IntCell()
{
  return;
}

//
//  get coordinates
//
vector< vector<double> > IntCell::NodalPosXiDomain() const
{
    dserror("no default implementation is given");
    vector<vector<double> > dummy;
    return dummy;
}

//
//  get coordinates in physical space
//
vector< vector<double> > IntCell::NodalPosXYZ(DRT::Element& ele) const
{
    dserror("no default implementation is given");
    vector<vector<double> > dummy;
    return dummy;
}

//
// virtual Print method
//
std::string IntCell::Print() const
{
  return "";
}


vector<vector<double> > IntCell::ComputePhysicalCoordinates(
        DRT::Element&  ele) const
{
    vector<vector<double> > physicalCoordinates;
    
    const int nsd = 3;
    const int nen_cell = DRT::UTILS::getNumberOfElementNodes(this->Shape());
    
    // coordinates
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        Epetra_SerialDenseVector funct(27);
        DRT::UTILS::shape_function_3D(
                funct,
                this->NodalPosXiDomain()[inen][0],
                this->NodalPosXiDomain()[inen][1],
                this->NodalPosXiDomain()[inen][2],
                ele.Shape());
    
        //interpolate position to x-space
        vector<double> x_interpol(nsd);
        for (int isd=0;isd < nsd;++isd ){
            x_interpol[isd] = 0.0;
        }
        
        DRT::Node** nodes = ele.Nodes();
        
        for (int inenparent = 0; inenparent < ele.NumNode();++inenparent)
        {
            const double* pos = nodes[inenparent]->X();
            for (int isd=0;isd < nsd;++isd )
            {
                x_interpol[isd] += pos[isd] * funct(inenparent);
            }
        }
        // store position
        physicalCoordinates.push_back(x_interpol);
    };
    return physicalCoordinates;
}

//
//  ctor
//
DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const vector< vector<double> >& domainCoordinates) :
            IntCell(distype),
            nodalpos_xi_domain_(domainCoordinates)
{
    return;
}

        
//
//  ctor for dummy cells
//
DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype) :
            IntCell(distype)
{
    nodalpos_xi_domain_ = GetDefaultCoordinates(distype);
    return;
}
        
/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
DomainIntCell::DomainIntCell(
        const DomainIntCell& old) :
            IntCell(old),
            nodalpos_xi_domain_(old.nodalpos_xi_domain_)
{
    return;   
}
     
string DomainIntCell::Print() const
{
    stringstream s;
    s << "DomainIntCell" << endl;
    MCONST_FOREACH(vector< vector<double> >, coordinate, nodalpos_xi_domain_)
    {
        s << "[";
        MPFOREACH(vector<double>, val, coordinate)
        {
            s << *val << " ";
        };
        s << "]" << endl;
    };
    return s.str();
}

// get volume in parameter space using Gauss integration
vector<double> DomainIntCell::modifyGaussRule3D(
        const bool standard_integration,
        const double& cell_e0,
        const double& cell_e1,
        const double& cell_e2) const
{   
    // return value
    vector<double> element_e(4);
    if (standard_integration) {
        // gauss coordinates of cell in element coordinates
        element_e[0] = cell_e0;
        element_e[1] = cell_e1;
        element_e[2] = cell_e2;
        element_e[3] = 1.0;
    } else {
    
        const DRT::Element::DiscretizationType celldistype = this->Shape();
        const int numnode = DRT::UTILS::getNumberOfElementNodes(celldistype);
        const int nsd = 3;
    
        // get node coordinates
        blitz::Array<double,2> xyze_cell(nsd,numnode,blitz::ColumnMajorArray<2>());
        for (int inode=0; inode<numnode; inode++)
            {
            xyze_cell(0,inode) = nodalpos_xi_domain_[inode][0];
            xyze_cell(1,inode) = nodalpos_xi_domain_[inode][1];
            xyze_cell(2,inode) = nodalpos_xi_domain_[inode][2];
            }    

        // init blitz indices
        blitz::firstIndex i;    // Placeholder for the first index
        blitz::secondIndex j;   // Placeholder for the second index
        blitz::thirdIndex k;    // Placeholder for the third index
    
        // create shape function vectors 
        blitz::Array<double,1> funct(numnode);
        blitz::Array<double,2> deriv(nsd,numnode,blitz::ColumnMajorArray<2>());
        DRT::UTILS::shape_function_3D(funct,cell_e0,cell_e1,cell_e2,celldistype);
        DRT::UTILS::shape_function_3D_deriv1(deriv,cell_e0,cell_e1,cell_e2,celldistype);

        // translate position into from cell coordinates to element coordinates
        const blitz::Array<double,1> e(blitz::sum(funct(j)*xyze_cell(i,j),j));


        // get Jacobian matrix and determinant
        // actually compute its transpose....
    /*
      +-            -+ T      +-            -+
      | dx   dx   dx |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dr   dr   dr |
      |              |        |              |
      | dy   dy   dy |        | dx   dy   dz |
      | --   --   -- |   =    | --   --   -- |
      | dr   ds   dt |        | ds   ds   ds |
      |              |        |              |
      | dz   dz   dz |        | dx   dy   dz |
      | --   --   -- |        | --   --   -- |
      | dr   ds   dt |        | dt   dt   dt |
      +-            -+        +-            -+
    */
        const blitz::Array<double,2> xjm(blitz::sum(deriv(i,k)*xyze_cell(j,k),k));
        const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                           xjm(0,1)*xjm(1,2)*xjm(2,0)+
                           xjm(0,2)*xjm(1,0)*xjm(2,1)-
                           xjm(0,2)*xjm(1,1)*xjm(2,0)-
                           xjm(0,0)*xjm(1,2)*xjm(2,1)-
                           xjm(0,1)*xjm(1,0)*xjm(2,2);
    
  
        // gauss coordinates of cell in element coordinates + gauss weight
        element_e[0] = e(0);
        element_e[1] = e(1);
        element_e[2] = e(2);
        element_e[3] = det;
    }  
    return element_e;
}

// set element nodal coordinates according to given distype
vector<vector<double> > DomainIntCell::GetDefaultCoordinates(
        const DRT::Element::DiscretizationType distype) const
{
    vector<vector<double> > coords;
    const int nsd = 3;
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    
    for(int j = 0; j < numnode; j++){
        vector<double> coord(nsd);
        switch (distype){
        case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        {
            for(int k = 0; k < nsd; k++){
                coord[k] = DRT::UTILS::eleNodeNumbering_hex27_nodes_reference[j][k];
                };
            break;
        }
        case DRT::Element::tet4: case DRT::Element::tet10:
        {
            for(int k = 0; k < nsd; k++){
                coord[k] = DRT::UTILS::eleNodeNumbering_tet10_nodes_reference[j][k];
                }
            break;
        }
        default:
            dserror("not supported in integrationcells. can be coded easily... ;-)");
        }
        coords.push_back(coord);  
    }
    return coords;
}

//
// return the center of the cell in physical coordinates
//
blitz::Array<double,1> DomainIntCell::GetPhysicalCenterPosition(DRT::Element& ele) const
{
    // number of space dimensions
    const int nsd = 3;
    
    //interpolate position to x-space
    blitz::Array<double, 1> x_interpol(nsd);
    x_interpol = 0.0;

    // physical positions of cell nodes
    const vector<vector<double> > physcoord = this->NodalPosXYZ(ele);
    
    // center in local coordinates
    const vector<double> localcenterpos = DRT::UTILS::getLocalCenterPosition(this->Shape());

    // shape functions
    blitz::Array<double, 1> funct(27);
    DRT::UTILS::shape_function_3D(
            funct,
            localcenterpos[0],
            localcenterpos[1],
            localcenterpos[2],
            this->Shape());
    
    for (int inen = 0; inen < this->NumNode();++inen)
    {
        for (int isd=0;isd < nsd;++isd )
        {
            x_interpol(isd) += physcoord[inen][isd] * funct(inen);
        }
    }
    // return position
    return x_interpol;
}


//
//  ctor
//
BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType    distype,
        const int                                 surface_ele_gid,
        const vector< vector<double> >&           domainCoordinates,
        const vector< vector<double> >&           boundaryCoordinates) :
            IntCell(distype),
            surface_ele_gid_(surface_ele_gid),
            nodalpos_xi_domain_(domainCoordinates),
            nodalpos_xi_boundary_(boundaryCoordinates)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            surface_ele_gid_(old.surface_ele_gid_),
            nodalpos_xi_domain_(old.nodalpos_xi_domain_),
            nodalpos_xi_boundary_(old.nodalpos_xi_boundary_)
{
    return;   
}
     
string BoundaryIntCell::Print() const
{
    stringstream s;
    s << "BoundaryIntCell" << endl;
    MCONST_FOREACH(vector< vector<double> >, coordinate, nodalpos_xi_domain_)
    {
        s << "[";
        MPFOREACH(vector<double>, val, coordinate)
        {
            s << *val << " ";
        };
        s << "]" << endl;
    };
    return s.str();
}

#endif  // #ifdef CCADISCRET
