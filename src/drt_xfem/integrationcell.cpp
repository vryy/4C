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

////
////  get coordinates
////
//vector< vector<double> >  IntCell::GetDomainCoord() const
//{
//    return domainCoordinates_;   
//}

//
// virtual Print method
//
std::string IntCell::Print() const
{
  return "";
}



//
//  ctor
//
DomainIntCell::DomainIntCell(
        const DRT::Element::DiscretizationType distype,
        const vector< vector<double> > domainCoordinates) :
            IntCell(DRT::Element::tet4),
            domainCoordinates_(domainCoordinates)
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
    SetDefaultCoordinates(distype);
    return;
}
        
/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
DomainIntCell::DomainIntCell(
        const DomainIntCell& old) :
            IntCell(old),
            domainCoordinates_(old.domainCoordinates_)
{
    return;   
}
     
//
//  get coordinates
//
vector< vector<double> >  DomainIntCell::GetDomainCoord() const
{
    return domainCoordinates_;   
}

string DomainIntCell::Print() const
{
    stringstream s;
    s << "DomainIntCell" << endl;
    MCONST_FOREACH(vector< vector<double> >, coordinate, domainCoordinates_)
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
const vector<double> DomainIntCell::modifyGaussRule3D(
        const double cell_e0,
        const double cell_e1,
        const double cell_e2) const
{   
    // return value
    vector<double> element_e(4);
    
    const DRT::Element::DiscretizationType celldistype = this->Shape();
    const int numnode = DRT::Utils::getNumberOfElementNodes(celldistype);
    const int nsd = 3;
    
    // get node coordinates
    blitz::Array<double,2> xyze_cell(nsd,numnode,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<numnode; inode++)
    {
        xyze_cell(0,inode) = domainCoordinates_[inode][0];
        xyze_cell(1,inode) = domainCoordinates_[inode][1];
        xyze_cell(2,inode) = domainCoordinates_[inode][2];
    }    

    // init blitz indices
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index
    
    // create shape function vectors 
    blitz::Array<double,1> funct(numnode);
    blitz::Array<double,2> deriv(nsd,numnode,blitz::ColumnMajorArray<2>());
    DRT::Utils::shape_function_3D(funct,cell_e0,cell_e1,cell_e2,celldistype);
    DRT::Utils::shape_function_3D_deriv1(deriv,cell_e0,cell_e1,cell_e2,celldistype);

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
    blitz::Array<double,2> xjm(nsd,nsd,blitz::ColumnMajorArray<2>());
    xjm = blitz::sum(deriv(i,k)*xyze_cell(j,k),k);
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                       xjm(0,1)*xjm(1,2)*xjm(2,0)+
                       xjm(0,2)*xjm(1,0)*xjm(2,1)-
                       xjm(0,2)*xjm(1,1)*xjm(2,0)-
                       xjm(0,0)*xjm(1,2)*xjm(2,1)-
                       xjm(0,1)*xjm(1,0)*xjm(2,2);
    
  
    // gauss coordinates of cell in element coordinates
    element_e[0] = e(0);
    element_e[1] = e(1);
    element_e[2] = e(2);
    element_e[3] = det;
    
    return element_e;
}

// set element nodal coordinates according to given distype
void DomainIntCell::SetDefaultCoordinates(
        const DRT::Element::DiscretizationType distype)
{
    
    const int nsd = 3;
    const int numnode = DRT::Utils::getNumberOfElementNodes(distype);
    
    domainCoordinates_.clear();
    for(int j = 0; j < numnode; j++){
        vector<double> coord(3);
        switch (distype){
        case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        {
            for(int k = 0; k < nsd; k++){
                coord[k] = DRT::Utils::eleNodeNumbering_hex27_nodes_reference[j][k];
                };
            break;
        }
        case DRT::Element::tet4: case DRT::Element::tet10:
        {
            for(int k = 0; k < nsd; k++){
                coord[k] = DRT::Utils::eleNodeNumbering_tet10_nodes_reference[j][k];
                }
            break;
        }
        default:
            dserror("not supported in integrationcells. can be coded easily... ;-)");
        }
        domainCoordinates_.push_back(coord);  
    }
    return;
}

//
//  ctor
//
BoundaryIntCell::BoundaryIntCell(
        const DRT::Element::DiscretizationType distype,
        const vector< vector<double> > domainCoordinates,
        const vector< vector<double> > boundaryCoordinates) :
            IntCell(distype),
            domainCoordinates_(domainCoordinates),
            boundaryCoordinates_(boundaryCoordinates_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
BoundaryIntCell::BoundaryIntCell(
        const BoundaryIntCell& old) :
            IntCell(old),
            domainCoordinates_(old.domainCoordinates_),
            boundaryCoordinates_(old.boundaryCoordinates_)
{
    return;   
}
     
//
//  get coordinates
//
vector< vector<double> >  BoundaryIntCell::GetDomainCoord() const
{
    return domainCoordinates_;   
}

//
//  get coordinates
//
vector< vector<double> >  BoundaryIntCell::GetBoundaryCoord() const
{
    return boundaryCoordinates_;   
}


string BoundaryIntCell::Print() const
{
    stringstream s;
    s << "BoundaryIntCell" << endl;
    MCONST_FOREACH(vector< vector<double> >, coordinate, domainCoordinates_)
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
