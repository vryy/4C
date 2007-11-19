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
            IntCell(distype),
            domainCoordinates_(domainCoordinates)
{
    volume_ = GetVolume();
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
    volume_ = DRT::Utils::getSizeInLocalCoordinates(distype);
    return;
}
        
/*----------------------------------------------------------------------*
 |  copy-ctor                                                mwgee 11/06|
 *----------------------------------------------------------------------*/
DomainIntCell::DomainIntCell(
        const DomainIntCell& old) :
            IntCell(old),
            domainCoordinates_(old.domainCoordinates_),
            volume_(old.volume_)
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

// get ratio between my parent element and myself 
// based on the "multiplication of determinants" rule
double DomainIntCell::VolumeRatio(
        const DRT::Element::DiscretizationType  parentdistype) const
{
    const double parentsize = DRT::Utils::getSizeInLocalCoordinates(parentdistype);
    return volume_/parentsize;
}

// get volume in parameter space using Gauss integration
double DomainIntCell::GetVolume() const
{    
    // worst case assumption: I am a tet10 element
    const DRT::Element::DiscretizationType distype = DRT::Element::tet10;
    const int numnode = 10;
    const int nsd = 3;
    
    // create shape function vectors 
    blitz::Array<double,1> funct(numnode);
    blitz::Array<double,2> deriv(nsd,numnode,blitz::ColumnMajorArray<2>());
    blitz::Array<double,2> xjm(nsd,nsd,blitz::ColumnMajorArray<2>());
    blitz::firstIndex i;    // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index
    blitz::thirdIndex k;    // Placeholder for the third index
    blitz::Range _  = blitz::Range::all();
    
    // get node coordinates
    blitz::Array<double,2> xyze_local(nsd,numnode,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<numnode; inode++)
    {
        xyze_local(0,inode) = domainCoordinates_[inode][0];
        xyze_local(1,inode) = domainCoordinates_[inode][1];
        xyze_local(2,inode) = domainCoordinates_[inode][2];
    }    
    
    const DRT::Utils::GaussRule3D gaussrule = DRT::Utils::intrule_tet_10point;
    
    // gaussian points
    const DRT::Utils::IntegrationPoints3D intpoints(gaussrule);

    double volume = 0.0;
    
    // integration loop
    for (int iquad=0; iquad<intpoints.nquad; ++iquad)
    {
        // coordiantes of the current integration point
        const double e1 = intpoints.qxg[iquad][0];
        const double e2 = intpoints.qxg[iquad][1];
        const double e3 = intpoints.qxg[iquad][2];

        // shape functions and their derivatives
        DRT::Utils::shape_function_3D(funct,e1,e2,e3,distype);
        DRT::Utils::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

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
        xjm = blitz::sum(deriv(i,k)*xyze_local(j,k),k);
        const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                           xjm(0,1)*xjm(1,2)*xjm(2,0)+
                           xjm(0,2)*xjm(1,0)*xjm(2,1)-
                           xjm(0,2)*xjm(1,1)*xjm(2,0)-
                           xjm(0,0)*xjm(1,2)*xjm(2,1)-
                           xjm(0,1)*xjm(1,0)*xjm(2,2);
        volume += intpoints.qwgt[iquad]*det;
    };
    return volume;
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
