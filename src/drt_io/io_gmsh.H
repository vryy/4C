/*!
\brief simple element print library for Gmsh (mostly for debuging, only sequential)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifndef GMSH_H
#define GMSH_H
#ifdef CCADISCRET

#include <blitz/array.h>
#include <string>
#include <vector>

#include "../drt_xfem/integrationcell.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_discret.H"

namespace IO
{
//! routines to transform results (elements, nodes, integrationcells, whole discretizations...) to Gmsh output
namespace GMSH
{

/*!
\brief  fill array with current nodal positions

\return array with element nodal positions
*/    
inline blitz::Array<double,2> getCurrentNodalPositions(
    const DRT::Element*                                 ele,                      ///< element with nodal pointers
    const std::map<int,blitz::TinyVector<double,3> >&   currentcutterpositions    ///< current positions of all cutter nodes
    )
{
    const int numnode = ele->NumNode();
    blitz::Array<double,2> xyze(3,numnode);
    const DRT::Node*const* nodes = ele->Nodes();
    for (int inode = 0; inode < numnode; ++inode)
    {
      const blitz::TinyVector<double,3> x = currentcutterpositions.find(nodes[inode]->Id())->second;
      xyze(0,inode) = x(0);
      xyze(1,inode) = x(1);
      xyze(2,inode) = x(2);
    }
    return xyze;
}
  
  
std::string CoordinatesToString(
        const vector<vector<double> >&         coord,
        const DRT::Element::DiscretizationType distype);
std::string CoordinatesToString(
        const blitz::Array<double,2>&          coord,
        const DRT::Element::DiscretizationType distype);

std::string ScalarToString(
        const double                           scalar,
        const DRT::Element::DiscretizationType distype);
std::string ScalarFieldToString(
        const vector<double>&                  scalarfield,
        const DRT::Element::DiscretizationType distype);
std::string VectorFieldToString(
        const vector<vector<double> >&         vectorfield,
        const DRT::Element::DiscretizationType distype);

std::string elementToString(
        const double                    scalar,
        const DRT::Element*             ele);
std::string elementToString(
        const double                            scalar,
        const DRT::Element::DiscretizationType  distype,
        const blitz::Array<double,2>&           xyze);
std::string elementToString(
        const vector<double>&           scalarfield,
        const DRT::Element*             ele);
std::string elementToString(
        const vector<vector<double> >&  vectorfield,
        const DRT::Element*             ele);

std::string cellToString(
        const blitz::Array<double,2>&           coord,
        const double                            scalar,
        const DRT::Element::DiscretizationType  distype);

//! print discretization in initial configuration (t = 0)
std::string disToString(
        const std::string&                     s,
        const double                           scalar,
        const RCP<DRT::Discretization>         dis);

//! print discretization in current configuration (t > 0)
std::string disToString(
        const std::string&                          s,
        const double                                scalar,
        const RCP<DRT::Discretization>              dis,
        std::map<int,blitz::TinyVector<double,3> >  currentpositions);

std::string disToString(
        const std::string&                     s,
        const double                           scalar,
        const RCP<DRT::Discretization>         dis,
        const map<int,XFEM::DomainIntCells >&  elementDomainIntCellsMap,
        const map<int,XFEM::BoundaryIntCells >&  elementBoundaryIntCellsMap);

std::string getConfigString(const int numview);

std::string distypeToGmshElementHeader(const DRT::Element::DiscretizationType distype);

int distypeToGmshNumNode(const DRT::Element::DiscretizationType distype);

std::string trifaceToString(
        const double scalar,
        const vector<vector<double> >& triface);

std::string XAABBToString(
        const double scalar,
        const vector<vector<double> >& XAABB);

}
}

#endif // #ifdef CCADISCRET
#endif /*GMSH_H*/
