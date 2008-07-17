/*!
\file io_gmsh.cpp

\brief simple element print library for Gmsh (debuging only)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include <string>
#include <blitz/array.h>

#include "io_gmsh.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"

std::string IO::GMSH::ScalarToString(const double scalar,
    const DRT::Element::DiscretizationType distype)
{
  std::stringstream pos_array_string;

  const int numnode = distypeToGmshNumNode(distype);

  // values
  pos_array_string << "{";
  for (int i = 0; i<numnode; ++i)
  {
    pos_array_string << scientific << scalar;
    if (i < numnode-1)
    {
      pos_array_string << ",";
    }
  };
  pos_array_string << "};";
  return pos_array_string.str();
}

std::string IO::GMSH::elementAtInitialPositionToString(const double scalar, const DRT::Element* ele)
{
  const DRT::Node*const* nodes = ele->Nodes();

  const DRT::Element::DiscretizationType distype = ele->Shape();
  const int numnode = distypeToGmshNumNode(distype);

  std::stringstream pos_array_string;
  pos_array_string << "S" << distypeToGmshElementHeader(distype) << "(";
  for (int i = 0; i<numnode; ++i)
  {
    const DRT::Node* node = nodes[i];
    const double* x = node->X();
    pos_array_string << scientific << x[0] << ",";
    pos_array_string << scientific << x[1] << ",";
    pos_array_string << scientific << x[2];
    if (i < numnode-1)
    {
      pos_array_string << ",";
    }
  };
  pos_array_string << ")";
  // values
  pos_array_string << ScalarToString(scalar, distype);

  return pos_array_string.str();
}


std::string IO::GMSH::elementAtCurrentPositionToString(
    const double scalar, 
    const DRT::Element* ele,
    const map<int,BlitzVec3>&      currentelepositions)
{

  const DRT::Element::DiscretizationType distype = ele->Shape();
  std::stringstream gmshfilecontent;
  
  blitz::Array<double,2> xyze(XFEM::getCurrentNodalPositions(ele,currentelepositions));
  gmshfilecontent << IO::GMSH::cellWithScalarToString(distype, scalar, xyze) << endl;
 
  return gmshfilecontent.str();
}



std::string IO::GMSH::disToString(
    const std::string& s,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis)
{
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" " << s << " Elements \" {" << endl;
  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    const DRT::Element* actele = dis->lColElement(i);
    gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(scalar, actele) << endl;
  };
  gmshfilecontent << "};" << endl;
  return gmshfilecontent.str();
}

std::string IO::GMSH::disToString(
    const std::string& s,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis,
    std::map<int,blitz::TinyVector<double,3> > currentpositions)
{
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" " << s << " Elements \" {" << endl;
  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    const DRT::Element* actele = dis->lColElement(i);
    blitz::Array<double,2> xyze(XFEM::getCurrentNodalPositions(actele,
        currentpositions));
    gmshfilecontent << IO::GMSH::cellWithScalarToString(actele->Shape(),
        scalar, xyze) << endl;
  };
  gmshfilecontent << "};" << endl;
  return gmshfilecontent.str();
}

std::string IO::GMSH::disToString(
    const std::string& s,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis,
    const std::map<int, XFEM::DomainIntCells >& elementDomainIntCellsMap,
    const std::map<int, XFEM::BoundaryIntCells >& elementBoundaryIntCellsMap)
{
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" " << s << " Elements and Integration Cells \" {"
      << endl;
  
  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    const DRT::Element* actele = dis->lColElement(i);
    const int id = actele->Id();
    // print integration cells, if available
    if (elementDomainIntCellsMap.count(id) > 0)
    {
      const XFEM::DomainIntCells cells = elementDomainIntCellsMap.find(id)->second;
      for (XFEM::DomainIntCells::const_iterator cell = cells.begin(); cell
          != cells.end(); ++cell)
      {
        static BlitzMat dxyz_ele(3, 27);
        cell->NodalPosXYZ(*actele, dxyz_ele);
        gmshfilecontent << IO::GMSH::cellWithScalarToString(cell->Shape(),
            scalar, dxyz_ele) << endl;
      }
      const XFEM::BoundaryIntCells bcells = elementBoundaryIntCellsMap.find(id)->second;
      for (XFEM::BoundaryIntCells::const_iterator bcell = bcells.begin(); bcell
          != bcells.end(); ++bcell)
      {
        static BlitzMat bxyz_ele(3, 9);
        bcell->NodalPosXYZ(*actele, bxyz_ele);
        gmshfilecontent << IO::GMSH::cellWithScalarToString(bcell->Shape(),
            scalar, bxyz_ele) << endl;
      }
    }
    // print element
    gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(scalar, actele) << endl;
  };
  gmshfilecontent << "};" << endl;
  return gmshfilecontent.str();
}

std::string IO::GMSH::getConfigString(const int numview)
{
  std::stringstream gmshfilecontent;
  for (int iview = 0; iview < numview; ++iview)
  {
    gmshfilecontent << "View["<<iview
        <<"].RangeType = 2;   // Value scale range type (1=default, 2=custom, 3=per time step)"
        << endl;
    gmshfilecontent << "View["<<iview
        <<"].CustomMax = 1.0; // User-defined maximum value to be displayed"
        << endl;
    gmshfilecontent << "View["<<iview
        <<"].CustomMin = 0.0; // User-defined minimum value to be displayed"
        << endl;
  }
  return gmshfilecontent.str();
}


#endif // #ifdef CCADISCRET
