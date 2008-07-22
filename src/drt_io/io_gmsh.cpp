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
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"

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
    const double                                      scalar, 
    const DRT::Element*                               ele,
    const map<int,blitz::TinyVector<double,3> >&      currentelepositions)
{

  const DRT::Element::DiscretizationType distype = ele->Shape();
  std::stringstream gmshfilecontent;
  
  blitz::Array<double,2> xyze(DRT::UTILS::getCurrentNodalPositions(ele,currentelepositions));
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
    blitz::Array<double,2> xyze(DRT::UTILS::getCurrentNodalPositions(actele,
        currentpositions));
    gmshfilecontent << IO::GMSH::cellWithScalarToString(actele->Shape(),
        scalar, xyze) << endl;
  };
  gmshfilecontent << "};" << endl;
  return gmshfilecontent.str();
}


#endif // #ifdef CCADISCRET
