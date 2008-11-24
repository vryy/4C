/*!
\file io_gmsh.cpp

\brief simple element print library for Gmsh (debugging only)

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#ifdef CCADISCRET

#include "io_gmsh.H"
#include "../drt_geometry/intersection_service.H"


std::string IO::GMSH::distypeToGmshElementHeader(
    const DRT::Element::DiscretizationType distype)
{
  switch (distype)
  {
    case DRT::Element::hex8:
      return "H";
      break;
    case DRT::Element::hex20:
      return "H";
      break;
    case DRT::Element::hex27:
      return "H";
      break;
    case DRT::Element::tet4:
      return "S";
      break;
    case DRT::Element::tet10:
      return "S";
      break;
    case DRT::Element::point1:
      return "P";
      break;
    case DRT::Element::quad4:
      return "Q";
      break;
    case DRT::Element::quad8:
      return "Q";
      break;
    case DRT::Element::quad9:
      return "Q";
      break;
    case DRT::Element::tri3:
      return "T";
      break;
    case DRT::Element::tri6:
      return "T";
      break;
    case DRT::Element::line2:
      return "L";
      break;
    case DRT::Element::line3:
      return "L2";
      break;
    default:
      dserror("distypeToGmshElementHeader: distype not supported for printout!");
  }
  return "xxx";
}

int IO::GMSH::distypeToGmshNumNode(
    const DRT::Element::DiscretizationType distype   ///< element shape
)
{
  switch (distype)
  {
    case DRT::Element::hex8:
      return 8;
      break;
    case DRT::Element::hex20:
      return 8;
      break;
    case DRT::Element::hex27:
      return 8;
      break;
    case DRT::Element::tet4:
      return 4;
      break;
    case DRT::Element::tet10:
      return 4;
      break;
    case DRT::Element::point1:
      return 1;
      break;
    case DRT::Element::quad4:
      return 4;
      break;
    case DRT::Element::quad8:
      return 4;
      break;
    case DRT::Element::quad9:
      return 4;
      break;
    case DRT::Element::tri3:
      return 3;
      break;
    case DRT::Element::tri6:
      return 3;
      break;
    case DRT::Element::line2:
      return 2;
      break;
    case DRT::Element::line3:
      return 3;
      break;
    default:
      dserror("distypeToGmshNumNode: distype not supported for printout!");
  }
  return -1;
}

std::string IO::GMSH::ScalarToString(
    const double scalar,
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
    const double                            scalar, 
    const DRT::Element*                     ele,
    const map<int,LINALG::Matrix<3,1> >&    currentelepositions)
{

  const DRT::Element::DiscretizationType distype = ele->Shape();
  
  std::stringstream gmshfilecontent;
  gmshfilecontent << IO::GMSH::cellWithScalarToString(
      distype, scalar, GEO::getCurrentNodalPositions(ele,currentelepositions)
      ) << "\n";
  return gmshfilecontent.str();
}



std::string text3dToString(
    const LINALG::Matrix<3,1>&            xyz,      ///< 3d Position of text
    const std::string                     text,     ///< text to be printed 
    const int                             fontsize  ///< font size
    )
{
  std::stringstream gmsh_ele_line;

  gmsh_ele_line << "T3";
  // coordinates
  gmsh_ele_line << "("<< scientific << xyz(0)<<",";
  gmsh_ele_line << scientific << xyz(1)<<",";
  gmsh_ele_line << scientific << xyz(2) <<",";             
  gmsh_ele_line << fontsize << ")";
  gmsh_ele_line <<"{\"" << text <<"\"};";

  return gmsh_ele_line.str();
}

std::string IO::GMSH::disToString(
    const std::string& s,
    const double scalar,
    const Teuchos::RCP<DRT::Discretization> dis)
{
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" " << s << " Elements \" {\n";
  for (int i=0; i<dis->NumMyRowElements(); ++i)
  {
    const DRT::Element* actele = dis->lRowElement(i);
    gmshfilecontent << IO::GMSH::elementAtInitialPositionToString(scalar, actele) << "\n";
  };
  gmshfilecontent << "};\n";
  return gmshfilecontent.str();
}

std::string IO::GMSH::disToString(
    const std::string&                      	s,
    const double                            	scalar,
    const Teuchos::RCP<DRT::Discretization> 	dis,
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions)
{
  std::stringstream gmshfilecontent;
  gmshfilecontent << "View \" " << s << " Elements \" {\n";

  for (int i=0; i<dis->NumMyColElements(); ++i)
  {
    const DRT::Element* actele = dis->lColElement(i);
    gmshfilecontent << IO::GMSH::cellWithScalarToString(actele->Shape(),
        scalar, GEO::getCurrentNodalPositions(actele,currentpositions) ) << "\n";
  };
  gmshfilecontent << "};\n";
  return gmshfilecontent.str();
}


#endif // #ifdef CCADISCRET
