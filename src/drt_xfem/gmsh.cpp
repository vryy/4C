/*!
\file gmsh.cpp

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

#include "gmsh.H"
#include "integrationcell.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"

using namespace std;
using namespace XFEM;

string GMSH::elementToString(const double scalar, DRT::Element* ele)
{
    DRT::Node** nodes = ele->Nodes();
    
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = distypeToGmshNumNode(distype);
    
    stringstream pos_array_string;
    pos_array_string << "S" << distypeToGmshElementHeader(distype) << "(";
    for (int i = 0; i<numnode;++i)
    {
        const DRT::Node* node = nodes[i];
        const double* x = node->X();
        pos_array_string << scientific << x[0] << ",";
        pos_array_string << scientific << x[1] << ",";
        pos_array_string << scientific << x[2];
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<numnode;++i)
    {
        pos_array_string << scientific << scalar;
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
}

string GMSH::elementToString(const vector<double> scalarfield, DRT::Element* ele)
{
    DRT::Node** nodes = ele->Nodes();
    
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = distypeToGmshNumNode(distype);
    
    if ((unsigned int)(scalarfield.size()) != (unsigned int)numnode) dserror("Size mismatch: No of Nodes vs Size of Scalarfield");
    
    stringstream pos_array_string;
    pos_array_string << "S" << distypeToGmshElementHeader(distype) << "(";
    for (int i = 0; i<numnode;++i)
    {
        const DRT::Node* node = nodes[i];
        const double* x = node->X();
        pos_array_string << scientific << x[0] << ",";
        pos_array_string << scientific << x[1] << ",";
        pos_array_string << scientific << x[2];
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<numnode;++i)
    {
        pos_array_string << scientific << scalarfield[i];
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
}

string GMSH::elementToString(const vector<vector<double> > vectorfield, DRT::Element* ele)
{
    DRT::Node** nodes = ele->Nodes();
    
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = distypeToGmshNumNode(distype);
    const int numcomp=3;
    
    if ((unsigned int)vectorfield.size() != (unsigned int)numnode) dserror("Size mismatch: No of Nodes vs Size of Vectorfield");
    if ((unsigned int)vectorfield[0].size() != 3) dserror("Size mismatch: Vector of Vectorfield must have length 3");
    
    stringstream pos_array_string;
    pos_array_string << "S" << distypeToGmshElementHeader(distype) << "(";
    for (int i = 0; i<numnode;++i)
    {
        const DRT::Node* node = nodes[i];
        const double* x = node->X();
        pos_array_string << scientific << x[0] << ",";
        pos_array_string << scientific << x[1] << ",";
        pos_array_string << scientific << x[2];
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<numnode;++i)
    {
        for (int component = 0; component < numcomp; ++component) {
          pos_array_string << scientific << vectorfield[i][component];
          if (component < numcomp-1) pos_array_string << ",";
        }
        if (i < numnode-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
}

string GMSH::cellToString(const double scalar, const DomainIntCell& cell, DRT::Element* ele)
{
    stringstream pos_array_string;
    pos_array_string << "SS(";
    
    const int nen = 4;
    
    // coordinates
    for (int inen = 0; inen < nen;++inen)
    {
        // print position in x-space
        pos_array_string << scientific << cell.GetPhysicalCoord(*ele)[inen][0] << ",";
        pos_array_string << scientific << cell.GetPhysicalCoord(*ele)[inen][1] << ",";
        pos_array_string << scientific << cell.GetPhysicalCoord(*ele)[inen][2];
        if (inen < nen-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<nen;++i)
    {
        pos_array_string << scientific << scalar;
        if (i < nen-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
};

string GMSH::cellToString(const double scalar, const BoundaryIntCell& cell, DRT::Element* ele)
{
    stringstream pos_array_string;
    pos_array_string << "SS2(";
    //pos_array_string << "SS(";
    
    const int nen = 10;
    
    // coordinates
    for (int inen = 0; inen < nen;++inen)
    {
        // print position in x-space
        pos_array_string << scientific << cell.GetPhysicalCoord(*ele)[inen][0] << ",";
        pos_array_string << scientific << cell.GetPhysicalCoord(*ele)[inen][1] << ",";
        pos_array_string << scientific << cell.GetPhysicalCoord(*ele)[inen][2];
        if (inen < nen-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<nen;++i)
    {
        pos_array_string << scientific << scalar;
        if (i < nen-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
};

string GMSH::disToString(const std::string s, const double scalar, 
        RefCountPtr<DRT::Discretization> dis)
{
    stringstream gmshfilecontent;
    gmshfilecontent << "View \" " << s << " Elements \" {" << endl;
    for (int i=0; i<dis->NumMyColElements(); ++i)
    {
        DRT::Element* actele = dis->lColElement(i);
        gmshfilecontent << GMSH::elementToString(scalar, actele) << endl;
    };
    gmshfilecontent << "};" << endl;
    return gmshfilecontent.str();
}

string GMSH::disToString(const std::string s, const double scalar, 
        RefCountPtr<DRT::Discretization> dis,
        map<int, DomainIntCells >& elementDomainIntCellsMap)
{
    stringstream gmshfilecontent;
    gmshfilecontent << "View \" " << s << " Elements and Integration Cells \" {" << endl;
    for (int i=0; i<dis->NumMyColElements(); ++i)
    {
        DRT::Element* actele = dis->lColElement(i);
        const int id = actele->Id();
        if (elementDomainIntCellsMap.count(id))
        {
            XFEM::DomainIntCells::const_iterator cell;
            for(cell = elementDomainIntCellsMap[id].begin(); cell != elementDomainIntCellsMap[id].end(); ++cell )
            {
                gmshfilecontent << GMSH::cellToString(scalar, *cell, actele) << endl;
            }
        }
        else
        {
            gmshfilecontent << GMSH::elementToString(scalar, actele) << endl;
        };
    };
    gmshfilecontent << "};" << endl;
    return gmshfilecontent.str();
}

std::string GMSH::getConfigString(const int numview)
{
    stringstream gmshfilecontent;
    for (int iview = 0; iview < numview; ++iview) {
        gmshfilecontent << "View["<<iview<<"].RangeType = 2;   // Value scale range type (1=default, 2=custom, 3=per time step)" << endl;
        gmshfilecontent << "View["<<iview<<"].CustomMax = 1.0; // User-defined maximum value to be displayed" << endl;
        gmshfilecontent << "View["<<iview<<"].CustomMin = 0.0; // User-defined minimum value to be displayed" << endl;  
    }
    return gmshfilecontent.str();
}

std::string GMSH::distypeToGmshElementHeader(const DRT::Element::DiscretizationType distype)
{
	switch (distype){
	case DRT::Element::point1: return "P";   break;
	case DRT::Element::quad4:  return "Q";   break;
	case DRT::Element::quad9:  return "Q";  break;
	case DRT::Element::tri3:   return "T";   break;
	case DRT::Element::tri6:   return "T";  break;
	case DRT::Element::hex8:   return "H";   break;
	case DRT::Element::hex27:  return "H";  break;
	case DRT::Element::tet4:   return "S";   break;
	case DRT::Element::tet10:  return "S";  break;
	default:
		dserror("distype not supported for printout!");
	}
	return "xxx";
}

int GMSH::distypeToGmshNumNode(const DRT::Element::DiscretizationType distype)
{
	switch (distype){
	case DRT::Element::point1: return 1;   break;
	case DRT::Element::quad4:  return 4;   break;
	case DRT::Element::quad9:  return 4;   break;
	case DRT::Element::tri3:   return 3;   break;
	case DRT::Element::tri6:   return 3;   break;
	case DRT::Element::hex8:   return 8;   break;
	case DRT::Element::hex20:  return 8;   break;
	case DRT::Element::hex27:  return 8;   break;
	case DRT::Element::tet4:   return 4;   break;
	case DRT::Element::tet10:  return 4;   break;
	default:
		dserror("distype not supported for printout!");
	}
	return -1;
}

#endif // #ifdef CCADISCRET
