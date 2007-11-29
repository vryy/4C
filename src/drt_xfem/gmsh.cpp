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
#include "../drt_lib/drt_utils.H"
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

string GMSH::cellToString(const double scalar, const DomainIntCell& cell, DRT::Element* ele)
{
    stringstream pos_array_string;
    pos_array_string << "SS(";
    
    const int nsd = 3;
    const int nen = 4;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    
    // coordinates
    for (int inen = 0; inen < nen;++inen)
    {
        // shape functions
        Epetra_SerialDenseVector funct(27);
        DRT::Utils::shape_function_3D(funct,cell.GetDomainCoord()[inen][0],
                                            cell.GetDomainCoord()[inen][1],
                                            cell.GetDomainCoord()[inen][2],
                                            distype);
        
        //interpolate position to x-space
        vector<double> x_interpol(nsd);
        for (int isd=0;isd < nsd;++isd ){
            x_interpol[isd] = 0.0;
        }
        for (int inenparent = 0; inenparent < ele->NumNode();++inenparent)
        {
            for (int isd=0;isd < nsd;++isd ){
                x_interpol[isd] += ele->Nodes()[inenparent]->X()[isd] * funct(inenparent);
            }
        }
        // print position in x-space
        pos_array_string << scientific << x_interpol[0] << ",";
        pos_array_string << scientific << x_interpol[1] << ",";
        pos_array_string << scientific << x_interpol[2];
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
    
    const int nsd = 3;
    const int nen = 10;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    
    // coordinates
    for (int inen = 0; inen < nen;++inen)
    {
        // shape functions
        Epetra_SerialDenseVector funct(27);
        DRT::Utils::shape_function_3D(funct,cell.GetDomainCoord()[inen][0],
                                            cell.GetDomainCoord()[inen][1],
                                            cell.GetDomainCoord()[inen][2],
                                            distype);
        
        //interpolate position to x-space
        vector<double> x_interpol(nsd);
        for (int isd=0;isd < nsd;++isd ){
            x_interpol[isd] = 0.0;
        }
        for (int inenparent = 0; inenparent < ele->NumNode();++inenparent)
        {
            for (int isd=0;isd < nsd;++isd ){
                x_interpol[isd] += ele->Nodes()[inenparent]->X()[isd] * funct(inenparent);
            }
        }
        // print position in x-space
        pos_array_string << scientific << x_interpol[0] << ",";
        pos_array_string << scientific << x_interpol[1] << ",";
        pos_array_string << scientific << x_interpol[2];
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
