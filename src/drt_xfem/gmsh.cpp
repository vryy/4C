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

using namespace std;
using namespace XFEM;

string GMSH::elementToString(const double scalar, DRT::Element* ele)
{
    //const int nen = ele->NumNode();
    DRT::Node** nodes = ele->Nodes();
    
    stringstream pos_array_string;
    pos_array_string << "SH(";
    int nendebug = 8;
    for (int i = 0; i<nendebug;++i)
    {
        const DRT::Node* node = nodes[i];
        const double* x = node->X();
        pos_array_string << scientific << x[0] << ",";
        pos_array_string << scientific << x[1] << ",";
        pos_array_string << scientific << x[2];
        if (i < nendebug-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << ")";
    // values
    pos_array_string << "{";
    for (int i = 0; i<nendebug;++i)
    {
        pos_array_string << scientific << scalar;
        if (i < nendebug-1){
            pos_array_string << ",";
        }
    };
    pos_array_string << "};";
    return pos_array_string.str();
}

string GMSH::cellToString(const double scalar, const DomainIntCell cell, DRT::Element* ele)
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
        DRT::Utils::shape_function_3D(funct,cell.GetCoord()[inen][0],
                                            cell.GetCoord()[inen][1],
                                            cell.GetCoord()[inen][2],
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

string GMSH::cellToString(const double scalar, const BoundaryIntCell cell, DRT::Element* ele)
{
    stringstream pos_array_string;
    pos_array_string << "SS(";
    
    const int nsd = 3;
    const int nen = 3;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    
    // coordinates
    for (int inen = 0; inen < nen;++inen)
    {
        // shape functions
        Epetra_SerialDenseVector funct(27);
        DRT::Utils::shape_function_3D(funct,cell.GetCoord()[inen][0],
                                            cell.GetCoord()[inen][1],
                                            cell.GetCoord()[inen][2],
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

#endif // #ifdef CCADISCRET
