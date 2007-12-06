 /*!
 \file drt_utils_local_connectivity_matrices.cpp

 \brief Provide a node numbering scheme together with a set of shape functions

 Provided are 1D, 2D and 3D shape functions

 The surface mapping gives the node numbers such that the 2D shape functions can be used
 Nodal mappings describe the relation between volume, surface and line node numbering.
 They should be used as the only reference for such relationships.
 The corresponding graphics and a detailed description can be found in the Baci guide in the Convention chapter.
 The numbering of lower order elements is included in the higher order element, such that
 e.g. the hex8 volume element uses only the first 8 nodes of the hex27 mapping

 !!!!
 The corresponding graphics and a detailed description can be found
 in the Baci guide in the Convention chapter.
 !!!!

 \author Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15236
 */
#ifdef CCADISCRET

#include "drt_element.H"
#include "drt_discret.H"
#include "drt_utils_local_connectivity_matrices.H"
#include "drt_dserror.H"

/*----------------------------------------------------------------------*
 |  returns the number of nodes                              a.ger 11/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/    
int DRT::Utils::getNumberOfElementNodes( 
    const DRT::Element::DiscretizationType&     distype)
{
    
    int numnodes = 0;
    
    switch(distype)
    {
        case DRT::Element::hex8:
            numnodes = 8;
            break;
        case DRT::Element::hex20:
            numnodes = 20;
            break;
        case DRT::Element::hex27:
            numnodes = 27;
            break;
        case DRT::Element::tet4:
            numnodes = 4;
            break;
        case DRT::Element::tet10:
            numnodes = 10;
            break;   
        default:
            dserror("discretization type not yet implemented");     
    }
    
    return numnodes;     
}   


/*----------------------------------------------------------------------*
 |  returns the number of corner nodes                       u.may 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/    
int DRT::Utils::getNumberOfElementCornerNodes( 
    const DRT::Element::DiscretizationType&     distype)
{
    
    int numCornerNodes = 0;
    
    switch(distype)
    {
        case DRT::Element::hex8:
            numCornerNodes = 8;
            break;
        case DRT::Element::hex20:
            numCornerNodes = 8;
            break;
        case DRT::Element::hex27:
            numCornerNodes = 8;
            break;
        case DRT::Element::tet4:
            numCornerNodes = 4;
            break;
        case DRT::Element::tet10:
            numCornerNodes = 4;
            break;   
        default:
            dserror("discretization type not yet implemented");     
    }
    
    return numCornerNodes;     
}   
    


/*----------------------------------------------------------------------*
 |  returns the number of lines                              a.ger 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/    
int DRT::Utils::getNumberOfElementLines( 
    const DRT::Element::DiscretizationType&     distype)
{
    int numLines = 0;
    switch(distype)
    {
        case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
            numLines = 12;
            break;
        case DRT::Element::wedge6: case DRT::Element::wedge15:
            numLines = 9;
            break;
        case DRT::Element::tet4: case DRT::Element::tet10:
            numLines = 6;
            break;
        case DRT::Element::quad4: case DRT::Element::quad8: case DRT::Element::quad9:
            numLines = 4;
            break;
        case DRT::Element::tri3: case DRT::Element::tri6:
            numLines = 3;
            break;
        default:
            dserror("discretization type not yet implemented");     
    }
    return numLines;
}  


/*----------------------------------------------------------------------*
 |  returns the number of lines                              a.ger 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/    
int DRT::Utils::getNumberOfElementSurfaces( 
    const DRT::Element::DiscretizationType&     distype)
{
    int numSurf = 0;
    switch(distype)
    {
        case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
            numSurf = 6;
            break;
        case DRT::Element::wedge6: case DRT::Element::wedge15:
            numSurf = 5;
            break;
        case DRT::Element::tet4: case DRT::Element::tet10:
            numSurf = 4;
            break;
        default:
            dserror("discretization type not yet implemented");     
    }
    return numSurf;
}  

/*----------------------------------------------------------------------*
 |  Fills a vector< vector<int> > with all nodes for         u.may 08/07|
 |  every surface for each discretization type                          |
 *----------------------------------------------------------------------*/    
vector< vector<int> > DRT::Utils::getEleNodeNumberingSurfaces(   
    const DRT::Element::DiscretizationType&     distype)
{   
    int nSurf;
    int nNode; 
    vector< vector<int> >   map;
   
    
    switch(distype)
    {
        case DRT::Element::hex8:
        {
            nSurf = 6;
            nNode = 4;
            
            vector<int> submap(nNode, 0);
            for(int i = 0; i < nSurf; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_hex27_surfaces[i][j];
            }
            break;
        }
        case DRT::Element::hex20:
        {
            nSurf = 6;
            nNode = 8;
           
            vector<int> submap(nNode, 0);          
            for(int i = 0; i < nSurf; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_hex27_surfaces[i][j];
            }
            break;
        }
        case DRT::Element::hex27:
        {
            nSurf = 6;
            nNode = 9;
                     
            vector<int> submap(nNode, 0);
            for(int i = 0; i < nSurf; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_hex27_surfaces[i][j];
            }
            break;
        }
        case DRT::Element::tet4:
        {
            nSurf = 4;
            nNode = 3;
          
            vector<int> submap(nNode, 0);           
            for(int i = 0; i < nSurf; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_tet10_surfaces[i][j];       
            }
            break;
        }
        case DRT::Element::tet10:
        {
            nSurf = 4;
            nNode = 6;
     
            vector<int> submap(nNode, 0);
            for(int i = 0; i < nSurf; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_tet10_surfaces[i][j];
            }
            break;
        }
        default: 
            dserror("discretizationtype is not yet implemented"); 
    }
    
    return map;
}
    
    
   
/*----------------------------------------------------------------------*
 |  Fills a vector< vector<int> > with all nodes for         u.may 08/07|
 |  every line for each discretization type                             |
 *----------------------------------------------------------------------*/      
vector< vector<int> > DRT::Utils::getEleNodeNumberingLines(  
    const DRT::Element::DiscretizationType&     distype)
{   
    int nLine;
    int nNode;
    
    vector< vector<int> >  map;
        
    switch(distype)
    {
        case DRT::Element::hex8:
        {
            nLine = 12;
            nNode = 2;          
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_hex27_lines[i][j];
            }
            break;
        }
        case DRT::Element::hex20:
        {
            nLine = 12;
            nNode = 3;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {   
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_hex27_lines[i][j];
            }
            break;
        }
        case DRT::Element::hex27:
        {
            nLine = 12;
            nNode = 3;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_hex27_lines[i][j];
            }
            break;
        }
        case DRT::Element::tet4:
        {
            nLine = 6;
            nNode = 2;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_tet10_lines[i][j];
            }
            break;
        }
        case DRT::Element::tet10:
        {
            nLine = 6;
            nNode = 3;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_tet10_lines[i][j];
            }
            break;
        }
        case DRT::Element::quad9:
        {
            nLine = 4;
            nNode = 3;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_quad9_lines[i][j];
            }
            break;
        } 
        case DRT::Element::quad4:
        {
            nLine = 4;
            nNode = 2;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_quad9_lines[i][j];
            }
            break;
        }
        case DRT::Element::tri6:
        {
            nLine = 3;
            nNode = 3;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_tri6_lines[i][j];
            }
            break;
        }
        case DRT::Element::tri3:
        {
            nLine = 3;
            nNode = 2;
            vector<int> submap(nNode, 0);
            
            for(int i = 0; i < nLine; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nNode; j++)
                    map[i][j] = eleNodeNumbering_tri6_lines[i][j];
            }
            break;
        }
        default: 
            dserror("discretizationtype is not yet implemented"); 
    }
    
    return map;
}


/*----------------------------------------------------------------------*
 |  Fills a vector< vector<int> > with all surfaces for      u.may 08/07|
 |  every line for each discretization type                             |
 *----------------------------------------------------------------------*/   
vector< vector<int> > DRT::Utils::getEleNodeNumbering_lines_surfaces(    
    const DRT::Element::DiscretizationType&     distype)
{   
    int nLine;
    int nSurf;
       
    vector< vector<int> > map;
        
    if(distype == DRT::Element::hex8 ||  distype == DRT::Element::hex20 || distype == DRT::Element::hex27)
    {
        nLine = 12;
        nSurf = 2;
        vector<int> submap(nSurf, 0);
        for(int i = 0; i < nLine; i++)
        { 
            map.push_back(submap);
            for(int j = 0; j < nSurf; j++)
                map[i][j] = eleNodeNumbering_hex27_lines_surfaces[i][j];
        }
    }
    else if(distype == DRT::Element::tet4 ||  distype == DRT::Element::tet10)
    {
        nLine = 6;
        nSurf = 2;
        vector<int> submap(nSurf, 0);
        for(int i = 0; i < nLine; i++)
        { 
            map.push_back(submap);
            for(int j = 0; j < nSurf; j++)
                map[i][j] = eleNodeNumbering_tet10_lines_surfaces[i][j];
        }
    }
    else
        dserror("discretizationtype not yet implemented");
        
        
    return map;
    
}
                                  
             
                                               
/*----------------------------------------------------------------------*
 |  Fills a vector< vector<int> > with all surfaces for      u.may 08/07|
 |  every node for each discretization type                             |
 *----------------------------------------------------------------------*/   
vector< vector<int> > DRT::Utils::getEleNodeNumbering_nodes_surfaces(    
    const DRT::Element::DiscretizationType&     distype)
{
    int nNode;
    int nSurf;
       
    vector< vector<int> >   map;
        
    if(distype == DRT::Element::hex8 ||  distype == DRT::Element::hex20 || distype == DRT::Element::hex27)
    {
        nNode = 8;
        nSurf = 3;
        vector<int> submap(nSurf, 0);
        for(int i = 0; i < nNode; i++)
        {
            map.push_back(submap);
            for(int j = 0; j < nSurf; j++)
                map[i][j] = eleNodeNumbering_hex27_nodes_surfaces[i][j];
        }
    }
    else if(distype == DRT::Element::tet4 ||  distype == DRT::Element::tet10)
    {
        nNode = 4;
        nSurf = 3;
        vector<int> submap(nSurf, 0);
        for(int i = 0; i < nNode; i++)
        {
            map.push_back(submap);
            for(int j = 0; j < nSurf; j++)
                map[i][j] = eleNodeNumbering_tet10_nodes_surfaces[i][j];
        }
    }
    else
        dserror("discretizationtype not yet implemented");
        
    return map; 
   
}



/*----------------------------------------------------------------------*
 |  Fills a vector< vector<int> > with all nodes for         u.may 08/07|
 |  every surface for each discretization type                          |
 *----------------------------------------------------------------------*/   
vector< vector<double> > DRT::Utils::getEleNodeNumbering_nodes_reference(   
    const DRT::Element::DiscretizationType&     distype)
{

    int nNode;
    int nCoord;
    
    vector< vector<double> >   map;
        
    switch(distype)
    {
        case DRT::Element::hex8:
        {
            nNode = 8;
            nCoord = 3;
            vector<double> submap(nCoord, 0);
            for(int i = 0; i < nNode; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nCoord; j++)
                    map[i][j] = eleNodeNumbering_hex27_nodes_reference[i][j];
            }
            break;
        }
        case DRT::Element::hex20:
        {
            nNode = 20;
            nCoord = 3;
            vector<double> submap(nCoord, 0);
            for(int i = 0; i < nNode; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nCoord; j++)
                    map[i][j] = eleNodeNumbering_hex27_nodes_reference[i][j];
            }
            break;
        }
        case DRT::Element::hex27:
        {
            nNode = 27;
            nCoord = 3;
            vector<double> submap(nCoord, 0);
            for(int i = 0; i < nNode; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nCoord; j++)
                    map[i][j] = eleNodeNumbering_hex27_nodes_reference[i][j];
            }
            break;
        }
        case DRT::Element::tet4:
        {
            nNode = 4;
            nCoord = 3;
            vector<double> submap(nCoord, 0);
            for(int i = 0; i < nNode; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nCoord; j++)
                    map[i][j] = eleNodeNumbering_tet10_nodes_reference[i][j];
            }
            break;
        }
        case DRT::Element::tet10:
        {
            nNode = 10;
            nCoord = 3;
            vector<double> submap(nCoord, 0);
            for(int i = 0; i < nNode; i++)
            {
                map.push_back(submap);
                for(int j = 0; j < nCoord; j++)
                    map[i][j] = eleNodeNumbering_tet10_nodes_reference[i][j];
            }
            break;
        }
        default:
            dserror("discretizationtype not yet implemented");
    }
    
    return map;
}


/*----------------------------------------------------------------------*
 |  Fills an array with surface ID s a point is lying on     u.may 08/07|
 |  for each discretization type                                        |
 *----------------------------------------------------------------------*/   
int DRT::Utils::getSurfaces(    
    const Epetra_SerialDenseVector&             rst,
    int*                                        surfaces,
    const DRT::Element::DiscretizationType&     distype)
{

    int countSurf = 0;
    double TOL = 1e-7;
    
    if(distype == DRT::Element::hex8 ||  distype == DRT::Element::hex20 || distype == DRT::Element::hex27)
    {
        if(fabs(rst[0]-1.0) < TOL)      surfaces[countSurf++] = 2;        
        if(fabs(rst[0]+1.0) < TOL)      surfaces[countSurf++] = 4;        
        if(fabs(rst[1]-1.0) < TOL)      surfaces[countSurf++] = 3;        
        if(fabs(rst[1]+1.0) < TOL)      surfaces[countSurf++] = 1;        
        if(fabs(rst[2]-1.0) < TOL)      surfaces[countSurf++] = 5;        
        if(fabs(rst[2]+1.0) < TOL)      surfaces[countSurf++] = 0;     
    }
    else if(distype == DRT::Element::tet4 ||  distype == DRT::Element::tet10 )
    {
        double tetcoord = rst[0]+rst[1]+rst[2];
        if(fabs(rst[1])         < TOL)  surfaces[countSurf++] = 0;        
        if(fabs(tetcoord-1.0)   < TOL)  surfaces[countSurf++] = 1;       
        if(fabs(rst[0])         < TOL)  surfaces[countSurf++] = 2;        
        if(fabs(rst[2])         < TOL)  surfaces[countSurf++] = 3;          
    }
    else
        dserror("discretization type not yet implemented");
        
    return countSurf;   
}


/*----------------------------------------------------------------------*
 |  Fills an array with coordinates in the reference         u.may 08/07|
 |  system of the cutter element                                        |
 |  according to the node ID for each discretization type               |
 *----------------------------------------------------------------------*/   
void DRT::Utils::getNodeCoordinates(    int                                         nodeId,
                                        double*                                     coord,
                                        const DRT::Element::DiscretizationType&     distype)
{

    if(distype == DRT::Element::quad4 ||  distype == DRT::Element::quad8 || distype == DRT::Element::quad9)
    {       
        switch(nodeId)
        {
            case 0:
            {
                coord[0] = -1.0; 
                coord[1] = -1.0;
                break;
            }
            case 1:
            {
                coord[0] =  1.0; 
                coord[1] = -1.0;
                break;
            }
            case 2:
            {
                coord[0] =  1.0; 
                coord[1] =  1.0;
                break;
            }
            case 3:
            {
                coord[0] = -1.0; 
                coord[1] =  1.0;
                break;
            }
            default:
                dserror("node number not correct");    
        }
        coord[2] = 0.0; 
    }
    else if(distype == DRT::Element::tri3 ||  distype == DRT::Element::tri6)
    {       
        switch(nodeId)
        {
            case 0:
            {
                coord[0] = 0.0; 
                coord[1] = 0.0;
                break;
            }
            case 1:
            {
                coord[0] =  1.0; 
                coord[1] =  0.0;
                break;
            }
            case 2:
            {
                coord[0] =  0.0; 
                coord[1] =  1.0;
                break;
            }
            default:
                dserror("node number not correct");    
        }
        coord[2] = 0.0; 
    }
    else dserror("discretizationtype is not yet implemented");
}



/*----------------------------------------------------------------------*
 |  Fills an array with coordinates in the reference         u.may 08/07|
 |  system of the cutter element                                        |
 |  according to the line ID for each discretization type               |
 *----------------------------------------------------------------------*/   
void DRT::Utils::getLineCoordinates(    
    int                                         lineId,
    double                                      lineCoord,
    double*                                     coord,
    const DRT::Element::DiscretizationType&     distype)
{

    if(distype == DRT::Element::quad4 ||  distype == DRT::Element::quad8 || distype == DRT::Element::quad9)
    {
        // change minus sign if you change the line numbering   
        switch(lineId)
        {
            case 0:
            {
                coord[0] = lineCoord; 
                coord[1] = -1.0;
                break;
            }
            case 1:
            {
                coord[0] = 1.0; 
                coord[1] = lineCoord;
                break;
            }
            case 2:
            {
                coord[0] =  -lineCoord;    
                coord[1] =  1.0;
                break;
            }
            case 3:
            {
                coord[0] = -1.0; 
                coord[1] = lineCoord;
                break;
            }
            default:
                dserror("node number not correct");   
                           
        }  
        coord[2] =  0.0;          
    }
    else if(distype == DRT::Element::tri3 ||  distype == DRT::Element::tri6)
    {
        // change minus sign if you change the line numbering   
        switch(lineId)
        {
            case 0:
            {
                coord[0] = (lineCoord+1)*0.5; 
                coord[1] = 0.0;
                break;
            }
            case 1:
            {
                coord[0] = 1.0; 
                coord[1] = (lineCoord+1)*0.5;
                break;
            }
            case 2:
            {
                coord[0] =  1.0 - (lineCoord+1)*0.5;    
                coord[1] =  (lineCoord+1)*0.5;
                break;
            }
            default:
                dserror("node number not correct");   
                           
        }  
        coord[2] =  0.0;          
    }
    else
        dserror("discretization type not yet implemented");
}



/*----------------------------------------------------------------------*
 |  returns the index of a higher order                      u.may 09/07|
 |  element node index lying between two specified corner               |
 |  node indices for each discretizationtype                            |
 *----------------------------------------------------------------------*/  
int DRT::Utils::getHigherOrderIndex(
    const int                                   index1, 
    const int                                   index2,
    const DRT::Element::DiscretizationType&     distype )
{

    int higherOrderIndex = 0;
    
    switch(distype)
    {
        case DRT::Element::tet10:
        {
            if     ( (index1 == 0 && index2 == 1) || (index1 == 1 && index2 == 0) )      higherOrderIndex = 4;
            else if( (index1 == 1 && index2 == 2) || (index1 == 2 && index2 == 1) )      higherOrderIndex = 5;
            else if( (index1 == 2 && index2 == 0) || (index1 == 0 && index2 == 2) )      higherOrderIndex = 6;
            else if( (index1 == 0 && index2 == 3) || (index1 == 3 && index2 == 0) )      higherOrderIndex = 7;
            else if( (index1 == 1 && index2 == 3) || (index1 == 3 && index2 == 1) )      higherOrderIndex = 8;
            else if( (index1 == 2 && index2 == 3) || (index1 == 3 && index2 == 2) )      higherOrderIndex = 9;
            else dserror("no valid tet10 edge found");
            break;
        }
        case DRT::Element::quad9:
        {
            if     ( (index1 == 0 && index2 == 1) || (index1 == 1 && index2 == 0) )      higherOrderIndex = 4;
            else if( (index1 == 1 && index2 == 2) || (index1 == 2 && index2 == 1) )      higherOrderIndex = 5;
            else if( (index1 == 2 && index2 == 3) || (index1 == 3 && index2 == 2) )      higherOrderIndex = 6;
            else if( (index1 == 3 && index2 == 0) || (index1 == 0 && index2 == 3) )      higherOrderIndex = 7;
            else dserror("no valid quad9 edge found");
            break;
        }
        case DRT::Element::tri6:
        {
            if     ( (index1 == 0 && index2 == 1) || (index1 == 1 && index2 == 0) )      higherOrderIndex = 3;
            else if( (index1 == 1 && index2 == 2) || (index1 == 2 && index2 == 1) )      higherOrderIndex = 4;
            else if( (index1 == 2 && index2 == 0) || (index1 == 0 && index2 == 2) )      higherOrderIndex = 5;
            else dserror("no valid tri6 edge found");
            break;
        }
        default:
            dserror("discretizationtype not yet implemented");  
    }          
    return higherOrderIndex;
}



/*----------------------------------------------------------------------*
 |  returns the dimension of the element parameter space     u.may 10/07|
 *----------------------------------------------------------------------*/  
int DRT::Utils::getDimension(
    DRT::Element*   element)
{
    int dim = 0;
    DRT::Element::DiscretizationType distype = element->Shape();
    
    switch(distype)
    {
        case DRT::Element::line2 :  case DRT::Element::line3 :
        {
            dim = 1;      
            break;
        }
        case DRT::Element::quad4 : case DRT::Element::quad8 : case DRT::Element::quad9 :
        case DRT::Element::tri3 : case DRT::Element::tri6 :
        {
            dim = 2;   
            break;
        }
        case DRT::Element::hex8 : case DRT::Element::hex20 : case DRT::Element::hex27 :
        case DRT::Element::tet4 : case DRT::Element::tet10 :
        {
            dim = 3;      
            break;
        }   
        default:
            dserror("discretization type is not yet implemented");
    }
    return dim;
}


/*----------------------------------------------------------------------*
 |  Returns the geometric center of the element in local coordinates    |     
 |                                                           a.ger 12/07|
 *----------------------------------------------------------------------*/
std::vector<double> DRT::Utils::getLocalCenterPosition(
        const DRT::Element::DiscretizationType   distype     ///< shape of the element
        )
{
    vector<double> pos(3);
    switch(distype)
    {
        case DRT::Element::line2 :  case DRT::Element::line3 :
        {
            pos[0] = 0.0;
            break;
        }
        case DRT::Element::quad4 : case DRT::Element::quad8 : case DRT::Element::quad9 :
        {
            pos[0] = 0.0;
            pos[1] = 0.0;
            break;
        }
        case DRT::Element::tri3 : case DRT::Element::tri6 :
        {
            pos[0] = 1.0/3.0;
            pos[1] = 1.0/3.0;
            break;
        }
        case DRT::Element::hex8 : case DRT::Element::hex20 : case DRT::Element::hex27 :
        {
            pos[0] = 0.0;
            pos[1] = 0.0;
            pos[2] = 0.0;
            break;
        }
        case DRT::Element::tet4 : case DRT::Element::tet10 :
        {
            pos[0] = 1.0/4.0;
            pos[1] = 1.0/4.0;
            pos[2] = 1.0/4.0;
            break;
        }   
        default:
            dserror("discretization type is not yet implemented");
    }
    return pos;
}



#endif  // #ifdef CCADISCRET
