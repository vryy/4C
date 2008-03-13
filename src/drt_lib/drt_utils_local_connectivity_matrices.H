/*!
 \file drt_utils_local_connectivity_matrices.H

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
#ifndef DRT_UTILS_LOCAL_CONNECTIVITY_MATRICES_H
#define DRT_UTILS_LOCAL_CONNECTIVITY_MATRICES_H

#include <blitz/array.h>

#include "drt_element.H"
#include "Epetra_Vector.h"

using namespace std;
using namespace Teuchos;

namespace DRT
{
namespace UTILS
{

  //! 6 Surfaces of a Hex27 element with 9 nodes per surface
  const int eleNodeNumbering_hex27_surfaces[6][9] = {
          {0,  3,  2,  1, 11, 10,  9,  8, 20},
          {0,  1,  5,  4,  8, 13, 16, 12, 21},
          {1,  2,  6,  5,  9, 14, 17, 13, 22},
          {2,  3,  7,  6, 10, 15, 18, 14, 23},
          {0,  4,  7,  3, 12, 19, 15, 11, 24},
          {4,  5,  6,  7, 16, 17, 18, 19, 25}};

  //! 12 Lines of a Hex27 element with 3 nodes per line
  const int eleNodeNumbering_hex27_lines[12][3] = {
          {0,  1,  8},
          {1,  2,  9},
          {2,  3, 10},
          {0,  3, 11},
          {0,  4, 12},
          {1,  5, 13},
          {2,  6, 14},
          {3,  7, 15},
          {4,  5, 16},
          {5,  6, 17},
          {6,  7, 18},
          {4,  7, 19}};

  //! 4 Surfaces of a Tet10 element with 6 nodes per surface
  const int eleNodeNumbering_tet10_surfaces[4][6] = {
          {0,  1,  3,  4,  8,  7},
          {1,  2,  3,  5,  9,  8},
          {0,  3,  2,  7,  9,  6},
          {0,  2,  1,  6,  5,  4}};

  //! 6 Lines of a Tet10 element with 3 nodes per line
  const int eleNodeNumbering_tet10_lines[6][3] = {
          {0,  1,  4},
          {1,  2,  5},
          {0,  2,  6},
          {0,  3,  7},
          {1,  3,  8},
          {2,  3,  9}};

  //! Triangular surfaces of wedge15 with 7 nodes per surface
  const int eleNodeNumbering_wedge15_trisurfaces[2][6] = {
          { 0, 1, 2,  6,  7,  8},
          { 3, 4, 5, 12, 13, 14}};

  //! Rectangular surfaces of wedge15 with 8 nodes per surface
  const int eleNodeNumbering_wedge15_quadsurfaces[3][8] = {
          {0, 1, 4, 3, 6, 10, 12,  9},
          {1, 2, 5, 4, 7, 11, 13, 10},
          {0, 3, 5, 2, 9, 14, 11,  8}};

  //! 9 lines of a wedge15 element with 3 nodes per line
  const int eleNodeNumbering_wedge15_lines[9][3] = {
          {0, 1, 6},
          {1, 2, 7},
          {2, 0, 8},
          {3, 4, 12},
          {4, 5, 13},
          {5, 3, 14},
          {1, 4, 10},
          {2, 5, 11},
          {0, 3, 9}};

  //! Rectangular surface of pyramid5
  const int eleNodeNumbering_pyramid5_quadsurfaces[1][4]= {
          {0, 1, 2, 3}};
  
  //! Triangular surfaces of pyramid5
  const int eleNodeNumbering_pyramid5_trisurfaces[4][3]={
    {0,1,4},
    {1,2,4},
    {2,3,4},
    {0,4,3}};

  //! 4 Lines of a Quad9 element with 3 nodes per line
  const int eleNodeNumbering_quad9_lines[4][3] = {
          {0,  1,  4},
          {1,  2,  5},
          {2,  3,  6},
          {3,  0,  7}};

  //! 3 Lines of a Tri6 element with 3 nodes per line
  const int eleNodeNumbering_tri6_lines[3][3] = {
          {0,  1,  3},
          {1,  2,  4},
          {2,  0,  5}};

  //! for each of the 12 lines, tell me the 2 neighbouring faces
  const int eleNodeNumbering_hex27_lines_surfaces[12][2] = {
          {0,  1},
          {0,  2},
          {0,  3},
          {0,  4},
          {1,  4},
          {1,  2},
          {2,  3},
          {3,  4},
          {1,  5},
          {2,  5},
          {3,  5},
          {4,  5}};


  //! for each of the 6 lines, tell me the 2 neighbouring faces
  const int eleNodeNumbering_tet10_lines_surfaces[6][2] = {
          {0,  3},
          {1,  3},
          {2,  3},
          {2,  0},
          {0,  1},
          {1,  2}};


  //! for each of the 8 corner nodes, tell me the 3 neighbouring faces
  const int eleNodeNumbering_hex27_nodes_surfaces[8][3] =  { {0,  1,  4},
                                                             {0,  1,  2},
                                                             {0,  2,  3},
                                                             {0,  3,  4},
                                                             {1,  4,  5},
                                                             {1,  2,  5},
                                                             {2,  3,  5},
                                                             {3,  4,  5}};


   //! for each of the 4 corner nodes, tell me the 3 neighbouring faces
  const int eleNodeNumbering_tet10_nodes_surfaces[4][3] =  { {0,  2,  3},
                                                             {0,  1,  3},
                                                             {1,  2,  3},
                                                             {0,  1,  2}};


  //! for each node the reference coordinates are stored for hex27
  const double eleNodeNumbering_hex27_nodes_reference[27][3] = {    {-1.0,  -1.0,  -1.0},
                                                                    { 1.0,  -1.0,  -1.0},
                                                                    { 1.0,   1.0,  -1.0},
                                                                    {-1.0,   1.0,  -1.0},

                                                                    {-1.0,  -1.0,   1.0},
                                                                    { 1.0,  -1.0,   1.0},
                                                                    { 1.0,   1.0,   1.0},
                                                                    {-1.0,   1.0,   1.0},

                                                                    { 0.0,  -1.0,  -1.0},
                                                                    { 1.0,   0.0,  -1.0},
                                                                    { 0.0,   1.0,  -1.0},
                                                                    {-1.0,   0.0,  -1.0},

                                                                    {-1.0,  -1.0,   0.0},
                                                                    { 1.0,  -1.0,   0.0},
                                                                    { 1.0,   1.0,   0.0},
                                                                    {-1.0,   1.0,   0.0},

                                                                    { 0.0,  -1.0,   1.0},
                                                                    { 1.0,   0.0,   1.0},
                                                                    { 0.0,   1.0,   1.0},
                                                                    {-1.0,   0.0,   1.0},

                                                                    { 0.0,   0.0,  -1.0},
                                                                    { 0.0,  -1.0,   0.0},
                                                                    { 1.0,   0.0,   0.0},
                                                                    { 0.0,   1.0,   0.0},

                                                                    {-1.0,   0.0,   0.0},
                                                                    { 0.0,   0.0,   1.0},
                                                                    { 0.0,   0.0,   0.0}};


  //! for each node the reference coordinates are stored for a tet10
  const double eleNodeNumbering_tet10_nodes_reference[10][3] = {    { 0.0,   0.0,   0.0},
                                                                    { 1.0,   0.0,   0.0},
                                                                    { 0.0,   1.0,   0.0},
                                                                    { 0.0,   0.0,   1.0},

                                                                    { 0.5,   0.0,   0.0},
                                                                    { 0.5,   0.5,   0.0},
                                                                    { 0.0,   0.5,   0.0},
                                                                    { 0.0,   0.0,   0.5},
                                                                    { 0.5,   0.0,   0.5},
                                                                    { 0.0,   0.5,   0.5}};

  /*!
 \brief Returns the number of nodes
        for each discretization type
        
 \return number of nodes of the discretization type
 */ 
 int getNumberOfElementNodes( 
   const DRT::Element::DiscretizationType&     distype
   );
  
  /*!
  \brief Returns the number of nodes
         for each discretization type
         - template version
       
  \return number of nodes of the discretization type
  */  
  template <DRT::Element::DiscretizationType DISTYPE>
  int getNumberOfElementNodes()
  {
     switch(DISTYPE)
     {
     case DRT::Element::dis_none:     return 0;    break;
     case DRT::Element::point1:       return 1;    break;
     case DRT::Element::line2:        return 2;    break;
     case DRT::Element::line3:        return 3;    break;
     case DRT::Element::tri3:         return 3;    break;
     case DRT::Element::tri6:         return 6;    break;
     case DRT::Element::hex8:         return 8;    break;
     case DRT::Element::hex20:        return 20;   break;
     case DRT::Element::hex27:        return 27;   break;
     case DRT::Element::tet4:         return 4;    break;
     case DRT::Element::tet10:        return 10;   break;   
     default:
         dserror("discretization type not yet implemented");    
         return -1;
     }
  }
 
   /*!
  \brief Returns the number of corner nodes
         for each discretization type
         
  \return number of corner nodes of the discretization type
  */ 
  int getNumberOfElementCornerNodes( 
    const DRT::Element::DiscretizationType&     distype
    );

  /*!
 \brief Returns the number of lines
        for each discretization type
 */
  int getNumberOfElementLines(
      const DRT::Element::DiscretizationType&     distype);

  /*!
 \brief Returns the number of surfaces
        for each discretization type
 */
  int getNumberOfElementSurfaces(
      const DRT::Element::DiscretizationType&     distype);


  /*!
  \brief Fills a vector< vector<int> > with all nodes for every surface
         for each discretization type
         
  \return map with all nodes for each surface
  */
  vector< vector<int> > getEleNodeNumberingSurfaces( 
    const DRT::Element::DiscretizationType&     distype     ///< discretization type
    );

  /*!
  \brief Fills a vector< vector<int> > with all nodes for every line
         for each discretization type
  
  \return map with all nodes for each line
  */
  vector< vector<int> > getEleNodeNumberingLines(    
    const DRT::Element::DiscretizationType&     distype     ///< discretization type
    );

  /*!
  \brief Fills a vector< vector<int> > with all surfaces for every line
         for each discretization type
         
  \return map with surfaces adjacent to each line
  */                                                                    
  vector< vector<int> > getEleNodeNumbering_lines_surfaces(  
    const DRT::Element::DiscretizationType&     distype     ///< discretization type
    );

  /*!
  \brief Fills a vector< vector<int> > with all surfaces for every node
         for each discretization type
         
  \return map with all surfaces adjacent to each node
  */                                           
  vector< vector<int> > getEleNodeNumbering_nodes_surfaces(  
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Fills a vector< vector<int> > with reference coordinates for each node
         for each discretization type
         
  \return map of reference coordinates for all nodes
  */
  vector< vector< double > > getEleNodeNumbering_nodes_reference( 
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Fills an array with the surface ID s a point in reference coordinates
  \      is lying on for each discretization type
  */
  int getSurfaces(
    const blitz::Array<double,1>&               rst,        ///< point in reference coordinates
    int*                                        surfaces,   ///< surfaces the point is lying on
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Fills an array with coordinates in reference system of the cutter element
  \      according to the node Id for each discretization type
  */
  void getNodeCoordinates(
    const int                                   nodeId,     ///< node ID
    double*                                     coord,      ///< coordinates
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Fills an array with coordinates in reference system of the cutter element
  \      according to the line Id for each discretization type
  */
  void getLineCoordinates(
    const int                                   lineId,     ///< line ID
    const double                                lineCoord,  ///< line coordinate
    double*                                     coord,      ///< coordinates
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );
    
  /*!
  \brief    Returns the index of a  higher order node lying 
            between two given corner nodes of an edge of an higher order element
            for each discretization type
            
  */   
  int getHigherOrderIndex(
    const int                                   index1,     ///< index of corner node
    const int                                   index2,     ///< index of corner node
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief    Returns the dimension of an element
            
  */   
  int getDimension(
    const DRT::Element*   element     ///< element
    );
  
  /*!
  \brief    Returns the dimension of an element based on its discretization-type
            
  */   
  int getDimension(
		  const DRT::Element::DiscretizationType distype///< element-discretization type
    );
  
  /*!
  \brief    Returns the dimension of an element based on its discretization-type
            Template version
            
  */  
  template <DRT::Element::DiscretizationType DISTYPE>
  inline int getDimension()
  {
      int dim = 0;
      
      switch(DISTYPE)
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

  /*!
   * \brief check, whether NodeIds need rewinding due to wrong input from mesher,
   * based on negative Jacobian determinant.
   * \return boolean rewind or not rewind
   */
  bool checkRewinding3D(const DRT::Element* ele);
    
  /*!
  \brief    Returns the geometric center of the element in local coordinates
            
  */   
  std::vector<double> getLocalCenterPosition(
    const DRT::Element::DiscretizationType   distype     ///< shape of the element
    );

} // namespace UTILS
} // namespace DRT

#endif  // #ifdef DRT_UTILS_LOCAL_CONNECTIVITY_MATRICES_H
#endif  // #ifdef CCADISCRET
