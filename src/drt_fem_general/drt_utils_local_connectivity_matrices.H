/*!----------------------------------------------------------------------
\file drt_utils_local_connectivity_matrices.H

\brief Provide a node numbering scheme together with a set of shape functions

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

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DRT_UTILS_LOCAL_CONNECTIVITY_MATRICES_H
#define DRT_UTILS_LOCAL_CONNECTIVITY_MATRICES_H

#include "../drt_lib/drt_element.H"

namespace DRT
{
namespace UTILS
{
  //! 2 surfaces on a Hex16 element with 8 nodes per surface
  const int eleNodeNumbering_hex16_surfaces_q8[2][8] = {
          {0,  3,  2,  1,  7,  6,  5,  4},
          {8,  9, 10, 11, 12, 13, 14, 15}};
  //! 4 surfaces on a Hex16 element with 6 nodes per surface
  const int eleNodeNumbering_hex16_surfaces_q6[4][6] = {
          {0,  1,  4,  8,  9, 12},
          {1,  2,  5,  9, 10, 13},
          {2,  3,  6, 10, 11, 14},
          {3,  0,  7, 11,  8, 15}};

  //! 8 Lines of a Hex16 element with 3 nodes per line
  const int eleNodeNumbering_hex16_lines_quad[8][3] = {
          { 0,  1,  4},
          { 1,  2,  5},
          { 2,  3,  6},
          { 0,  3,  7},
          { 8,  9, 12},
          { 9, 10, 13},
          {10, 11, 14},
          { 8, 11, 15}};
  //! 4 Lines in a Hex16 element with 2 nodes per line
  const int eleNodeNumbering_hex16_lines_lin[4][2] = {
          {0, 8},
          {1, 9},
          {2,10},
          {3,11}};

  //! 2 surfaces on a Hex18 element with 9 nodes per surface
  const int eleNodeNumbering_hex18_surfaces_q9[2][9] = {
          {0,  3,  2,  1,  7,  6,  5,  4,  8},
          {9, 10, 11, 12, 13, 14, 15, 16, 17}};
  //! 4 surfaces on a Hex18 element with 6 nodes per surface
  const int eleNodeNumbering_hex18_surfaces_q6[4][6] = {
          {0,  1,  4,  9, 10, 13},
          {1,  2,  5, 10, 11, 14},
          {2,  3,  6, 11, 12, 15},
          {3,  0,  7, 12,  9, 16}};

  //! 8 Lines of a Hex18 element with 3 nodes per line
  const int eleNodeNumbering_hex18_lines_quad[8][3] = {
          { 0,  1,  4},
          { 1,  2,  5},
          { 2,  3,  6},
          { 0,  3,  7},
          { 9, 10, 13},
          {10, 11, 14},
          {11, 12, 15},
          { 9, 12, 16}};
  //! 4 Lines in a Hex18 element with 2 nodes per line
  const int eleNodeNumbering_hex18_lines_lin[4][2] = {
          {0, 9},
          {1,10},
          {2,11},
          {3,12}};

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
          { 2, 1, 0, 7, 6, 8},
          { 5, 3, 4, 14, 12, 13}};

  //! Rectangular surfaces of wedge15 with 8 nodes per surface
  const int eleNodeNumbering_wedge15_quadsurfaces[3][8] = {
          {0, 1, 4, 3, 6, 10, 12,  9},
          {1, 2, 5, 4, 7, 11, 13, 10},
          {2, 0, 3, 5, 8, 9, 14, 11}};

  //! 9 lines of a wedge15 element with 3 nodes per line
  const int eleNodeNumbering_wedge15_lines[9][3] = {
          {0, 1, 6},
          {1, 2, 7},
          {2, 0, 8},
          {3, 4, 12},
          {4, 5, 13},
          {5, 3, 14},
          {0, 3, 9},
          {1, 4, 10},
          {2, 5, 11}};

  //! Rectangular surface of pyramid5
  const int eleNodeNumbering_pyramid5_quadsurfaces[1][4]= {
          {0, 3, 2, 1}};

  //! Triangular surfaces of pyramid5
  const int eleNodeNumbering_pyramid5_trisurfaces[4][3]={
    {0,1,4},
    {1,2,4},
    {2,3,4},
    {0,4,3}};

  //! 9 lines of a wedge15 element with 3 nodes per line
  const int eleNodeNumbering_pyramid5_lines[8][2] = {
          {0, 1},
          {1, 2},
          {2, 3},
          {3, 0},
          {0, 4},
          {1, 4},
          {2, 4},
          {3, 4}};


  //! 4 Lines of a Quad9 element with 3 nodes per line
  const int eleNodeNumbering_quad9_lines[4][3] = {
          {0,  1,  4},
          {1,  2,  5},
          {2,  3,  6},
          {3,  0,  7}};

  //! 4 Lines of a nurbs4 element with 2 control points per line
  //! we apply the nurbs specific cartesian numbering
  //!
  //! \verbatim
  //!
  //!          2   (2)   3
  //!           +-------+
  //!           |       |
  //!           |       |(1)
  //!        (3)|       |
  //!           |       |
  //!           +-------+
  //!          0   (0)   1
  //!
  //! \endverbatim

  const int eleNodeNumbering_nurbs4_lines[4][2] = {
          {0,  1},
          {1,  2},
          {2,  3},
          {3,  0}};


  //! 4 Lines of a nurbs9 element with 3 control points per line
  //! we apply the nurbs specific cartesian numbering.
  //!
  //! \verbatim
  //!
  //!            (2) -->
  //!          6    7    8
  //!           +---+---+
  //!           |       | (1)
  //!       (3) |       |
  //!           +   +   +  ^
  //!        ^ 3|   4   |5 |
  //!        |  |       |
  //!           +---+---+
  //!          0    1    2
  //!            (0) -->
  //!
  //! \endverbatim

  const int eleNodeNumbering_nurbs9_lines[4][3] = {
          {0,  1,  2},
          {2,  5,  8},
          {6,  7,  8},
          {0,  3,  6}};

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


  //! for each of the 8 corner nodes, tell me the 3 neighbouring lines
  const int eleNodeNumbering_hex27_nodes_lines[8][3] = { {0,  3,  4},
                                                         {0,  1,  5},
                                                         {1,  2,  6},
                                                         {2,  3,  7},
                                                         {4,  8,  11},
                                                         {5,  8,  9},
                                                         {6,  9,  10},
                                                         {7,  10, 11}};


   //! for each of the 4 corner nodes, tell me the 3 neighbouring faces
  const int eleNodeNumbering_tet10_nodes_surfaces[4][3] =  { {0,  2,  3},
                                                             {0,  1,  3},
                                                             {1,  2,  3},
                                                             {0,  1,  2}};


  //! for each of the 4 corner nodes, tell me the 3 neighbouring lines
 const int eleNodeNumbering_tet10_nodes_lines[4][3] =  { {0,  2,  3},
                                                         {0,  1,  4},
                                                         {1,  2,  5},
                                                         {3,  4,  5}};

  //! for each node the reference coordinates are stored for quad9
  const double eleNodeNumbering_quad9_nodes_reference[9][2] = { {-1.0, -1.0},
                                                                { 1.0, -1.0},
                                                                { 1.0,  1.0},
                                                                {-1.0,  1.0},

                                                                { 0.0, -1.0},
                                                                { 1.0,  0.0},
                                                                { 0.0,  1.0},
                                                                {-1.0,  0.0},

                                                                { 0.0,  0.0} };


  //! 6 Surfaces of a nurbs8 element with 4 control points per
  //! surface we apply the nurbs specific cartesian numbering,
  //! only interpolated surfaces are valid!
  //!
  //! \verbatim
  //!                         v
  //!                        /
  //!                       /
  //!          X-------------------X
  //!         /|6         /       /|7
  //!  w     / |         /       / |
  //!  ^    /  |        /       /  |
  //!  |   /   |               /   |
  //!  |  /    |              /    |
  //!  | /     |             /     |
  //!   /      |            /      |
  //!  X-------------------X       |
  //!  |4      |           |5      |
  //!  |       |           |       |
  //!  |       |           |       |
  //!  |       X-----------|-------X
  //!  |      / 2          |      / 3
  //!  |     /             |     /
  //!  |    /              |    /
  //!  |   /               |   /
  //!  |  /                |  /
  //!  | /                 | /
  //!  |/                  |/
  //!  X-------------------X ----->u
  //!   0                   1
  //!
  //! \endverbatim

  const int eleNodeNumbering_nurbs8_surfaces[6][4] = {
          {0,1,2,3},
          {4,5,6,7},
          {0,1,4,5},
          {2,3,6,7},
          {0,2,4,6},
          {1,3,5,7}};


  //! 6 Surfaces of a nurbs27 element with 9 control points per surface
  //! we apply the nurbs specific cartesian numbering, only interpolated surfaces are valid!
  //!
  //! \verbatim
  //!                          v
  //!                         /
  //!                        /
  //!          X---------X---------X
  //!         /|24      /|25      /|26
  //!  w     / |       / |       / |
  //!  ^    /  |      /  |      /  |
  //!  |   X---------X---------X   |
  //!  |  /|21 |    /|22 |    /|23 |
  //!  | / |   X---/-|---X---/-|---X
  //!   /  |  /|15/  |  /|16/  |  /|17
  //!  X---------X---------X   | / |
  //!  |18 |/  | |19 |/  | |20 |/  |
  //!  |   X-----|---X-----|---X   |
  //!  |  /|12 | |  /|13 | |  /|14 |
  //!  | / |   X-|-/-|---X-|-/-|---X
  //!  |/  |  /6 |/  |  /7 |/  |  /8
  //!  X---------X---------X   | /
  //!  |9  |/    |10 |/    |11 |/
  //!  |   X-----|---X-----|---X
  //!  |  /3     |  /4     |  /5
  //!  | /       | /       | /
  //!  |/        |/        |/
  //!  X---------X---------X ----->u
  //!   0         1         2
  //!
  //! \endverbatim

  const int eleNodeNumbering_nurbs27_surfaces[6][9] = {
          { 0, 1, 2, 3, 4, 5, 6, 7, 8},
          {18,19,20,21,22,23,24,25,26},
          { 0, 1, 2, 9,10,11,18,19,20},
          { 6, 7, 8,15,16,17,24,25,26},
          { 2, 5, 8,11,14,17,20,23,26},
          { 0, 3, 6, 9,12,15,18,21,24}};

  //! 12 Lines of a Nurbs27 element with 3 control points per line
  //!                          v
  //!                         /
  //!                   10   /
  //!          X---------X---------X
  //!         /|                  /|
  //!  w     / |                 / |
  //!  ^    /  |                /  |
  //!  | 11X   |               X9  |
  //!  |  /    |              /    |
  //!  | /     X7            /     X6
  //!   /      |            /      |
  //!  X---------X---------X       |
  //!  |       | 8         |       |
  //!  |       |           |       |
  //!  |       |         2 |       |
  //!  |       X---------X-|-------X
  //!  |      /            |      /
  //! 4X     /             X5    /
  //!  |    /              |    /
  //!  |  3X               |   X1
  //!  |  /                |  /
  //!  | /                 | /
  //!  |/                  |/
  //!  X---------X---------X ----->u
  //!            0
  //!
  const int eleNodeNumbering_nurbs27_lines[12][3] = {
         { 0,  1,  2},
         { 2,  5,  8},
         { 6,  7,  8},
         { 0,  3,  6},
         { 0,  9, 18},
         { 2, 11, 20},
         { 8, 17, 26},
         { 6, 15, 24},
         {18, 19, 20},
         {20, 23, 26},
         {24, 25, 26},
         {18, 21, 24}};


  //! for each node the reference coordinates are stored for tri6
  const double eleNodeNumbering_tri6_nodes_reference[6][2] = { {0.0, 0.0},
                                                               {1.0, 0.0},
                                                               {0.0, 1.0},

                                                               {0.5, 0.0},
                                                               {0.5, 0.5},
                                                               {0.0, 0.5} };


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
  const double eleNodeNumbering_line3_nodes_reference[3][1] = {     {-1.0},
                                                                    { 1.0},
                                                                    { 0.0}};
  /*!
 \brief Returns the number of nodes
        for each discretization type

 \return number of nodes of the discretization type
 */
 int getNumberOfElementNodes(
   const DRT::Element::DiscretizationType&     distype
   );


  //! Template Meta Programming version of switch over discretization type
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToNumNodePerEle {};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::point1>   {static const int numNodePerElement =  1;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::line2>    {static const int numNodePerElement =  2;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::line3>    {static const int numNodePerElement =  3;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::nurbs2>   {static const int numNodePerElement =  2;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::nurbs3>   {static const int numNodePerElement =  3;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::quad4>    {static const int numNodePerElement =  4;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::quad8>    {static const int numNodePerElement =  8;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::quad9>    {static const int numNodePerElement =  9;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::nurbs4>   {static const int numNodePerElement =  4;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::nurbs9>   {static const int numNodePerElement =  9;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::tri3>     {static const int numNodePerElement =  3;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::tri6>     {static const int numNodePerElement =  6;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::hex8>     {static const int numNodePerElement =  8;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::hex16>    {static const int numNodePerElement = 16;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::hex18>    {static const int numNodePerElement = 18;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::hex20>    {static const int numNodePerElement = 20;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::hex27>    {static const int numNodePerElement = 27;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::nurbs8>   {static const int numNodePerElement =  8;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::nurbs27>  {static const int numNodePerElement = 27;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::tet4>     {static const int numNodePerElement =  4;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::tet10>    {static const int numNodePerElement = 10;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::wedge6>   {static const int numNodePerElement =  6;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::wedge15>  {static const int numNodePerElement = 15;};
  template<> struct DisTypeToNumNodePerEle<DRT::Element::pyramid5> {static const int numNodePerElement =  5;};


  /*!
  \brief Returns the number of corner nodes
         for each discretization type

  \return number of corner nodes of the discretization type
  */
  int getNumberOfElementCornerNodes(
    const DRT::Element::DiscretizationType&     distype
    );


  /*!
  \brief  Returns the number of corner nodes
          for each surface element of a given volume discretization type

  \return   number of corner nodes for each surface element of a given volume
             of the discretization type
   */
  std::vector<int> getNumberOfSurfaceElementCornerNodes(
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
 \brief Returns the number of Volumes
        for each discretization type
   */
  int getNumberOfElementVolumes(
      const DRT::Element::DiscretizationType&     distype);

  /*!
 \brief Returns the number of faces
        for each discretization type (as opposed to Lines/Surfaces, this is always dim-1)
   */
  int getNumberOfElementFaces(
      const DRT::Element::DiscretizationType&     distype);

  /*!
  \brief Returns the shape type of a face of a given element
   */
  DRT::Element::DiscretizationType getEleFaceShapeType (
      const DRT::Element::DiscretizationType&     distype,     ///< discretization type
      const unsigned int face = 0                              ///< which face?
      );

  /*!
  \brief Fills a vector< std::vector<int> > with all nodes for every face
         (surface in 3D, line in 2D) for each discretization type

  \return map with all nodes for each face
   */
  std::vector< std::vector<int> > getEleNodeNumberingFaces(
      const DRT::Element::DiscretizationType&     distype     ///< discretization type
  );

  /*!
  \brief Fills a vector< std::vector<int> > with all nodes for every surface
         for each discretization type

  \return map with all nodes for each surface
   */
  std::vector< std::vector<int> > getEleNodeNumberingSurfaces(
      const DRT::Element::DiscretizationType&     distype     ///< discretization type
  );

  /*!
  \brief Fills a vector< std::vector<int> > with all nodes for every line
         for each discretization type

  \return map with all nodes for each line
   */
  std::vector< std::vector<int> > getEleNodeNumberingLines(
      const DRT::Element::DiscretizationType&     distype     ///< discretization type
  );

  /*!
  \brief Fills a vector< std::vector<int> > with all surfaces for every line
         for each discretization type

  \return map with surfaces adjacent to each line
  */
  std::vector< std::vector<int> > getEleNodeNumbering_lines_surfaces(
    const DRT::Element::DiscretizationType&     distype     ///< discretization type
    );

  /*!
  \brief Fills a vector< std::vector<int> > with all lines for every node
         for each discretization type

  \return map with all lines adjacent to each node
  */
  std::vector< std::vector<int> > getEleNodeNumbering_nodes_lines(
      const DRT::Element::DiscretizationType      distype
      );

  /*!
  \brief Fills a vector< std::vector<int> > with all surfaces for every node
         for each discretization type

  \return map with all surfaces adjacent to each node
  */
  std::vector< std::vector<int> > getEleNodeNumbering_nodes_surfaces(
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Fills a LINALG::SerialDenseMatrix with parameter space coordinates for each node
         for each discretization type

  \return map of parameter space coordinates for all nodes
  */
  LINALG::SerialDenseMatrix getEleNodeNumbering_nodes_paramspace(
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Returns a vector with the surface ID s a point in 3d reference coordinates
  \      is lying on for each discretization type
  */
  std::vector<int> getSurfaces(
    const LINALG::Matrix<3,1>&                  rst,        ///< point in reference coordinates
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );


  /*!
   \brief Returns a vector with the line ID s a point in 3d reference coordinates
   \      is lying on for each discretization type
   */
  std::vector<int> getLines(
      const LINALG::Matrix<3,1>&                  rst,      ///< point in reference coordinates
      const DRT::Element::DiscretizationType      distype   ///< discretization type
      );

   /*!
   \brief Returns the node ID s a point in 3d reference coordinates
   \      is lying on for each discretization type
   */
  int getNode(
      const LINALG::Matrix<3,1>&                  rst,      ///< point in reference coordinates
      const DRT::Element::DiscretizationType      distype   ///< discretization type
      );


  /*!
  \brief Returns a vector with coordinates in reference system of the cutter element
  \      according to the node Id for each discretization type
  */
  LINALG::Matrix<3,1> getNodeCoordinates(
    const int                                   nodeId,     ///< node Id
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief Returns a with coordinates in reference system of the cutter element
  \      according to the line Id for each discretization type
  */
  LINALG::Matrix<3,1> getLineCoordinates(
    const int                                   lineId,     ///< line ID
    const double                                lineCoord,  ///< line coordinate
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
  \brief    Returns the indices of the element corner nodes lying
            adjacent to a specified higher order (edge) node index
            for each discretization type

  */
  void getCornerNodeIndices(
    int&                                        index1,     ///< index of corner node
    int&                                        index2,     ///< index of corner node
    const int&                                  hoindex,    ///< index of edge node
    const DRT::Element::DiscretizationType      distype     ///< discretization type
    );

  /*!
  \brief    Returns the dimension of an element based on its discretization-type

  */
  int getDimension(
      const DRT::Element::DiscretizationType distype  ///< element-discretization type
    );

  /*!
  \brief    Returns the dimension of an element based on its discretization-type
            Template version
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToDim {};
  template<> struct DisTypeToDim<DRT::Element::point1>   {static const int dim = 0;};
  template<> struct DisTypeToDim<DRT::Element::line2>    {static const int dim = 1;};
  template<> struct DisTypeToDim<DRT::Element::line3>    {static const int dim = 1;};
  template<> struct DisTypeToDim<DRT::Element::nurbs2>   {static const int dim = 1;};
  template<> struct DisTypeToDim<DRT::Element::nurbs3>   {static const int dim = 1;};
  template<> struct DisTypeToDim<DRT::Element::quad4>    {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::quad6>    {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::quad8>    {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::quad9>    {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::nurbs4>   {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::nurbs9>   {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::tri3>     {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::tri6>     {static const int dim = 2;};
  template<> struct DisTypeToDim<DRT::Element::hex8>     {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::hex16>    {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::hex18>    {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::hex20>    {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::hex27>    {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::tet4>     {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::tet10>    {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::wedge6>   {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::wedge15>  {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::pyramid5> {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::nurbs8>   {static const int dim = 3;};
  template<> struct DisTypeToDim<DRT::Element::nurbs27>  {static const int dim = 3;};

  /*!
  \brief    Returns the degree of an element based on its discretization-type
            Template version
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToDegree {};
  template<> struct DisTypeToDegree<DRT::Element::point1>   {static const int degree =0;};
  template<> struct DisTypeToDegree<DRT::Element::line2>    {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::line3>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::nurbs2>   {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::nurbs3>   {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::quad4>    {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::quad6>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::quad8>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::quad9>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::nurbs4>   {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::nurbs9>   {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::tri3>     {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::tri6>     {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::hex8>     {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::hex20>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::hex27>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::tet4>     {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::tet10>    {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::wedge6>   {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::wedge15>  {static const int degree =2;};
  template<> struct DisTypeToDegree<DRT::Element::pyramid5> {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::nurbs8>   {static const int degree =1;};
  template<> struct DisTypeToDegree<DRT::Element::nurbs27>  {static const int degree =2;};

  /*!
  \brief    Returns the number of independent components for 2nd derivatives per discretization-type
            Template version

            Number of components necessary to store second derivatives:
            in general:
              (n*(n+1))/2 components for nsd=n;
            in practice:
              1 component  for nsd=1:  (N,xx)
              3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
              6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToNumDeriv2 {};
  template<> struct DisTypeToNumDeriv2<DRT::Element::point1>   {static const int numderiv2 = 0;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::line2>    {static const int numderiv2 = 1;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::line3>    {static const int numderiv2 = 1;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::nurbs2>   {static const int numderiv2 = 1;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::nurbs3>   {static const int numderiv2 = 1;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::quad4>    {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::quad6>    {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::quad8>    {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::quad9>    {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::nurbs4>   {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::nurbs9>   {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::tri3>     {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::tri6>     {static const int numderiv2 = 3;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::hex8>     {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::hex16>    {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::hex18>    {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::hex20>    {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::hex27>    {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::tet4>     {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::tet10>    {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::wedge6>   {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::wedge15>  {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::pyramid5> {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::nurbs8>   {static const int numderiv2 = 6;};
  template<> struct DisTypeToNumDeriv2<DRT::Element::nurbs27>  {static const int numderiv2 = 6;};


  /*!
  \brief    Returns the shape of face elements of a given element at compile time
  (invalid for wedge and pyramid as it depends on the face number)
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToFaceShapeType {static const DRT::Element::DiscretizationType shape = DRT::Element::max_distype;};
  //template<> struct DisTypeToFaceShapeType<DRT::Element::point1>   {}; // invalid
  template<> struct DisTypeToFaceShapeType<DRT::Element::line2>    {static const DRT::Element::DiscretizationType shape = DRT::Element::point1;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::line3>    {static const DRT::Element::DiscretizationType shape = DRT::Element::point1;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::nurbs2>   {static const DRT::Element::DiscretizationType shape = DRT::Element::point1;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::nurbs3>   {static const DRT::Element::DiscretizationType shape = DRT::Element::point1;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::quad4>    {static const DRT::Element::DiscretizationType shape = DRT::Element::line2;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::quad6>    {static const DRT::Element::DiscretizationType shape = DRT::Element::line3;}; // invalid
  template<> struct DisTypeToFaceShapeType<DRT::Element::quad8>    {static const DRT::Element::DiscretizationType shape = DRT::Element::line3;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::quad9>    {static const DRT::Element::DiscretizationType shape = DRT::Element::line3;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::nurbs4>   {static const DRT::Element::DiscretizationType shape = DRT::Element::nurbs2;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::nurbs9>   {static const DRT::Element::DiscretizationType shape = DRT::Element::nurbs3;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::tri3>     {static const DRT::Element::DiscretizationType shape = DRT::Element::line2;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::tri6>     {static const DRT::Element::DiscretizationType shape = DRT::Element::line3;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::hex8>     {static const DRT::Element::DiscretizationType shape = DRT::Element::quad4;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::hex16>    {static const DRT::Element::DiscretizationType shape = DRT::Element::quad8;}; // invalid
  template<> struct DisTypeToFaceShapeType<DRT::Element::hex18>    {static const DRT::Element::DiscretizationType shape = DRT::Element::quad9;}; // invalid
  template<> struct DisTypeToFaceShapeType<DRT::Element::hex20>    {static const DRT::Element::DiscretizationType shape = DRT::Element::quad8;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::hex27>    {static const DRT::Element::DiscretizationType shape = DRT::Element::quad9;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::nurbs8>   {static const DRT::Element::DiscretizationType shape = DRT::Element::nurbs4;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::nurbs27>  {static const DRT::Element::DiscretizationType shape = DRT::Element::nurbs9;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::tet4>     {static const DRT::Element::DiscretizationType shape = DRT::Element::tri3;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::tet10>    {static const DRT::Element::DiscretizationType shape = DRT::Element::tri6;};
  template<> struct DisTypeToFaceShapeType<DRT::Element::wedge6>   {static const DRT::Element::DiscretizationType shape = DRT::Element::tri3;}; //invalid
  template<> struct DisTypeToFaceShapeType<DRT::Element::wedge15>  {static const DRT::Element::DiscretizationType shape = DRT::Element::tri6;}; //invalid
  template<> struct DisTypeToFaceShapeType<DRT::Element::pyramid5> {static const DRT::Element::DiscretizationType shape = DRT::Element::tri3;}; //invalid



  /*!
  \brief    Returns the shape of face elements of a given element
  (invalid for wedge and pyramid as it depends on the face number)
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToNumFaces {static const int numFaces = 0;};
  //template<> struct DisTypeToNumFaces<DRT::Element::point1>   {}; // invalid
  template<> struct DisTypeToNumFaces<DRT::Element::line2>    {static const int numFaces = 2;};
  template<> struct DisTypeToNumFaces<DRT::Element::line3>    {static const int numFaces = 2;};
  template<> struct DisTypeToNumFaces<DRT::Element::nurbs2>   {static const int numFaces = 2;};
  template<> struct DisTypeToNumFaces<DRT::Element::nurbs3>   {static const int numFaces = 2;};
  template<> struct DisTypeToNumFaces<DRT::Element::quad4>    {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::quad6>    {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::quad8>    {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::quad9>    {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::nurbs4>   {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::nurbs9>   {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::tri3>     {static const int numFaces = 3;};
  template<> struct DisTypeToNumFaces<DRT::Element::tri6>     {static const int numFaces = 3;};
  template<> struct DisTypeToNumFaces<DRT::Element::hex8>     {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::hex16>    {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::hex18>    {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::hex20>    {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::hex27>    {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::nurbs8>   {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::nurbs27>  {static const int numFaces = 6;};
  template<> struct DisTypeToNumFaces<DRT::Element::tet4>     {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::tet10>    {static const int numFaces = 4;};
  template<> struct DisTypeToNumFaces<DRT::Element::wedge6>   {static const int numFaces = 5;};
  template<> struct DisTypeToNumFaces<DRT::Element::wedge15>  {static const int numFaces = 5;};
  template<> struct DisTypeToNumFaces<DRT::Element::pyramid5> {static const int numFaces = 5;};


  /*!
  \brief    Returns the number of nodes on face elements
  (invalid for wedge and pyramid as it depends on the face number)
  */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToNumNodePerFace { static const int numNodePerFace = DisTypeToNumNodePerEle<DisTypeToFaceShapeType<DISTYPE>::shape>::numNodePerElement; };


  /*!
    \brief    Returns the order of an element edges based on its discretization-type
              Template version
   */
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToEdgeOrder {};
  template<> struct DisTypeToEdgeOrder<DRT::Element::point1>   {static const int order = 0;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::line2>    {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::line3>    {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::nurbs2>   {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::nurbs3>   {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::quad4>    {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::quad8>    {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::quad9>    {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::nurbs4>   {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::nurbs9>   {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::tri3>     {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::tri6>     {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::hex8>     {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::hex20>    {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::hex27>    {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::tet4>     {static const int order = 1;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::tet10>    {static const int order = 2;};
  template<> struct DisTypeToEdgeOrder<DRT::Element::pyramid5> {static const int order = 1;};


  /*!
    \brief    returns the order of the element-edges
   */
  int getOrder(
      const DRT::Element::DiscretizationType distype///< element-discretization type
  );

  /*!
    \brief    returns the order of the element
   */
  int getDegree(
      const DRT::Element::DiscretizationType distype///< element-discretization type
  );

  /*!
    \brief    returns local node number in the parent element
              from a given face number and local node number in the face element
   */
  int getParentNodeNumberFromFaceNodeNumber(
      const DRT::Element::DiscretizationType parent_distype,
      const int faceId, const int faceNodeId);

  /*!
   * Returns the geometric center of the element in local coordinates
   *
   * \tparam   dim  dimension of the element space (e.g. dim=2 for quad4)
   */
  template<int dim>
  LINALG::Matrix<dim,1> getLocalCenterPosition(
    const DRT::Element::DiscretizationType   distype     ///< shape of the element
    )
  {
    LINALG::Matrix<dim,1> pos;
    switch(distype)
    {
        case DRT::Element::line2 :  case DRT::Element::line3 :
        {
            pos = 0.0;
            break;
        }
        case DRT::Element::quad4 : case DRT::Element::quad6 : case DRT::Element::quad8 : case DRT::Element::quad9 :
        {
            pos = 0.0;
            break;
        }
        case DRT::Element::tri3 : case DRT::Element::tri6 :
        {
            pos = 1.0/3.0;
            break;
        }
        case DRT::Element::hex8 : case DRT::Element::hex16 : case DRT::Element::hex18 : case DRT::Element::hex20 : case DRT::Element::hex27 :
        {
            pos = 0.0;
            break;
        }
        case DRT::Element::tet4 : case DRT::Element::tet10 :
        {
            pos = 1.0/4.0;
            break;
        }
        case DRT::Element::wedge6 : case DRT::Element::wedge15 :
        {
            pos = 1.0/3.0;
            pos( 2, 0 ) = 0;
            break;
        }
        case DRT::Element::pyramid5 :
        {
            pos = 0.0;
            pos( 2, 0 ) = 0.25;
            break;
        }
        default:
            dserror("discretization type %s not yet implemented", (DRT::DistypeToString(distype)).c_str()); break;
    }
    return pos;
  }

  /*!
   *  returns the element size in local coordinates
   */
  double getSizeInLocalCoordinates(
      const DRT::Element::DiscretizationType     distype);

  /*!
   *  returns the shape of a given boundary element (is used in various element classes
   *  for implementing the Shape() method of boundary elements)
   */
  DRT::Element::DiscretizationType getShapeOfBoundaryElement(
      const int nen, ///< number of nodes of the boundary element
      const DRT::Element::DiscretizationType parentshape ///< shape of parent element
  );

  /*!
   * returns the power of a number as a compile time constant
   */
  template <int a, int N>
  struct FixedPower
  {
    static const int value = a * FixedPower<a,N-1>::value;
  };

  /*!
   * Base case for fixed power: A number raised to power 0 is 1.
   */
  template <int a>
  struct FixedPower<a,0>
  {
    static const int value = 1;
  };


} // namespace UTILS
} // namespace DRT

#endif  // #ifdef DRT_UTILS_LOCAL_CONNECTIVITY_MATRICES_H
