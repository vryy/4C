/*!----------------------------------------------------------------------
\file drt_lib_fem_shapefunction.H

\brief Provide a node numbering scheme together with a set of shape functions

Provided are 1D, 2D and 3D shape functions

The surface mapping gives the node numbers such that the 2D shape functions can be used



<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef DRT_UTILS_FEM_SHAPEFUNCTION
#define DRT_UTILS_FEM_SHAPEFUNCTION

#include "drt_element.H"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
namespace Utils
{

  /*!
  \brief Nodal mappings describing the relation between volume, surface and line node numbering

  They should be used as the only reference for such relationships.

  The corresponding graphics and a detailed description can be found
  in the Baci guide in the Convention chapter.

  The numbering of lower order elements is included in the higher order element, such that
  e.g. the hex8 volume element uses only the first 8 nodes of the hex27 mapping
  */

  // 6 Surfaces of a Hex27 element with 9 nodes per surface
  const int     eleNodeNumbering_hex27_surfaces[6][9] = {{0,  3,  2,  1, 11, 10,  9,  8, 20},
                                                         {0,  1,  5,  4,  8, 13, 16, 12, 21},
                                                         {1,  2,  6,  5,  9, 14, 17, 13, 22},
                                                         {2,  3,  7,  6, 10, 15, 18, 14, 23},
                                                         {0,  4,  7,  3, 12, 19, 15, 11, 24},
                                                         {4,  5,  6,  7, 16, 17, 18, 19, 25}};

  // 12 Lines of a Hex27 element with 3 nodes per line
  const int     eleNodeNumbering_hex27_lines[12][3] = {{0,  1,  8},
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

  // 4 Surfaces of a Tet10 element with 6 nodes per surface
  const int     eleNodeNumbering_tet10_surfaces[4][6] = {{0,  1,  3,  4,  8,  7},
                                                         {1,  2,  3,  5,  9,  8},
                                                         {0,  3,  2,  7,  9,  6},
                                                         {0,  2,  1,  6,  5,  4}};

  // 6 Lines of a Tet10 element with 3 nodes per line
  const int     eleNodeNumbering_tet10_lines[6][3] = {{0,  1,  4},
                                                      {1,  2,  5},
                                                      {0,  2,  6},
                                                      {0,  3,  7},
                                                      {1,  3,  8},
                                                      {2,  3,  9}};

   //Triangular surfaces of weg15 with 7 nodes per surface
  const int     eleNodeNumbering_weg15_trisurfaces[2][6]={{0, 1 ,2, 6, 7, 8},
                                                          {3, 4, 5, 12, 13, 14}};

  //Rectangular surfaces of weg15 with 8 nodes per surface
  const int     eleNodeNumbering_weg15_quadsurfaces[3][8]={{0, 1, 4, 3, 6, 10, 12, 9},
                                                           {1, 2, 5, 4, 7, 11, 13, 10},
                                                           {0, 3, 5, 2, 9, 14, 11, 8}};

  //9 lines of a Weg15 element with 3 nodes per line
  const int     eleNodeNumbering_weg15_lines[9][3]={{0, 1, 6},
                                                    {1, 2, 7},
                                                    {2, 0, 8},
                                                    {3, 4, 12},
                                                    {4, 5, 13},
                                                    {5, 3, 14},
                                                    {1, 4, 10},
                                                    {2, 5, 11},
                                                    {0, 3, 9}};

  // 4 Lines of a Quad9 element with 3 nodes per line
  const int     eleNodeNumbering_quad9_lines[4][3] = {{0,  1,  4},
                                                      {1,  2,  5},
                                                      {2,  3,  6},
                                                      {0,  3,  7}};

  // 3 Lines of a Tri6 element with 3 nodes per line
  const int     eleNodeNumbering_tri6_lines[3][3] = {{0,  1,  3},
                                                     {1,  2,  4},
                                                     {2,  0,  5}};


  // for each of the 12 lines, tell me the 2 neighbouring faces
  const int	eleNodeNumbering_hex27_lines_surfaces[12][2] = { {0,  1},
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

  /*!
  \brief Fill a Epetra_SerialDenseVector with with 3D shape function

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   s        DOUBLE   (i)    coordinate
  \param   t        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_3D( Epetra_SerialDenseVector&                  funct,
                          const double&                              r,
                          const double&                              s,
                          const double&                              t,
                          const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseMatrix with with 3D first shape function derivatives

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   s        DOUBLE   (i)    coordinate
  \param   t        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_3D_deriv1( Epetra_SerialDenseMatrix&                  deriv1,
                                 const double&                              r,
                                 const double&                              s,
                                 const double&                              t,
                                 const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseMatrix with with 3D second shape function derivatives

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   s        DOUBLE   (i)    coordinate
  \param   t        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_3D_deriv2( Epetra_SerialDenseMatrix&                  deriv2,
                                 const double&                              r,
                                 const double&                              s,
                                 const double&                              t,
                                 const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseVector with 2D shape function

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   s        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_2D( Epetra_SerialDenseVector&                  funct,
                          const double&                              r,
                          const double&                              s,
                          const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseVector with first 2D shape function derivative

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   s        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_2D_deriv1( Epetra_SerialDenseMatrix&                  deriv1,
                                 const double&                              r,
                                 const double&                              s,
                                 const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseVector with second 2D shape function derivative

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   s        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_2D_deriv2( Epetra_SerialDenseMatrix&                  deriv2,
                                 const double&                              r,
                                 const double&                              s,
                                 const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseVector with 1D shape function

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_1D( Epetra_SerialDenseVector&                  funct,
                          const double&                              r,
                          const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseVector with first 1D shape function derivative

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_1D_deriv1( Epetra_SerialDenseMatrix&                  deriv1,
                                 const double&                              r,
                                 const DRT::Element::DiscretizationType&    distype);

  /*!
  \brief Fill a Epetra_SerialDenseVector with second 1D shape function derivative

  \param  *funct    DOUBLE   (o)    shape functions
  \param   r        DOUBLE   (i)    coordinate
  \param   distype           (i)    distuinguish between DiscretizationType
  */
  void shape_function_1D_deriv2( Epetra_SerialDenseMatrix&                  deriv2,
                                 const double&                              r,
                                 const DRT::Element::DiscretizationType&    distype);

} // namespace Utils
} // namespace DRT

#endif  // #ifdef DRT_UTILS_FEM_SHAPEFUNCTION
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET


