/*!----------------------------------------------------------------------
\file fluid3_xfem_nodemap.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fluid3_xfem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"


//Doku: reference on guide or pictures in here

// update!!
const int DRT::Elements::XFluid3::hex27_surfaces_[6][9] = {{0,  3,  2,  1, 11, 10,  9,  8, 20},
                                                           {4,  5,  6,  7, 12, 13, 14, 15, 21},
                                                           {0,  1,  5,  4,  8, 17, 12, 16, 22},
                                                           {2,  3,  7,  6, 10, 19, 14, 18, 23},
                                                           {0,  4,  7,  3, 16, 15, 19, 11, 24},
                                                           {1,  2,  6,  5,  9, 18, 13, 17, 25}};

// update!!
const int DRT::Elements::XFluid3::hex27_lines_[12][3] = {{0,  1,  8},
                                                         {1,  2,  9},
                                                         {2,  3, 10},
                                                         {3,  0, 11},
                                                         {0,  4, 16},
                                                         {1,  5, 17},
                                                         {2,  6, 18},
                                                         {3,  7, 19},
                                                         {4,  5, 12},
                                                         {5,  6, 13},
                                                         {6,  7, 14},
                                                         {7,  4, 15}};

// update!!
const int DRT::Elements::XFluid3::tet10_surfaces_[4][6] = {{0,  1,  8,  0,  1,  8},
                                                           {1,  2,  9,  0,  1,  8},
                                                           {2,  3, 10,  0,  1,  8},
                                                           {3,  0, 11,  0,  1,  8}};

// update!!
const int DRT::Elements::XFluid3::tet10_lines_[6][3] = {{0,  1,  8},
                                                        {1,  2,  9},
                                                        {2,  3, 10},
                                                        {3,  0, 11},
                                                        {0,  4, 16},
                                                        {1,  5, 17}};



/*
 *
 *                  
 *                  edge2     
 *             (line2)
 *           
 *              3     6     2
 *               X----o----X
 *               |         |
 *               |         |   
 *   edge3      7o    O    o5      edge1      
 *  (line3)      |    8    |      (line1)
 *               |         |
 *               X----o----X
 *              0     4     1
 *
 *                  edge0
 *                 (line0)
 *
 *
 *      X nodes for quad
 *      o addotional nodes for quad8
 *      O addotional nodes for full quadratic quad9
 */
const int DRT::Elements::XFluid3Surface::quad9_lines_[4][3] = {{0,  1,  4},
                                                               {1,  2,  5},
                                                               {2,  3,  6},
                                                               {3,  0,  7}};


/*
 *
 *               2
 *                X
 *                |\
 *                | \       edge1
 *   edge2       6o  o5    (line1)
 *  (line2)       |   \
 *                |    \
 *                X--o--X
 *               0   4   1
 *
 *             edge0
 *            (line0)
 *
 *
 *      X nodes for tri3
 *      o addotional nodes for tri6
 */
const int DRT::Elements::XFluid3Surface::tri6_lines_[3][3] = {{0,  1,  3},
                                                              {1,  2,  4},
                                                              {2,  0,  5}};

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
