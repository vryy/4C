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


//!!! Doku: reference on guide, chapter conventions !!!
const int DRT::Elements::XFluid3::hex27_surfaces_[6][9] = {{0,  3,  2,  1, 11, 10,  9,  8, 20},
                                                           {0,  1,  5,  4,  8, 13, 16, 12, 21},
                                                           {1,  2,  6,  5,  9, 14, 17, 13, 22},
                                                           {2,  3,  7,  6, 10, 15, 18, 14, 23},
                                                           {0,  4,  7,  3, 12, 19, 15, 11, 24},
                                                           {4,  5,  6,  7, 16, 17, 18, 19, 25}};

const int DRT::Elements::XFluid3::hex27_lines_[12][3] = {{0,  1,  8},
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

const int DRT::Elements::XFluid3::tet10_surfaces_[4][6] = {{0,  1,  3,  4,  8,  7},
                                                           {1,  2,  3,  5,  9,  8},
                                                           {0,  3,  2,  7,  9,  6},
                                                           {0,  2,  1,  6,  5,  4}};

const int DRT::Elements::XFluid3::tet10_lines_[6][3] = {{0,  1,  4},
                                                        {1,  2,  5},
                                                        {0,  2,  6},
                                                        {0,  3,  7},
                                                        {1,  3,  8},
                                                        {2,  3,  9}};

const int DRT::Elements::XFluid3Surface::quad9_lines_[4][3] = {{0,  1,  4},
                                                               {1,  2,  5},
                                                               {2,  3,  6},
                                                               {0,  3,  7}};

const int DRT::Elements::XFluid3Surface::tri6_lines_[3][3] = {{0,  1,  3},
                                                              {1,  2,  4},
                                                              {2,  0,  5}};

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
