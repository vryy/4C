/*----------------------------------------------------------------------*/
/*!
\file aero_tfsi_delaunay.cpp

\brief Helper class for coupled simulations (INCA - BACI)
       --> aero_tfsi.H & aero_tfsi.cpp
       The Delaunay triangulation makes use of software as stated below
       which is distributed under LGPL licence.

  Purpose:

    DTRIS2 constructs a Delaunay triangulation of 2D vertices.

  Discussion:

    The routine constructs the Delaunay triangulation of a set of 2D vertices
    using an incremental approach and diagonal edge swaps.  Vertices are
    first sorted in lexicographically increasing (X,Y) order, and
    then are inserted one at a time from outside the convex hull.

  Licensing:

    This code is distributed under the GNU LGPL license.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Modified:

    03 June 2009

  Author:

    Original FORTRAN77 version by Barry Joe.
    C++ version by John Burkardt.

  Reference:

    Barry Joe,
    GEOMPACK - a software package for the generation of meshes
    using geometric algorithms,
    Advances in Engineering Software,
    Volume 13, pages 325-331, 1991.

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              ghamm 12/11 |
 *----------------------------------------------------------------------*/
#ifndef AERO_TFSI_DELAUNAY_H
#define AERO_TFSI_DELAUNAY_H


/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include <map>
#include <vector>
#include "../linalg/linalg_fixedsizematrix.H"


namespace FS3I
{
  namespace UTILS
  {
    /// Helper class for coupled simulations INCA-BACI,
    /// Does Delaunay triangulation for needed for building mortar matrices
    /*!
      \author ghamm
      \date 11/12
     */

    // special delaunay triangulation which first projects aerocoords into plane of the structure
    // and then triangulates the points
	  void DelaunayTriangulation(std::map<int, LINALG::Matrix<3,1> > aerocoords, std::map<int, LINALG::Matrix<3,1> > structure, std::map<int, std::vector<int> >& tris);

	  // DIAEDG chooses a diagonal edge.
	  int diaedg(double x0, double y0, double x1, double y1, double x2,
				double y2, double x3, double y3);

	  // DTRIS2 constructs a Delaunay triangulation of 2D vertices.
		int dtris2(int point_num, int base, double point_xy[], int *tri_num,
				int tri_vert[], int tri_nabe[]);

    // I4_MODP returns the nonnegative remainder of I4 division.
    int i4_modp ( int i, int j );

    // I4_SIGN returns the sign of an I4.
    int i4_sign( int i );

    // I4_MAX returns the maximum of two I4's
    int i4_max ( int i1, int i2 );

    // I4_MIN returns the minimum of two I4's
    int i4_min ( int i1, int i2 );

		// I4_WRAP forces an I4 to lie between given limits by wrapping.
		int i4_wrap ( int ival, int ilo, int ihi );

    // I4VEC_MIN returns the value of the minimum element in an I4VEC.
    int i4vec_min ( int n, int a[] );

		// LRLINE determines where a point lies in relation to a directed line.
		int lrline ( double xu, double yu, double xv1, double yv1, double xv2, double yv2, double dv );

		// PERM_CHECK checks that a vector represents a permutation.
		bool perm_check ( int n, int p[], int base );

		// PERM_INVERSE inverts a permutation "in place".
		void perm_inverse ( int n, int p[] );

		// R8_MAX returns the maximum of two R8's.
		double r8_max ( double x, double y );

		// R82VEC_PERMUTE permutes an R82VEC in place.
		void r82vec_permute ( int n, int p[], int base, double a[] );

		// R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
		int *r82vec_sort_heap_index_a ( int n, int base, double a[] );

		// SWAPEC swaps diagonal edges until all triangles are Delaunay.
		int swapec ( int i, int *top, int *btri, int *bedg, int point_num,
					double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
					int stack[] );

		// VBEDG determines which boundary edges are visible to a point.
		void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num,
					int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg );
  }
}

#endif
