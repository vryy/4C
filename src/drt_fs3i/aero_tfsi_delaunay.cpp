/*----------------------------------------------------------------------*/
/*!
 \file aero_tfsi_delaunay.cpp

 \brief Helper class for coupled simulations (INCA - BACI)

 <pre>
 Maintainer: Georg Hammerl
 hammerl@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15237
 </pre>
 */

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 12/11 |
 *----------------------------------------------------------------------*/

#include "aero_tfsi_delaunay.H"
#include "../headers/definitions.h"
#include "../drt_io/io_pstream.H"


void FS3I::UTILS::DelaunayTriangulation
(
  std::map<int, LINALG::Matrix<3, 1> > aerocoords,
  std::map<int,	LINALG::Matrix<3, 1> > structure,
  std::map<int,	std::vector<int> >& tris
)
{
	// Calculation of vectors forming structure element
	std::vector<double> vec1(3);
	std::vector<double> vec2(3);

	for (int i = 0; i < 3; i++)
	{
		vec1[i] = structure[1](i) - structure[0](i);
		vec2[i] = structure[2](i) - structure[0](i);
	}

	//cross-product of vec1 & vec2 to obtain structural element's normal
	std::vector<double> normal(3);
	normal[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
	normal[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
	normal[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

  //length of normal vector
  double nlength = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
#ifdef DEBUG
  //checking if structure element is a straight line
  if (nlength < EPS10)
    dserror("Structure element is 1-Dimensional - No Projection possible!");
#endif

  //unit normal vector
  normal[0] *= 1 / nlength;
  normal[1] *= 1 / nlength;
  normal[2] *= 1 / nlength;

  if(normal[1] <= 0.0)
    dserror("Either change the simple geometry so that y-dir goes along with structural normal OR do OPTIMIZATION");

  // it is important to have always a similar basis to ensure that the resulting 2D coords lead to the same ordering of the aero coords
  // hence: direction (pointing in pos. x-direction) is given by cut of structural element plane and x-y-plane --> vec
  // --> new coordinate system: [ vec1  normalxvec1  normal ]
  LINALG::Matrix<3,3> basis;
  // assumption: v in x-y-plane; given p (point) and n (normal); sought-after: v (second basis of tripod)
  // (p_0+1)*n_0 + v_1*n_1 + p_2*n_2 = p_0*n_0 + p_1*n_1 + p_2*n_2
  // arbitrarily: v_0 = p0+1 and vec_2 = p2 --> in x-y-plane
  // double v_1 = ( normal[1]*p1 - normal[0] ) / normal[1];
  // not to forget: p - v = sought-after direction (=vec1)

  vec1[0] = 1.0;
  vec1[1] = -normal[0]/normal[1];
  vec1[2] = 0.0;

  vec2[0] = normal[1] * vec1[2] - normal[2] * vec1[1];
  vec2[1] = normal[2] * vec1[0] - normal[0] * vec1[2];
  vec2[2] = normal[0] * vec1[1] - normal[1] * vec1[0];
  //length of ortho-normal basis vectors
  double vec1length = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] /* + 0.0*0.0 */);
  double vec2length = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2]);

  basis(0, 0) = vec1[0] / vec1length;
  basis(1, 0) = vec1[1] / vec1length;
  basis(2, 0) = vec1[2] / vec1length;
  basis(0, 1) = vec2[0] / vec2length;
  basis(1, 1) = vec2[1] / vec2length;
  basis(2, 1) = vec2[2] / vec2length;
  basis(0, 2) = normal[0];
  basis(1, 2) = normal[1];
  basis(2, 2) = normal[2];

  //only if matrix a is not singular
	double det = basis.Determinant();
  if (0.999 < det and det < 1.001)
    basis.Invert();
  else
    dserror("determinant of ortho-normal basis is unequal unity");

  //first point of aerocoords with normal creates plane
  const LINALG::Matrix<3,1>& firstcoord = aerocoords.begin()->second;
  double p0 = firstcoord(0);
  double p1 = firstcoord(1);
  double p2 = firstcoord(2);

  //input for delaunay triangulation
  double pointxy[2 * aerocoords.size()];
  //input for normal estimation of triangles
  double pointxyz[3 * aerocoords.size()];

  //projection of aerocoords along the normal using hesseform of plane and transformation of basis
  int i=0;
  std::vector<int> gidmapping(aerocoords.size());
  std::map<int, LINALG::Matrix<3,1> >::iterator iter;
  for(iter=aerocoords.begin(); iter!= aerocoords.end(); ++iter)
  {
    int gid = iter->first;
    LINALG::Matrix<3,1> coord = iter->second;

    //calculating distance of point to plane
    double distance = normal[0] * (coord(0) - p0) + normal[1] * (coord(1) - p1) + normal[2] * (coord(2) - p2);

    for (int j = 0; j < 3; j++)
      coord(j) -= distance * normal[j];

    // transform basis to obtain 2D problem
    LINALG::Matrix<3, 1> newcoords(true);
    newcoords.Multiply(basis, coord);

    pointxy[2 * i] = newcoords(0);
    pointxy[2 * i + 1] = newcoords(1);
    pointxyz[3 * i] = newcoords(0);
    pointxyz[3 * i + 1] = newcoords(1);
    pointxyz[3 * i + 2] = newcoords(2);

    //mapping of input data to consecutive numbering for triangulation
    gidmapping[i] = gid;

    i++;
  }

	//several parameters

	//    Input, int BASE, the base for the indexing of TRI_VERT.
	int base = 0;

	int triangle_num;
	int triangle_order = 3;
	int length_triangle = triangle_order * 3 * aerocoords.size();
	int triangle_node[length_triangle];
	int triangle_neighbor[length_triangle];

	dtris2(aerocoords.size(), base, pointxy, &triangle_num, triangle_node, triangle_neighbor);

#ifdef DEBUG
	//map for saving the normal vector of the triangles
	std::map<int, std::vector<double> > normal_triangles;
#endif
  //Check triangle normal; write output triangle_node into tris-output
	for (int i = 0; i < triangle_num; i++)
	{
		// Construction of the triangle normal
		std::vector<double> vec_test1(3);
		std::vector<double> vec_test2(3);

		int node0 = triangle_node[3 * i];
		int node1 = triangle_node[3 * i + 1];
		int node2 = triangle_node[3 * i + 2];

		for (int j = 0; j < 3; j++)
		{
			vec_test1[j] = pointxyz[3 * node1 + j] - pointxyz[3 * node0 + j];
			vec_test2[j] = pointxyz[3 * node2 + j] - pointxyz[3 * node0 + j];
		}

		std::vector<double> normal_test(3);
		normal_test[0] = (vec_test1[1] * vec_test2[2] - vec_test1[2] * vec_test2[1]);
		normal_test[1] = (vec_test1[2] * vec_test2[0] - vec_test1[0] * vec_test2[2]);
		normal_test[2] = (vec_test1[0] * vec_test2[1] - vec_test1[1] * vec_test2[0]);

		double scalar = normal[0] * normal_test[0] + normal[1] * normal_test[1]	+ normal[2] * normal_test[2];

		// normals of triangles must be oriented opposed to the structural element normal
		if (scalar < 0.0)
		{
			tris[i].push_back( gidmapping[node0] );
			tris[i].push_back( gidmapping[node1] );
			tris[i].push_back( gidmapping[node2] );
		}
		else
		{
			tris[i].push_back( gidmapping[node0] );
			tris[i].push_back( gidmapping[node2] );
			tris[i].push_back( gidmapping[node1] );

			for (int k = 0; k < 3; k++)
				normal_test[k] *= -1.0;
		}
#ifdef DEBUG
      //saving the normals of the triangles for later testing
      normal_triangles[i] = normal_test;
#endif
	}

#ifdef DEBUG
	//testing the normals of the triangle -> if parallel: n1 scalar product n2 == length_n1*length_n2
	for (int i = 0; i < triangle_num; i++)
	{
		for (int j = 0; j < triangle_num; j++)
		{
			if (i != j)
			{
				double length_i = sqrt(normal_triangles[i][0]	* normal_triangles[i][0]
				    + normal_triangles[i][1] * normal_triangles[i][1]
				    + normal_triangles[i][2] * normal_triangles[i][2]);
				double length_j = sqrt(normal_triangles[j][0]	* normal_triangles[j][0]
				    + normal_triangles[j][1] * normal_triangles[j][1]
				    + normal_triangles[j][2] * normal_triangles[j][2]);
				double absolut_product = length_i * length_j;
				double scalar_product = abs(normal_triangles[i][0]
						* normal_triangles[j][0] + normal_triangles[i][1]
						* normal_triangles[j][1] + normal_triangles[i][2]
						* normal_triangles[j][2]);
				double test = absolut_product - scalar_product; // has to be zero for parallel vectors!
				if (abs(test) > EPS10)
					dserror("Normals are NOT parallel - Fatal Error!");
			}
		}
	}
#endif

	return;
}

//****************************************************************************80

int FS3I::UTILS::diaedg(double x0, double y0, double x1, double y1, double x2,
		double y2, double x3, double y3)

//****************************************************************************80
//
//  Purpose:
//
//    DIAEDG chooses a diagonal edge.
//
//  Discussion:
//
//    The routine determines whether 0--2 or 1--3 is the diagonal edge
//    that should be chosen, based on the circumcircle criterion, where
//    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
//    quadrilateral in counterclockwise order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
//    vertices of a quadrilateral, given in counter clockwise order.
//
//    Output, int DIAEDG, chooses a diagonal:
//    +1, if diagonal edge 02 is chosen;
//    -1, if diagonal edge 13 is chosen;
//     0, if the four vertices are cocircular.
//
{
	double ca;
	double cb;
	double dx10;
	double dx12;
	double dx30;
	double dx32;
	double dy10;
	double dy12;
	double dy30;
	double dy32;
	double s;
	double tol;
	double tola;
	double tolb;
	int value;

  tol = 100.0 * EPS14;

	dx10 = x1 - x0;
	dy10 = y1 - y0;
	dx12 = x1 - x2;
	dy12 = y1 - y2;
	dx30 = x3 - x0;
	dy30 = y3 - y0;
	dx32 = x3 - x2;
	dy32 = y3 - y2;

	tola = tol * r8_max(fabs(dx10), r8_max(fabs(dy10), r8_max(fabs(dx30), fabs(
			dy30))));

	tolb = tol * r8_max(fabs(dx12), r8_max(fabs(dy12), r8_max(fabs(dx32), fabs(
			dy32))));

	ca = dx10 * dx30 + dy10 * dy30;
	cb = dx12 * dx32 + dy12 * dy32;

	if (tola < ca && tolb < cb) {
		value = -1;
	} else if (ca < -tola && cb < -tolb) {
		value = 1;
	} else {
		tola = r8_max(tola, tolb);
		s = (dx10 * dy30 - dx30 * dy10) * cb + (dx32 * dy12 - dx12 * dy32) * ca;

		if (tola < s) {
			value = -1;
		} else if (s < -tola) {
			value = 1;
		} else {
			value = 0;
		}

	}

	return value;
}
//****************************************************************************80

int FS3I::UTILS::dtris2(int point_num, int base, double point_xy[],
		int *tri_num, int tri_vert[], int tri_nabe[])

//****************************************************************************80
//
//  Purpose:
//
//    DTRIS2 constructs a Delaunay triangulation of 2D vertices.
//
//  Discussion:
//
//    The routine constructs the Delaunay triangulation of a set of 2D vertices
//    using an incremental approach and diagonal edge swaps.  Vertices are
//    first sorted in lexicographically increasing (X,Y) order, and
//    then are inserted one at a time from outside the convex hull.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int POINT_NUM, the number of vertices.
//
//    Input, int BASE, the base for the indexing of TRI_VERT.
//    0, use 0-based indexing.
//    1, use 1-based indexing.
//
//    Input/output, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//    On output, the vertices have been sorted into dictionary order.
//
//    Output, int *TRI_NUM, the number of triangles in the triangulation;
//    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
//    of boundary vertices.
//
//    Output, int TRI_VERT[TRI_NUM*3], the nodes that make up each triangle.
//    The elements are indices of POINT_XY.  The vertices of the triangles are
//    in counter clockwise order.
//
//    Output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list.
//    Positive elements are indices of TIL; negative elements are used for links
//    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
//    where I, J = triangle, edge index; TRI_NABE[I,J] refers to
//    the neighbor along edge from vertex J to J+1 (mod 3).
//
//    Output, int DTRIS2, is 0 for no error.
{
	double cmax;
	int e;
	int error;
	int i;
	int *indx;
	int j;
	int k;
	int l;
	int ledg;
	int lr;
	int ltri;
	int m;
	int m1;
	int m2;
	int n;
	int redg;
	int rtri;
	int *stack;
	int t;
	double tol;
	int top;

	stack = new int[point_num];

  tol = 100.0 * EPS14;
	//
	//  Sort the vertices by increasing (x,y).
	//
	indx = r82vec_sort_heap_index_a(point_num, base, point_xy);

	r82vec_permute(point_num, indx, base, point_xy);
	//
	//  Make sure that the data points are "reasonably" distinct.
	//
	m1 = 1;

	for (i = 2; i <= point_num; i++) {
		m = m1;
		m1 = i;

		k = -1;

		for (j = 0; j <= 1; j++) {
			cmax = r8_max(fabs(point_xy[2 * (m - 1) + j]), fabs(point_xy[2
					* (m1 - 1) + j]));

			if (tol * (cmax + 1.0) < fabs(point_xy[2 * (m - 1) + j]
					- point_xy[2 * (m1 - 1) + j])) {
				k = j;
				break;
			}

		}

		if (k == -1) {
			std::cout << "\n";
			std::cout << "DTRIS2 - Fatal error!\n";
			std::cout << "  Fails for point number I = " << i << "\n";
			std::cout << "  M =  " << m << "\n";
			std::cout << "  M1 = " << m1 << "\n";
			std::cout << "  X,Y(M)  = " << point_xy[2 * (m - 1) + 0] << "  "
					<< point_xy[2 * (m - 1) + 1] << "\n";
			std::cout << "  X,Y(M1) = " << point_xy[2 * (m1 - 1) + 0] << "  "
					<< point_xy[2 * (m1 - 1) + 1] << "\n";
			delete[] stack;
			return 224;
		}

	}
	//
	//  Starting from points M1 and M2, search for a third point M that
	//  makes a "healthy" triangle (M1,M2,M)
	//
	m1 = 1;
	m2 = 2;
	j = 3;

	for (;;) {
		if (point_num < j) {
			std::cout << "\n";
			std::cout << "DTRIS2 - Fatal error!\n";
			delete[] stack;
			return 225;
		}

		m = j;

		lr = lrline(point_xy[2 * (m - 1) + 0], point_xy[2 * (m - 1) + 1],
				point_xy[2 * (m1 - 1) + 0], point_xy[2 * (m1 - 1) + 1],
				point_xy[2 * (m2 - 1) + 0], point_xy[2 * (m2 - 1) + 1], 0.0);

		if (lr != 0) {
			break;
		}

		j = j + 1;

	}
	//
	//  Set up the triangle information for (M1,M2,M), and for any other
	//  triangles you created because points were collinear with M1, M2.
	//
	*tri_num = j - 2;

	if (lr == -1) {
		tri_vert[3 * 0 + 0] = m1;
		tri_vert[3 * 0 + 1] = m2;
		tri_vert[3 * 0 + 2] = m;
		tri_nabe[3 * 0 + 2] = -3;

		for (i = 2; i <= *tri_num; i++) {
			m1 = m2;
			m2 = i + 1;
			tri_vert[3 * (i - 1) + 0] = m1;
			tri_vert[3 * (i - 1) + 1] = m2;
			tri_vert[3 * (i - 1) + 2] = m;
			tri_nabe[3 * (i - 1) + 0] = -3 * i;
			tri_nabe[3 * (i - 1) + 1] = i;
			tri_nabe[3 * (i - 1) + 2] = i - 1;

		}

		tri_nabe[3 * (*tri_num - 1) + 0] = -3 * (*tri_num) - 1;
		tri_nabe[3 * (*tri_num - 1) + 1] = -5;
		ledg = 2;
		ltri = *tri_num;
	} else {
		tri_vert[3 * 0 + 0] = m2;
		tri_vert[3 * 0 + 1] = m1;
		tri_vert[3 * 0 + 2] = m;
		tri_nabe[3 * 0 + 0] = -4;

		for (i = 2; i <= *tri_num; i++) {
			m1 = m2;
			m2 = i + 1;
			tri_vert[3 * (i - 1) + 0] = m2;
			tri_vert[3 * (i - 1) + 1] = m1;
			tri_vert[3 * (i - 1) + 2] = m;
			tri_nabe[3 * (i - 2) + 2] = i;
			tri_nabe[3 * (i - 1) + 0] = -3 * i - 3;
			tri_nabe[3 * (i - 1) + 1] = i - 1;
		}

		tri_nabe[3 * (*tri_num - 1) + 2] = -3 * (*tri_num);
		tri_nabe[3 * 0 + 1] = -3 * (*tri_num) - 2;
		ledg = 2;
		ltri = 1;
	}
	//
	//  Insert the vertices one at a time from outside the convex hull,
	//  determine visible boundary edges, and apply diagonal edge swaps until
	//  Delaunay triangulation of vertices (so far) is obtained.
	//
	top = 0;

	for (i = j + 1; i <= point_num; i++) {
		m = i;
		m1 = tri_vert[3 * (ltri - 1) + ledg - 1];

		if (ledg <= 2) {
			m2 = tri_vert[3 * (ltri - 1) + ledg];
		} else {
			m2 = tri_vert[3 * (ltri - 1) + 0];
		}

		lr = lrline(point_xy[2 * (m - 1) + 0], point_xy[2 * (m - 1) + 1],
				point_xy[2 * (m1 - 1) + 0], point_xy[2 * (m1 - 1) + 1],
				point_xy[2 * (m2 - 1) + 0], point_xy[2 * (m2 - 1) + 1], 0.0);

		if (0 < lr) {
			rtri = ltri;
			redg = ledg;
			ltri = 0;
		} else {
			l = -tri_nabe[3 * (ltri - 1) + ledg - 1];
			rtri = l / 3;
			redg = (l % 3) + 1;
		}

		vbedg(point_xy[2 * (m - 1) + 0], point_xy[2 * (m - 1) + 1], point_num,
				point_xy, *tri_num, tri_vert, tri_nabe, &ltri, &ledg, &rtri,
				&redg);

		n = *tri_num + 1;
		l = -tri_nabe[3 * (ltri - 1) + ledg - 1];

		for (;;) {
			t = l / 3;
			e = (l % 3) + 1;
			l = -tri_nabe[3 * (t - 1) + e - 1];
			m2 = tri_vert[3 * (t - 1) + e - 1];

			if (e <= 2) {
				m1 = tri_vert[3 * (t - 1) + e];
			} else {
				m1 = tri_vert[3 * (t - 1) + 0];
			}

			*tri_num = *tri_num + 1;
			tri_nabe[3 * (t - 1) + e - 1] = *tri_num;
			tri_vert[3 * (*tri_num - 1) + 0] = m1;
			tri_vert[3 * (*tri_num - 1) + 1] = m2;
			tri_vert[3 * (*tri_num - 1) + 2] = m;
			tri_nabe[3 * (*tri_num - 1) + 0] = t;
			tri_nabe[3 * (*tri_num - 1) + 1] = *tri_num - 1;
			tri_nabe[3 * (*tri_num - 1) + 2] = *tri_num + 1;
			top = top + 1;

			if (point_num < top) {
				std::cout << "\n";
				std::cout << "DTRIS2 - Fatal error!\n";
				std::cout << "  Stack overflow.\n";
				delete[] stack;
				return 8;
			}

			stack[top - 1] = *tri_num;

			if (t == rtri && e == redg) {
				break;
			}

		}

		tri_nabe[3 * (ltri - 1) + ledg - 1] = -3 * n - 1;
		tri_nabe[3 * (n - 1) + 1] = -3 * (*tri_num) - 2;
		tri_nabe[3 * (*tri_num - 1) + 2] = -l;
		ltri = n;
		ledg = 2;

		error = swapec(m, &top, &ltri, &ledg, point_num, point_xy, *tri_num,
				tri_vert, tri_nabe, stack);

		if (error != 0) {
			std::cout << "\n";
			std::cout << "DTRIS2 - Fatal error!\n";
			std::cout << "  Error return from SWAPEC.\n";
			delete[] stack;
			return error;
		}

	}
	//
	//  Now account for the sorting that we did.
	//
	for (i = 0; i < 3; i++) {
		for (j = 0; j < *tri_num; j++) {
			tri_vert[i + j * 3] = indx[tri_vert[i + j * 3] - 1];
		}
	}

	perm_inverse(point_num, indx);

	r82vec_permute(point_num, indx, base, point_xy);

	delete[] indx;
	delete[] stack;

	return 0;
}

//****************************************************************************80

int FS3I::UTILS::i4_modp(int i, int j)

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
	int value;

	if (j == 0) {
		std::cout << "\n";
		std::cout << "I4_MODP - Fatal error!\n";
		std::cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
		exit(1);
	}

	value = i % j;

	if (value < 0) {
		value = value + abs(j);
	}

	return value;
}
//****************************************************************************80

int FS3I::UTILS::i4_sign(int i)

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
	int value;

	if (i < 0) {
		value = -1;
	} else {
		value = 1;
	}
	return value;
}
//****************************************************************************80

int FS3I::UTILS::i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int FS3I::UTILS::i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int FS3I::UTILS::i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

int FS3I::UTILS::i4vec_min(int n, int a[])

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MIN returns the value of the minimum element in an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MIN, the value of the minimum element.  This
//    is set to 0 if N <= 0.
//
{
	int i;
	int value;

	if (n <= 0) {
		return 0;
	}

	value = a[0];

	for (i = 1; i < n; i++) {
		if (a[i] < value) {
			value = a[i];
		}
	}
	return value;
}
//****************************************************************************80

int FS3I::UTILS::lrline(double xu, double yu, double xv1, double yv1,
		double xv2, double yv2, double dv)

//****************************************************************************80
//
//  Purpose:
//
//    LRLINE determines where a point lies in relation to a directed line.
//
//  Discussion:
//
//    LRLINE determines whether a point is to the left of, right of,
//    or on a directed line parallel to a line through given points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
//    directed line is parallel to and at signed distance DV to the left of
//    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
//    which the position relative to the directed line is to be determined.
//
//    Input, double DV, the signed distance, positive for left.
//
//    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
//    to the right of, on, or left of the directed line.  LRLINE is 0 if
//    the line degenerates to a point.
//
{
	double dx;
	double dxu;
	double dy;
	double dyu;
	double t;
	double tol;
	double tolabs;
	int value;

  tol = 100.0 * EPS14;

	dx = xv2 - xv1;
	dy = yv2 - yv1;
	dxu = xu - xv1;
	dyu = yu - yv1;

	tolabs = tol * r8_max(fabs(dx), r8_max(fabs(dy), r8_max(fabs(dxu), r8_max(
			fabs(dyu), fabs(dv)))));

	t = dy * dxu - dx * dyu + dv * sqrt(dx * dx + dy * dy);

	if (tolabs < t)
		value = 1;
	else if (-tolabs <= t)
		value = 0;
	else // (t < -tolabs)
		value = -1;

	return value;
}
//****************************************************************************80

bool FS3I::UTILS::perm_check(int n, int p[], int base)

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from BASE to
//    to BASE+N-1 occurs among the N entries of the permutation.
//
//    Set the input quantity BASE to 0, if P is a 0-based permutation,
//    or to 1 if P is a 1-based permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Input, int BASE, the index base.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
	bool found;
	int i;
	int seek;

	for (seek = base; seek < base + n; seek++) {
		found = false;

		for (i = 0; i < n; i++) {
			if (p[i] == seek) {
				found = true;
				break;
			}
		}

		if (!found) {
			return false;
		}

	}

	return true;
}
//****************************************************************************80

void FS3I::UTILS::perm_inverse(int n, int p[])

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INVERSE inverts a permutation "in place".
//
//  Discussion:
//
//    This algorithm assumes that the entries in the permutation vector are
//    strictly positive.  In particular, the value 0 must no occur.
//
//    When necessary, this function shifts the data temporarily so that
//    this requirement is satisfied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
	int base;
	int i;
	int i0;
	int i1;
	int i2;
	int is;
	int p_min;

	if (n <= 0) {
		std::cout << "\n";
		std::cout << "PERM_INVERSE - Fatal error!\n";
		std::cout << "  Input value of N = " << n << "\n";
		exit(1);
	}
	//
	//  Find the least value, and shift data so it begins at 1.
	//
	p_min = i4vec_min(n, p);
	base = 1;

	for (i = 0; i < n; i++) {
		p[i] = p[i] - p_min + base;
	}
	//
	//  Now we can safely check the permutation.
	//
	if (!perm_check(n, p, base)) {
		std::cerr << "\n";
		std::cerr << "PERM_INVERSE - Fatal error!\n";
		std::cerr << "  PERM_CHECK rejects this permutation.\n";
		exit(1);
	}
	//
	//  Now we can invert the permutation.
	//
	is = 1;

	for (i = 1; i <= n; i++) {
		i1 = p[i - 1];

		while (i < i1) {
			i2 = p[i1 - 1];
			p[i1 - 1] = -i2;
			i1 = i2;
		}

		is = -i4_sign(p[i - 1]);
		p[i - 1] = i4_sign(is) * abs(p[i - 1]);
	}

	for (i = 1; i <= n; i++) {
		i1 = -p[i - 1];

		if (0 <= i1) {
			i0 = i;

			for (;;) {
				i2 = p[i1 - 1];
				p[i1 - 1] = i0;

				if (i2 < 0) {
					break;
				}
				i0 = i1;
				i1 = i2;
			}
		}
	}
	//
	//  Now we can restore the permutation.
	//
	for (i = 0; i < n; i++) {
		p[i] = p[i] + p_min - base;
	}

	return;
}

//****************************************************************************80

double FS3I::UTILS::r8_max(double x, double y)

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
	double value;

	if (y < x) {
		value = x;
	} else {
		value = y;
	}
	return value;
}
//****************************************************************************80

void FS3I::UTILS::r82vec_permute(int n, int p[], int base, double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PERMUTE permutes an R82VEC in place.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type double precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.
//
//    Input, int BASE, is 0 for a 0-based permutation and 1 for a 1-based permutation.
//
//    Input/output, double A[2*N], the array to be permuted.
//
{
	double a_temp[2];
	int i;
	int iget;
	int iput;
	int istart;

	if (!perm_check(n, p, base)) {
		std::cerr << "\n";
		std::cerr << "R82VEC_PERMUTE - Fatal error!\n";
		std::cerr << "  PERM_CHECK rejects this permutation.\n";
		exit(1);
	}
	//
	//  In order for the sign negation trick to work, we need to assume that the
	//  entries of P are strictly positive.  Presumably, the lowest number is BASE.
	//  So temporarily add 1-BASE to each entry to force positivity.
	//
	for (i = 0; i < n; i++) {
		p[i] = p[i] + 1 - base;
	}
	//
	//  Search for the next element of the permutation that has not been used.
	//
	for (istart = 1; istart <= n; istart++) {
		if (p[istart - 1] < 0) {
			continue;
		} else if (p[istart - 1] == istart) {
			p[istart - 1] = -p[istart - 1];
			continue;
		} else {
			a_temp[0] = a[0 + (istart - 1) * 2];
			a_temp[1] = a[1 + (istart - 1) * 2];
			iget = istart;
			//
			//  Copy the new value into the vacated entry.
			//
			for (;;) {
				iput = iget;
				iget = p[iget - 1];

				p[iput - 1] = -p[iput - 1];

				if (iget < 1 || n < iget) {
					std::cout << "\n";
					std::cout << "R82VEC_PERMUTE - Fatal error!\n";
					std::cout << "  Entry IPUT = " << iput
							<< " of the permutation has\n";
					std::cout << "  an illegal value IGET = " << iget << ".\n";
					exit(1);
				}

				if (iget == istart) {
					a[0 + (iput - 1) * 2] = a_temp[0];
					a[1 + (iput - 1) * 2] = a_temp[1];
					break;
				}
				a[0 + (iput - 1) * 2] = a[0 + (iget - 1) * 2];
				a[1 + (iput - 1) * 2] = a[1 + (iget - 1) * 2];
			}
		}
	}
	//
	//  Restore the signs of the entries.
	//
	for (i = 0; i < n; i++) {
		p[i] = -p[i];
	}
	//
	//  Restore the base of the entries.
	//
	for (i = 0; i < n; i++) {
		p[i] = p[i] - 1 + base;
	}
	return;
}
//****************************************************************************80

int *FS3I::UTILS::r82vec_sort_heap_index_a(int n, int base, double a[])

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
//
//  Discussion:
//
//    An R82VEC is a vector whose entries are R82's.
//    An R82 is a vector of type double precision with two entries.
//    An R82VEC may be stored as a 2 by N array.
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      a(*,indx(*))
//
//    or explicitly, by the call
//
//      r82vec_permute ( n, indx, base, a )
//
//    after which a(*,*) is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int BASE, the desired indexing for the sort index:
//    0 for 0-based indexing,
//    1 for 1-based indexing.
//
//    Input, double A[2*N], an array to be index-sorted.
//
//    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
//    I-th element of the sorted array is A(0:1,R8VEC_SORT_HEAP_INDEX_A(I)).
//
{
	double aval[2];
	int i;
	int *indx;
	int indxt;
	int ir;
	int j;
	int l;

	if (n < 1) {
		return NULL;
	}

	indx = new int[n];

	for (i = 0; i < n; i++) {
		indx[i] = i;
	}

	if (n == 1) {
		indx[0] = indx[0] + base;
		return indx;
	}

	l = n / 2 + 1;
	ir = n;

	for (;;) {
		if (1 < l) {
			l = l - 1;
			indxt = indx[l - 1];
			aval[0] = a[0 + indxt * 2];
			aval[1] = a[1 + indxt * 2];
		} else {
			indxt = indx[ir - 1];
			aval[0] = a[0 + indxt * 2];
			aval[1] = a[1 + indxt * 2];
			indx[ir - 1] = indx[0];
			ir = ir - 1;

			if (ir == 1) {
				indx[0] = indxt;
				break;
			}
		}
		i = l;
		j = l + l;

		while (j <= ir) {
			if (j < ir) {
				if (a[0 + indx[j - 1] * 2] < a[0 + indx[j] * 2] || (a[0
						+ indx[j - 1] * 2] == a[0 + indx[j] * 2] && a[1
						+ indx[j - 1] * 2] < a[1 + indx[j] * 2])) {
					j = j + 1;
				}
			}

			if (aval[0] < a[0 + indx[j - 1] * 2] || (aval[0] == a[0 + indx[j
					- 1] * 2] && aval[1] < a[1 + indx[j - 1] * 2])) {
				indx[i - 1] = indx[j - 1];
				i = j;
				j = j + j;
			} else {
				j = ir + 1;
			}
		}
		indx[i - 1] = indxt;
	}
	//
	//  Take care of the base.
	//
	for (i = 0; i < n; i++) {
		indx[i] = indx[i] + base;
	}

	return indx;
}
//****************************************************************************80

int FS3I::UTILS::swapec(int i, int *top, int *btri, int *bedg, int point_num,
		double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
		int stack[])

//****************************************************************************80
//
//  Purpose:
//
//    SWAPEC swaps diagonal edges until all triangles are Delaunay.
//
//  Discussion:
//
//    The routine swaps diagonal edges in a 2D triangulation, based on
//    the empty circumcircle criterion, until all triangles are Delaunay,
//    given that I is the index of the new vertex added to the triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int I, the index of the new vertex.
//
//    Input/output, int *TOP, the index of the top of the stack.
//    On output, TOP is zero.
//
//    Input/output, int *BTRI, *BEDG; on input, if positive, are the
//    triangle and edge indices of a boundary edge whose updated indices
//    must be recorded.  On output, these may be updated because of swaps.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the points.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input/output, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//    May be updated on output because of swaps.
//
//    Input/output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list;
//    negative values are used for links of the counter-clockwise linked
//    list of boundary edges;  May be updated on output because of swaps.
//
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
//    contain the indices of initial triangles (involving vertex I)
//    put in stack; the edges opposite I should be in interior;  entries
//    TOP+1 through MAXST are used as a stack.
//
//    Output, int SWAPEC, is set to 8 for abnormal return.
//
{
	int a;
	int b;
	int c;
	int e;
	int ee;
	int em1;
	int ep1;
	int f;
	int fm1;
	int fp1;
	int l;
	int r;
	int s;
	int swap;
	int t;
	int tt;
	int u;
	double x;
	double y;
	//
	//  Determine whether triangles in stack are Delaunay, and swap
	//  diagonal edge of convex quadrilateral if not.
	//
	x = point_xy[2 * (i - 1) + 0];
	y = point_xy[2 * (i - 1) + 1];

	for (;;) {
		if (*top <= 0) {
			break;
		}

		t = stack[(*top) - 1];
		*top = *top - 1;

		if (tri_vert[3 * (t - 1) + 0] == i) {
			e = 2;
			b = tri_vert[3 * (t - 1) + 2];
		} else if (tri_vert[3 * (t - 1) + 1] == i) {
			e = 3;
			b = tri_vert[3 * (t - 1) + 0];
		} else {
			e = 1;
			b = tri_vert[3 * (t - 1) + 1];
		}

		a = tri_vert[3 * (t - 1) + e - 1];
		u = tri_nabe[3 * (t - 1) + e - 1];

		if (tri_nabe[3 * (u - 1) + 0] == t) {
			f = 1;
			c = tri_vert[3 * (u - 1) + 2];
		} else if (tri_nabe[3 * (u - 1) + 1] == t) {
			f = 2;
			c = tri_vert[3 * (u - 1) + 0];
		} else {
			f = 3;
			c = tri_vert[3 * (u - 1) + 1];
		}

		swap = diaedg(x, y, point_xy[2 * (a - 1) + 0],
				point_xy[2 * (a - 1) + 1], point_xy[2 * (c - 1) + 0],
				point_xy[2 * (c - 1) + 1], point_xy[2 * (b - 1) + 0],
				point_xy[2 * (b - 1) + 1]);

		if (swap == 1) {
			em1 = i4_wrap(e - 1, 1, 3);
			ep1 = i4_wrap(e + 1, 1, 3);
			fm1 = i4_wrap(f - 1, 1, 3);
			fp1 = i4_wrap(f + 1, 1, 3);

			tri_vert[3 * (t - 1) + ep1 - 1] = c;
			tri_vert[3 * (u - 1) + fp1 - 1] = i;
			r = tri_nabe[3 * (t - 1) + ep1 - 1];
			s = tri_nabe[3 * (u - 1) + fp1 - 1];
			tri_nabe[3 * (t - 1) + ep1 - 1] = u;
			tri_nabe[3 * (u - 1) + fp1 - 1] = t;
			tri_nabe[3 * (t - 1) + e - 1] = s;
			tri_nabe[3 * (u - 1) + f - 1] = r;

			if (0 < tri_nabe[3 * (u - 1) + fm1 - 1]) {
				*top = *top + 1;
				stack[(*top) - 1] = u;
			}

			if (0 < s) {
				if (tri_nabe[3 * (s - 1) + 0] == u) {
					tri_nabe[3 * (s - 1) + 0] = t;
				} else if (tri_nabe[3 * (s - 1) + 1] == u) {
					tri_nabe[3 * (s - 1) + 1] = t;
				} else {
					tri_nabe[3 * (s - 1) + 2] = t;
				}

				*top = *top + 1;

				if (point_num < *top) {
					return 8;
				}

				stack[(*top) - 1] = t;
			} else {
				if (u == *btri && fp1 == *bedg) {
					*btri = t;
					*bedg = e;
				}

				l = -(3 * t + e - 1);
				tt = t;
				ee = em1;

				while (0 < tri_nabe[3 * (tt - 1) + ee - 1]) {
					tt = tri_nabe[3 * (tt - 1) + ee - 1];

					if (tri_vert[3 * (tt - 1) + 0] == a) {
						ee = 3;
					} else if (tri_vert[3 * (tt - 1) + 1] == a) {
						ee = 1;
					} else {
						ee = 2;
					}
				}
				tri_nabe[3 * (tt - 1) + ee - 1] = l;
			}

			if (0 < r) {
				if (tri_nabe[3 * (r - 1) + 0] == t) {
					tri_nabe[3 * (r - 1) + 0] = u;
				} else if (tri_nabe[3 * (r - 1) + 1] == t) {
					tri_nabe[3 * (r - 1) + 1] = u;
				} else {
					tri_nabe[3 * (r - 1) + 2] = u;
				}
			} else {
				if (t == *btri && ep1 == *bedg) {
					*btri = u;
					*bedg = f;
				}

				l = -(3 * u + f - 1);
				tt = u;
				ee = fm1;

				while (0 < tri_nabe[3 * (tt - 1) + ee - 1]) {
					tt = tri_nabe[3 * (tt - 1) + ee - 1];

					if (tri_vert[3 * (tt - 1) + 0] == b) {
						ee = 3;
					} else if (tri_vert[3 * (tt - 1) + 1] == b) {
						ee = 1;
					} else {
						ee = 2;
					}

				}
				tri_nabe[3 * (tt - 1) + ee - 1] = l;
			}
		}
	}

	return 0;
}
//****************************************************************************80

void FS3I::UTILS::vbedg(double x, double y, int point_num, double point_xy[],
		int tri_num, int tri_vert[], int tri_nabe[], int *ltri, int *ledg,
		int *rtri, int *redg)

//****************************************************************************80
//
//  Purpose:
//
//    VBEDG determines which boundary edges are visible to a point.
//
//  Discussion:
//
//    The point (X,Y) is assumed to be outside the convex hull of the
//    region covered by the 2D triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2008
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Modified:
//
//    02 September 2003
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point outside the convex hull
//    of the current triangulation.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//
//    Input, int TRI_NABE[TRI_NUM*3], the triangle neighbor list; negative
//    values are used for links of a counter clockwise linked list of boundary
//    edges;
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
//    assumed to be already computed and are not changed, else they are updated.
//    On output, LTRI is the index of boundary triangle to the left of the
//    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
//    edge of triangle LTRI to the left of the leftmost boundary edge visible
//    from (X,Y).  1 <= LEDG <= 3.
//
//    Input/output, int *RTRI.  On input, the index of the boundary triangle
//    to begin the search at.  On output, the index of the rightmost boundary
//    triangle visible from (X,Y).
//
//    Input/output, int *REDG, the edge of triangle RTRI that is visible
//    from (X,Y).  1 <= REDG <= 3.
//
{
	int a;
	double ax;
	double ay;
	int b;
	double bx;
	double by;
	bool done;
	int e;
	int l;
	int lr;
	int t;
	//
	//  Find the rightmost visible boundary edge using links, then possibly
	//  leftmost visible boundary edge using triangle neighbor information.
	//
	if (*ltri == 0) {
		done = false;
		*ltri = *rtri;
		*ledg = *redg;
	} else {
		done = true;
	}

	for (;;) {
		l = -tri_nabe[3 * ((*rtri) - 1) + (*redg) - 1];
		t = l / 3;
		e = 1 + l % 3;
		a = tri_vert[3 * (t - 1) + e - 1];

		if (e <= 2) {
			b = tri_vert[3 * (t - 1) + e];
		} else {
			b = tri_vert[3 * (t - 1) + 0];
		}

		ax = point_xy[2 * (a - 1) + 0];
		ay = point_xy[2 * (a - 1) + 1];

		bx = point_xy[2 * (b - 1) + 0];
		by = point_xy[2 * (b - 1) + 1];

		lr = lrline(x, y, ax, ay, bx, by, 0.0);

		if (lr <= 0) {
			break;
		}

		*rtri = t;
		*redg = e;

	}

	if (done) {
		return;
	}

	t = *ltri;
	e = *ledg;

	for (;;) {
		b = tri_vert[3 * (t - 1) + e - 1];
		e = i4_wrap(e - 1, 1, 3);

		while (0 < tri_nabe[3 * (t - 1) + e - 1]) {
			t = tri_nabe[3 * (t - 1) + e - 1];

			if (tri_vert[3 * (t - 1) + 0] == b) {
				e = 3;
			} else if (tri_vert[3 * (t - 1) + 1] == b) {
				e = 1;
			} else {
				e = 2;
			}

		}

		a = tri_vert[3 * (t - 1) + e - 1];
		ax = point_xy[2 * (a - 1) + 0];
		ay = point_xy[2 * (a - 1) + 1];

		bx = point_xy[2 * (b - 1) + 0];
		by = point_xy[2 * (b - 1) + 1];

		lr = lrline(x, y, ax, ay, bx, by, 0.0);

		if (lr <= 0) {
			break;
		}

	}

	*ltri = t;
	*ledg = e;

	return;
}

