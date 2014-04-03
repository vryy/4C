/*!----------------------------------------------------------------------
\file volmortar_cell.cpp

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "volmortar_cell.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

/*---------------------------------------------------------------------*
 | constructor                                             farah 01/14 |
 *---------------------------------------------------------------------*/
VOLMORTAR::Cell::Cell(int id,
                      int nvertices,
                      const Epetra_SerialDenseMatrix& coords,
                      const DRT::Element::DiscretizationType& shape) :
id_(id),
coords_(coords),
shape_(shape)
{
  /* get the matrix of the coordinates of edges needed to compute the volume,
  ** which is used here as detJ in the quadrature rule.
  ** ("Jacobian matrix") for the quadrarture rule:
  **             [  1    1    1    1  ]
  ** jac_coord = [ x_1  x_2  x_3  x_4 ]
  **             [ y_1  y_2  y_3  y_4 ]
  **             [ z_1  z_2  z_3  z_4 ]
  */

  LINALG::Matrix<4,4> jac;
  for (int i=0; i<4; i++)  jac(0,i)=1;
  for (int row=0;row<3;row++)
    for (int col=0;col<4;col++)
      jac(row+1,col)= coords_(row,col);

  // volume of the element
  vol_ = jac.Determinant()/6.0;
  if (vol_ <= 0.0)
    dserror("Element volume %10.5e <= 0.0",vol_);

  return;
}


/*---------------------------------------------------------------------*
 | mapping from parameter space to global space            farah 01/14 |
 *---------------------------------------------------------------------*/
void VOLMORTAR::Cell::LocalToGlobal(double* local, double* global)
{
  // check input
  if (!local) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!global) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  static const int n    = 4;
  static const int ndim = 3;

  for (int i=0;i<ndim;++i)
    global[i]=0.0;

  LINALG::Matrix<n,1>      val;
  DRT::UTILS::shape_function_3D(val,local[0],local[1],local[2],shape_);

  for (int i=0;i<n;++i)
  {
    for(int j=0;j<ndim;++j)
    {
      // use shape function values for interpolation
      global[j] += val(i)*coords_(j,i);
    }
  }


  return;
}

/*---------------------------------------------------------------------*
 | output                                                  farah 03/14 |
 *---------------------------------------------------------------------*/
void VOLMORTAR::Cell::Print()
{
  for(int i=0;i<4;++i)
    std::cout << "coords= " << coords_(0,i) << " "<< coords_(1,i) << " "<< coords_(2,i) << std::endl;

  return;
}
