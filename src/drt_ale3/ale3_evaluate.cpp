//-----------------------------------------------------------------------
/*!
\file ale3_evaluate.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include <mpi.h>
#endif

#include "ale3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/stvenantkirchhoff.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale3::Evaluate(ParameterList& params,
                                  DRT::Discretization&      discretization,
                                  vector<int>&              lm,
                                  Epetra_SerialDenseMatrix& elemat1,
                                  Epetra_SerialDenseMatrix& elemat2,
                                  Epetra_SerialDenseVector& elevec1,
                                  Epetra_SerialDenseVector& elevec2,
                                  Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Ale3::ActionType act = Ale3::none;

  // get the action required
  string action = params.get<string>("action","none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_ale_lin_stiff")
    act = Ale3::calc_ale_lin_stiff;
  else if (action == "calc_ale_spring")
    act = Ale3::calc_ale_spring;
  else
    dserror("Unknown type of action for Ale3");

  // get the material
  RefCountPtr<MAT::Material> mat = Material();

  switch(act)
  {
  case calc_ale_lin_stiff:
  {
    //RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    //vector<double> my_dispnp(lm.size());
    //DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

    static_ke(lm,&elemat1,&elevec1,mat,params);

    break;
  }

  case calc_ale_spring:
  {
    RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    vector<double> my_dispnp(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

    static_ke_spring(&elemat1,my_dispnp);

    break;
  }


  default:
    dserror("Unknown type of action for Ale3");
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the ale elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

void DRT::ELEMENTS::Ale3::ale3_edge_geometry(int i, int j, const Epetra_SerialDenseMatrix& xyze,
                                             double* length,
                                             double* edge_x,
                                             double* edge_y,
                                             double* edge_z)
{
  double delta_x, delta_y, delta_z;
  /*---------------------------------------------- x-, y- and z-difference ---*/
  delta_x = xyze(0,j)-xyze(0,i);
  delta_y = xyze(1,j)-xyze(1,i);
  delta_z = xyze(2,j)-xyze(2,i);
  /*------------------------------- determine distance between i and j ---*/
  *length = sqrt( delta_x * delta_x
                + delta_y * delta_y
                + delta_z * delta_z);
  if (*length < (1.0E-14)) dserror("edge or diagonal of element has zero length");
  /*--------------------------------------- determine direction of edge i-j ---*/
  *edge_x = delta_x / *length;
  *edge_y = delta_y / *length;
  *edge_z = delta_z / *length;
  /*----------------------------------------------------------------------*/
}

double DRT::ELEMENTS::Ale3::ale3_area_tria(const Epetra_SerialDenseMatrix& xyze,
                                           int i, int j, int k)
{
  double a, b, c;  /* geometrical values */
  double area;  /* triangle area */
  /*----------------------------------------------------------------------*/

  a = (xyze(0,i)-xyze(0,j))*(xyze(0,i)-xyze(0,j))
     +(xyze(1,i)-xyze(1,j))*(xyze(1,i)-xyze(1,j)); /* line i-j squared */
  b = (xyze(0,j)-xyze(0,k))*(xyze(0,j)-xyze(0,k))
     +(xyze(1,j)-xyze(1,k))*(xyze(1,j)-xyze(1,k)); /* line j-k squared */
  c = (xyze(0,k)-xyze(0,i))*(xyze(0,k)-xyze(0,i))
     +(xyze(1,k)-xyze(1,i))*(xyze(1,k)-xyze(1,i)); /* line k-i squared */
  area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
  return area;
}

void DRT::ELEMENTS::Ale3::ale3_torsional(int i, int j, int k,
                                         const Epetra_SerialDenseMatrix& xyze,
                                         Epetra_SerialDenseMatrix* k_torsion)
{
/*
                           k
                           *
	                  / \
       y,v ^	    l_ki /   \  l_jk
           |	       	/     \
	    --->     i *-------* j
            x,u	        l_ij
*/

  double x_ij, x_jk, x_ki;  /* x-differences between nodes */
  double y_ij, y_jk, y_ki;  /* y-differences between nodes */
  double l_ij, l_jk, l_ki;  /* side lengths */
  double a_ij, a_jk, a_ki;  /* auxiliary values same as in Farhat et al. */
  double b_ij, b_jk, b_ki;  /*                  - " -                    */
  double area;              /* area of the triangle */


  Epetra_SerialDenseMatrix R(3,6);   /* rotation matrix same as in Farhat et al. */
  Epetra_SerialDenseMatrix C(3,3);   /* torsion stiffness matrix same as in Farhat et al. */
  Epetra_SerialDenseMatrix A(6,3);   /* auxiliary array of intermediate results */


/*--------------------------------- determine basic geometric values ---*/
  x_ij = xyze(0,j) - xyze(0,i);
  x_jk = xyze(0,k) - xyze(0,j);
  x_ki = xyze(0,i) - xyze(0,k);
  y_ij = xyze(1,j) - xyze(1,i);
  y_jk = xyze(1,k) - xyze(1,j);
  y_ki = xyze(1,i) - xyze(1,k);

  l_ij = sqrt( x_ij*x_ij + y_ij*y_ij );
  l_jk = sqrt( x_jk*x_jk + y_jk*y_jk );
  l_ki = sqrt( x_ki*x_ki + y_ki*y_ki );

/*----------------------------------------------- check edge lengths ---*/
  if (l_ij < (1.0E-14)) dserror("edge or diagonal of element has zero length");
  if (l_jk < (1.0E-14)) dserror("edge or diagonal of element has zero length");
  if (l_ki < (1.0E-14)) dserror("edge or diagonal of element has zero length");

/*-------------------------------------------- fill auxiliary values ---*/
  a_ij = x_ij / (l_ij*l_ij);
  a_jk = x_jk / (l_jk*l_jk);
  a_ki = x_ki / (l_ki*l_ki);
  b_ij = y_ij / (l_ij*l_ij);
  b_jk = y_jk / (l_jk*l_jk);
  b_ki = y_ki / (l_ki*l_ki);

/*--------------------------------------------------- determine area ---*/
  area = ale3_area_tria(xyze,i,j,k);

/*---------------------------------- determine torsional stiffnesses ---*/
  C(0,0) = l_ij*l_ij * l_ki*l_ki / (4.0*area*area);
  C(1,1) = l_ij*l_ij * l_jk*l_jk / (4.0*area*area);
  C(2,2) = l_ki*l_ki * l_jk*l_jk / (4.0*area*area);

/*--------------------------------------- fill transformation matrix ---*/
  R(0,0) = - b_ki - b_ij;
  R(0,1) = a_ij + a_ki;
  R(0,2) = b_ij;
  R(0,3) = - a_ij;
  R(0,4) = b_ki;
  R(0,5) = - a_ki;

  R(1,0) = b_ij;
  R(1,1) = - a_ij;
  R(1,2) = - b_ij - b_jk;
  R(1,3) = a_jk + a_ij;
  R(1,4) = b_jk;
  R(1,5) = - a_jk;

  R(2,0) = b_ki;
  R(2,1) = - a_ki;
  R(2,2) = b_jk;
  R(2,3) = - a_jk;
  R(2,4) = - b_jk - b_ki;
  R(2,5) = a_ki + a_jk;

/*----------------------------------- perform matrix multiplications ---*/


  int err = A.Multiply('T','N',1,R,C,0);	// A = R^t * C
  if (err!=0)
    dserror("Multiply failed");
  err = k_torsion->Multiply('N','N',1,A,R,0);	// stiff = A * R
  if (err!=0)
    dserror("Multiply failed");

}

void DRT::ELEMENTS::Ale3::ale3_add_tria_stiffness(int node_p, int node_q, int node_r, int node_s,
                                                  Epetra_SerialDenseMatrix* sys_mat,
                                                  const Epetra_SerialDenseMatrix& xyze)
{
  const DiscretizationType distype = this->Shape();       // Discretization type of this element
  int numcnd;                                             // number of corner nodes

  //Positions and local rotational stiffness matrix for dynamic triangle (2D)
  //sequence: s,j,q
  Epetra_SerialDenseMatrix xyze_dyn_tria(2,3);
  Epetra_SerialDenseMatrix k_dyn_tria(6,6);

  //auxiliary matrices
  Epetra_SerialDenseMatrix A(9,6);
  Epetra_SerialDenseMatrix B(12,9);
  Epetra_SerialDenseMatrix S(12,9);

  //transformation matrix from the plane of the triangle to the
  //three-dimensional global frame
  Epetra_SerialDenseMatrix trans_matrix(6,9);

  //rotational stiffness matrix for dynamic triangle in global frame (3D)
  Epetra_SerialDenseMatrix k_dyn_tria_global(9,9);

  //rotational stiffness matrix for tetrahedron with given dynamic triangle
  Epetra_SerialDenseMatrix k_dyn_tet(12,12);

  //local x,y in the plane of the dynamic triangle
  blitz::Array<double,1> local_x(3);
  blitz::Array<double,1> local_y(3);
  local_x = 0.0;
  local_y = 0.0;

  //auxiliary vectors
  blitz::Array<double,1> sq(3);
  blitz::Array<double,1> pq(3);
  blitz::Array<double,1> rp(3);
  blitz::Array<double,1> pj(3);
  blitz::Array<double,1> res(3);
  sq = 0.0;
  pq = 0.0;
  rp = 0.0;
  pj = 0.0;
  res = 0.0;

  double length, factor, ex, ey, ez;

  //number of corner nodes
  switch (distype)
  {
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    numcnd = 4;
    break;
  case DRT::Element::pyramid5:
    numcnd = 5;
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    numcnd = 6;
    break;
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    numcnd = 8;
    break;
  default:
    numcnd = 0;
    dserror("distype unkown");
    break;
  }

  ale3_edge_geometry(node_s,node_q,xyze,&length,&sq(0),&sq(1),&sq(2));
  sq = length*sq;
  ale3_edge_geometry(node_p,node_q,xyze,&length,&pq(0),&pq(1),&pq(2));
  pq = length*pq;
  ale3_edge_geometry(node_r,node_p,xyze,&length,&rp(0),&rp(1),&rp(2));
  rp = length*rp;


  //local_x := normal vector of face pqr = pq x rp
  local_x(0) = pq(1)*rp(2)-pq(2)*rp(1);  //just an intermediate step
  local_x(1) = pq(2)*rp(0)-pq(0)*rp(2);  //just an intermediate step
  local_x(2) = pq(0)*rp(1)-pq(1)*rp(0);  //just an intermediate step
  length = sqrt(local_x(0)*local_x(0)+local_x(1)*local_x(1)+local_x(2)*local_x(2));
  local_x = (1.0/length)*local_x;

  //local x-value of s xyze_dyn_tria(0,0) := (-1.0)*sq*local_x (local origin lies on plane pqr)
  //local y-value of s xyze_dyn_tria(1,0 := 0.0
  xyze_dyn_tria(0,0)=(-1.0)*(sq(0)*local_x(0)+sq(1)*local_x(1)+sq(2)*local_x(2));
  xyze_dyn_tria(1,0)=0.0;

  //local_y = (sq + xyze_dyn_tria(0,0)*local_x)/|(sq + xyze_dyn_tria(0,0)*local_x)|
  //xyze_dyn_tria(1,2) = |sq + xyze_dyn_tria(0,0)*local_x|, xyze_dyn_tria(0,2) := 0.0
  local_y = sq + xyze_dyn_tria(0,0)*local_x;  //just an intermediate step
  xyze_dyn_tria(1,2) = sqrt(local_y(0)*local_y(0)+local_y(1)*local_y(1)+local_y(2)*local_y(2));
  xyze_dyn_tria(0,2) = 0.0;

  //if s lies directly above q calculate stiffness of lineal spring s-q
  if (fabs(xyze_dyn_tria(1,2)/xyze_dyn_tria(0,0)) < 1e-4)
  {
    ale3_edge_geometry(node_s,node_q,xyze,&length,&ex,&ey,&ez);
    factor = 1.0 / length;
    //put values in 'element stiffness'
    //rows for node_s
    (*sys_mat)(node_s*3,  node_s*3  ) += ex*ex * factor;
    (*sys_mat)(node_s*3,  node_s*3+1) += ex*ey * factor;
    (*sys_mat)(node_s*3,  node_s*3+2) += ex*ez * factor;

    (*sys_mat)(node_s*3,  node_q*3  ) -= ex*ex * factor;
    (*sys_mat)(node_s*3,  node_q*3+1) -= ex*ey * factor;
    (*sys_mat)(node_s*3,  node_q*3+2) -= ex*ez * factor;

    (*sys_mat)(node_s*3+1,  node_s*3  ) += ex*ey * factor;
    (*sys_mat)(node_s*3+1,  node_s*3+1) += ey*ey * factor;
    (*sys_mat)(node_s*3+1,  node_s*3+2) += ey*ez * factor;

    (*sys_mat)(node_s*3+1,  node_q*3  ) -= ex*ey * factor;
    (*sys_mat)(node_s*3+1,  node_q*3+1) -= ey*ey * factor;
    (*sys_mat)(node_s*3+1,  node_q*3+2) -= ey*ez * factor;

    (*sys_mat)(node_s*3+2,  node_s*3  ) += ex*ez * factor;
    (*sys_mat)(node_s*3+2,  node_s*3+1) += ey*ez * factor;
    (*sys_mat)(node_s*3+2,  node_s*3+2) += ez*ez * factor;

    (*sys_mat)(node_s*3+2,  node_q*3  ) -= ex*ez * factor;
    (*sys_mat)(node_s*3+2,  node_q*3+1) -= ey*ez * factor;
    (*sys_mat)(node_s*3+2,  node_q*3+2) -= ez*ez * factor;

    //rows for node_q
    (*sys_mat)(node_q*3,  node_s*3  ) -= ex*ex * factor;
    (*sys_mat)(node_q*3,  node_s*3+1) -= ex*ey * factor;
    (*sys_mat)(node_q*3,  node_s*3+2) -= ex*ez * factor;

    (*sys_mat)(node_q*3,  node_q*3  ) += ex*ex * factor;
    (*sys_mat)(node_q*3,  node_q*3+1) += ex*ey * factor;
    (*sys_mat)(node_q*3,  node_q*3+2) += ex*ez * factor;

    (*sys_mat)(node_q*3+1,  node_s*3  ) -= ex*ey * factor;
    (*sys_mat)(node_q*3+1,  node_s*3+1) -= ey*ey * factor;
    (*sys_mat)(node_q*3+1,  node_s*3+2) -= ey*ez * factor;

    (*sys_mat)(node_q*3+1,  node_q*3  ) += ex*ey * factor;
    (*sys_mat)(node_q*3+1,  node_q*3+1) += ey*ey * factor;
    (*sys_mat)(node_q*3+1,  node_q*3+2) += ey*ez * factor;

    (*sys_mat)(node_q*3+2,  node_s*3  ) -= ex*ez * factor;
    (*sys_mat)(node_q*3+2,  node_s*3+1) -= ey*ez * factor;
    (*sys_mat)(node_q*3+2,  node_s*3+2) -= ez*ez * factor;

    (*sys_mat)(node_q*3+2,  node_q*3  ) += ex*ez * factor;
    (*sys_mat)(node_q*3+2,  node_q*3+1) += ey*ez * factor;
    (*sys_mat)(node_q*3+2,  node_q*3+2) += ez*ez * factor;
  }
  else
  {
    local_y = (1.0/xyze_dyn_tria(1,2))*local_y;

    //local x,y-values of j, using pO + Oj + jp = 0
    //(O is local origin on plane pqr)
    xyze_dyn_tria(0,1) = 0.0;
    res = xyze_dyn_tria(1,2)*local_y - pq;

    double numerator = 0.0;
    double denominator = 1.0;
    double f, check = 0.0;
    int solved = 0;

    //solve linear system of equations, case differentiation
    if (!solved and fabs(rp(0))>1e-14 and fabs(local_y(1))>1e-14)
    {
      numerator = (res(1)-(rp(1)/rp(0))*res(0));
      denominator = (local_y(1)-(rp(1)/rp(0))*local_y(0));

      //check if result (value j_y) is valid, using third equation
      if (fabs(denominator)>1e-6)
      {
        f = (res(0)-(numerator/denominator)*local_y(0))/rp(0);
        check = (local_y(2)*(numerator/denominator)+rp(2)*f)-res(2);
        if (fabs(check)<1e-6 and fabs(numerator/denominator)<1e10)
          solved = 1;
      }
    }
    if (!solved and fabs(rp(0))>1e-14 and fabs(local_y(2))>1e-14)
    {
      numerator = (res(2)-(rp(2)/rp(0))*res(0));
      denominator = (local_y(2)-(rp(2)/rp(0))*local_y(0));

      //check if result (value j_y) is valid, using second equation
      if (fabs(denominator)>1e-6)
      {
        f = (res(0)-(numerator/denominator)*local_y(0))/rp(0);
        check = (local_y(1)*(numerator/denominator)+rp(1)*f)-res(1);
        if (fabs(check)<1e-6 and fabs(numerator/denominator)<1e10)
          solved = 1;
      }
    }
    if (!solved and fabs(rp(1))>1e-14 and fabs(local_y(0))>1e-14)
    {
      numerator = (res(0)-(rp(0)/rp(1))*res(1));
      denominator = (local_y(0)-(rp(0)/rp(1))*local_y(1));

      //check if result (value j_y) is valid, using third equation
      if (fabs(denominator)>1e-6)
      {
        f = (res(1)-(numerator/denominator)*local_y(1))/rp(1);
        check = (local_y(2)*(numerator/denominator)+rp(2)*f)-res(2);
        if (fabs(check)<1e-6 and fabs(numerator/denominator)<1e10)
          solved = 1;
      }
    }
    if (!solved and fabs(rp(1))>1e-14 and fabs(local_y(2))>1e-14)
    {
      numerator = (res(2)-(rp(2)/rp(1))*res(1));
      denominator = (local_y(2)-(rp(2)/rp(1))*local_y(1));

      //check if result (value j_y) is valid, using first equation
      if (fabs(denominator)>1e-6)
      {
        f = (res(1)-(numerator/denominator)*local_y(1))/rp(1);
        check = (local_y(0)*(numerator/denominator)+rp(0)*f)-res(0);
        if (fabs(check)<1e-6 and fabs(numerator/denominator)<1e10)
          solved = 1;
      }
    }
    if (!solved and fabs(rp(2))>1e-14 and fabs(local_y(0))>1e-14)
    {
      numerator = (res(0)-(rp(0)/rp(2))*res(2));
      denominator = (local_y(0)-(rp(0)/rp(2))*local_y(2));

      //check if result (value j_y) is valid, using second equation
      if (fabs(denominator)>1e-6)
      {
        f = (res(2)-(numerator/denominator)*local_y(2))/rp(2);
        check = (local_y(1)*(numerator/denominator)+rp(1)*f)-res(1);
        if (fabs(check)<1e-6 and fabs(numerator/denominator)<1e10)
          solved = 1;
      }
    }
    if (!solved and fabs(rp(2))>1e-14 and fabs(local_y(1))>1e-14)
    {
      numerator = (res(1)-(rp(1)/rp(2))*res(2));
      denominator = (local_y(1)-(rp(1)/rp(2))*local_y(2));

      //check if result (value j_y) is valid, using first equation
      if (fabs(denominator)>1e-6)
      {
        f = (res(2)-(numerator/denominator)*local_y(2))/rp(2);
        check = (local_y(0)*(numerator/denominator)+rp(0)*f)-res(0);
        if (fabs(check)<1e-6 and fabs(numerator/denominator)<1e10)
          solved = 1;
      }
    }

    xyze_dyn_tria(1,1) = (numerator/denominator);

    if (solved)
    {
      //evaluate torsional stiffness of dynamic triangle
      ale3_torsional(0,1,2,xyze_dyn_tria,&k_dyn_tria);


      //transformation matrix from the three-dimensional global frame to the
      //plane of the dynamic triangle
      int row,col;
      for (row=0; row<3; row++)
      {
        col = row;
        trans_matrix(row*2,col*3)  = local_x(0);
        trans_matrix(row*2,col*3+1)= local_x(1);
        trans_matrix(row*2,col*3+2)= local_x(2);

        trans_matrix(row*2+1,col*3)  = local_y(0);
        trans_matrix(row*2+1,col*3+1)= local_y(1);
        trans_matrix(row*2+1,col*3+2)= local_y(2);
      }


      int err = A.Multiply('T','N',1,trans_matrix,k_dyn_tria,0);    	// A = trans_matrix^t * k_dyn_tria
      if (err!=0)
        dserror("Multiply failed");

      err = k_dyn_tria_global.Multiply('N','N',1,A,trans_matrix,0);     // k_dyn_tria_global = A * trans_matrix
      if (err!=0)
        dserror("Multiply failed");


      //S transfers elastic forces at s,q,j to cornernodes p,q,r,s of the
      //tetrahedron
      double lambda;
      pj = pq + (xyze_dyn_tria(1,1)-xyze_dyn_tria(1,2))*local_y;
      if ((pj(0)*rp(0)+pj(1)*rp(1)+pj(2)*rp(2))<0.0)
      {
        lambda = sqrt((pj(0)*pj(0)+pj(1)*pj(1)+pj(2)*pj(2))/
                      (rp(0)*rp(0)+rp(1)*rp(1)+rp(2)*rp(2)));
        lambda = min(1.0,lambda);
      }
      else
        lambda = 0.0;

      S(0,3)=S(1,4)=S(2,5)  = (1.0-lambda);
      S(3,6)=S(4,7)=S(5,8)  =  1.0;
      S(6,3)=S(7,4)=S(8,5)  =  lambda;
      S(9,0)=S(10,1)=S(11,2)=  1.0;


      err = B.Multiply('N','N',1,S,k_dyn_tria_global,0);	// B = S * k_dyn_tria_global
      if (err!=0)
        dserror("Multiply failed");

      err = k_dyn_tet.Multiply('N','T',1,B,S,0);    	        // k_dyn_tet = B * S^t
      if (err!=0)
        dserror("Multiply failed");


      //Sort values in element's sys_mat
      for (int elem_node_i=0; elem_node_i<numcnd; elem_node_i++)
        for (int elem_node_j=0; elem_node_j<numcnd; elem_node_j++)
          if (((elem_node_i == node_p)||(elem_node_i == node_q)||(elem_node_i == node_r)||(elem_node_i == node_s))and((elem_node_j == node_p)||(elem_node_j == node_q)||(elem_node_j == node_r)||(elem_node_j == node_s)))
          {
            int elem_node_i_tetID = ((0*((int)(elem_node_i == node_p)))+(1*((int)(elem_node_i == node_q)))+(2*((int)(elem_node_i == node_r)))+(3*((int)(elem_node_i == node_s))));
            int elem_node_j_tetID = ((0*((int)(elem_node_j == node_p)))+(1*((int)(elem_node_j == node_q)))+(2*((int)(elem_node_j == node_r)))+(3*((int)(elem_node_j == node_s))));

            (*sys_mat)(elem_node_i*3  ,elem_node_j*3  ) += k_dyn_tet(elem_node_i_tetID*3  ,elem_node_j_tetID*3  );
            (*sys_mat)(elem_node_i*3  ,elem_node_j*3+1) += k_dyn_tet(elem_node_i_tetID*3  ,elem_node_j_tetID*3+1);
            (*sys_mat)(elem_node_i*3  ,elem_node_j*3+2) += k_dyn_tet(elem_node_i_tetID*3  ,elem_node_j_tetID*3+2);

            (*sys_mat)(elem_node_i*3+1,elem_node_j*3  ) += k_dyn_tet(elem_node_i_tetID*3+1,elem_node_j_tetID*3  );
            (*sys_mat)(elem_node_i*3+1,elem_node_j*3+1) += k_dyn_tet(elem_node_i_tetID*3+1,elem_node_j_tetID*3+1);
            (*sys_mat)(elem_node_i*3+1,elem_node_j*3+2) += k_dyn_tet(elem_node_i_tetID*3+1,elem_node_j_tetID*3+2);

            (*sys_mat)(elem_node_i*3+2,elem_node_j*3  ) += k_dyn_tet(elem_node_i_tetID*3+2,elem_node_j_tetID*3  );
            (*sys_mat)(elem_node_i*3+2,elem_node_j*3+1) += k_dyn_tet(elem_node_i_tetID*3+2,elem_node_j_tetID*3+1);
            (*sys_mat)(elem_node_i*3+2,elem_node_j*3+2) += k_dyn_tet(elem_node_i_tetID*3+2,elem_node_j_tetID*3+2);
          }
    }
  }
}

void DRT::ELEMENTS::Ale3::ale3_add_tetra_stiffness(int tet_0, int tet_1, int tet_2, int tet_3,
                                                   Epetra_SerialDenseMatrix* sys_mat,
                                                   const Epetra_SerialDenseMatrix& xyze)
{
  //according to Farhat et al.
  //twelve-triangle configuration

  int node_s;
  int node_p;
  int node_q;
  int node_r;

  blitz::Array<int,1> nodeID(4);
  nodeID(0) = tet_0;
  nodeID(1) = tet_1;
  nodeID(2) = tet_2;
  nodeID(3) = tet_3;

  for (node_p=0; node_p<3; node_p++)
    for (node_r = node_p+1; node_r<4; node_r++)
      for (node_q=0; node_q<4; node_q++)
        if ((node_q != node_p)and(node_q != node_r))
          for (node_s=0; node_s<4; node_s++)
            if ((node_s != node_p)and(node_s != node_q)and(node_s != node_r))
              ale3_add_tria_stiffness(nodeID(node_p),nodeID(node_q),nodeID(node_r),nodeID(node_s),sys_mat,xyze);
}

void DRT::ELEMENTS::Ale3::ale3_tors_spring_tet4(Epetra_SerialDenseMatrix* sys_mat,
                                               const Epetra_SerialDenseMatrix& xyze)
{
  ale3_add_tetra_stiffness(0,1,2,3,sys_mat,xyze);
}

void DRT::ELEMENTS::Ale3::ale3_tors_spring_pyramid5(Epetra_SerialDenseMatrix* sys_mat,
                                                    const Epetra_SerialDenseMatrix& xyze)
{
  ale3_add_tetra_stiffness(0,1,3,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(0,1,2,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(1,2,3,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(0,2,3,4,sys_mat,xyze);
}

void DRT::ELEMENTS::Ale3::ale3_tors_spring_wedge6(Epetra_SerialDenseMatrix* sys_mat,
                                                  const Epetra_SerialDenseMatrix& xyze)
{
  int tet_0, tet_1, tet_2, tet_3;
  for (tet_0=0; tet_0<3; tet_0++)
    for (tet_1=tet_0+1; tet_1<4; tet_1++)
      for (tet_2=tet_1+1; tet_2<5; tet_2++)
        for (tet_3=tet_2+1; tet_3<6; tet_3++)
        {
          if (!((tet_0==0 and tet_1==1 and tet_2==3 and tet_3==4)||
                (tet_0==0 and tet_1==2 and tet_2==3 and tet_3==5)||
                (tet_0==1 and tet_1==2 and tet_2==4 and tet_3==5)))
            ale3_add_tetra_stiffness(tet_2, tet_0, tet_1, tet_3, sys_mat, xyze);
        }
}

void DRT::ELEMENTS::Ale3::ale3_tors_spring_hex8(Epetra_SerialDenseMatrix* sys_mat,
                                               const Epetra_SerialDenseMatrix& xyze)
{
//Use all valid tetrahedra to prevent node-face-penetration. This is not working well
//   int tet_0, tet_1, tet_2, tet_3;

//   for (tet_0=0; tet_0<5; tet_0++)
//     for (tet_1=tet_0+1; tet_1<6; tet_1++)
//       for (tet_2=tet_1+1; tet_2<7; tet_2++)
//         for (tet_3=tet_2+1; tet_3<8; tet_3++)
//         {
//           //if the 4 nodes don't belong to a face or diagonal plane of the hexahedron they form
//           //a tetrahedron
//           if (!((tet_0==0 and tet_1==1 and tet_2==2 and tet_3==3)||
//                 (tet_0==0 and tet_1==1 and tet_2==4 and tet_3==5)||
//                 (tet_0==0 and tet_1==3 and tet_2==4 and tet_3==7)||
//                 (tet_0==0 and tet_1==1 and tet_2==6 and tet_3==7)||
//                 (tet_0==0 and tet_1==3 and tet_2==5 and tet_3==6)||
//                 (tet_0==0 and tet_1==2 and tet_2==4 and tet_3==6)||
//                 (tet_0==1 and tet_1==2 and tet_2==5 and tet_3==6)||
//                 (tet_0==1 and tet_1==2 and tet_2==4 and tet_3==7)||
//                 (tet_0==1 and tet_1==3 and tet_2==5 and tet_3==7)||
//                 (tet_0==2 and tet_1==3 and tet_2==6 and tet_3==7)||
//                 (tet_0==2 and tet_1==3 and tet_2==4 and tet_3==5)||
//                 (tet_0==4 and tet_1==5 and tet_2==6 and tet_3==7)))
//             ale3_add_tetra_stiffness(tet_2, tet_0, tet_1, tet_3, sys_mat, xyze);
//         }

  //Use 8 tetrahedra to prevent node-face-penetration
  ale3_add_tetra_stiffness(0,1,3,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(0,1,2,5,sys_mat,xyze);
  ale3_add_tetra_stiffness(1,2,3,6,sys_mat,xyze);
  ale3_add_tetra_stiffness(0,2,3,7,sys_mat,xyze);

  ale3_add_tetra_stiffness(0,4,5,7,sys_mat,xyze);
  ale3_add_tetra_stiffness(1,4,5,6,sys_mat,xyze);
  ale3_add_tetra_stiffness(2,5,6,7,sys_mat,xyze);
  ale3_add_tetra_stiffness(3,4,6,7,sys_mat,xyze);
}

void DRT::ELEMENTS::Ale3::static_ke_spring(Epetra_SerialDenseMatrix* sys_mat,
                                           vector<double>& displacements)
{
  const int iel = NumNode();                              // numnp to this element
  const DiscretizationType distype = this->Shape();
  int numcnd;                                             // number of corner nodes
  int node_i, node_j;                                     // end nodes of actual spring
  double length;                                          // length of actual edge
  double ex, ey, ez;                                      // direction of actual edge
  double factor;


  //number of corner nodes
  switch (distype)
  {
  case DRT::Element::tet4:
  case DRT::Element::tet10:
    numcnd = 4;
    break;
  case DRT::Element::pyramid5:
    numcnd = 5;
    break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
    numcnd = 6;
    break;
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
    numcnd = 8;
    break;
  default:
    numcnd = 0;
    dserror("distype unkown");
    break;
  }

  // get actual node coordinates
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0] + displacements[i*3];
    xyze(1,i)=Nodes()[i]->X()[1] + displacements[i*3+1];
    xyze(2,i)=Nodes()[i]->X()[2] + displacements[i*3+2];
  }

//lineal springs from all corner nodes to all corner nodes
//loop over all edges and diagonals of the element
  for (node_i=0; node_i<(numcnd-1); node_i++)
  {
    for (node_j=node_i+1; node_j<numcnd; node_j++)
    {
      ale3_edge_geometry(node_i,node_j,xyze,&length,&ex,&ey,&ez);
      factor = (1.0/length);
      //put values in 'element stiffness'
      //rows for node_i
      (*sys_mat)(node_i*3,  node_i*3  ) += (ex*ex * factor);
      (*sys_mat)(node_i*3,  node_i*3+1) += (ex*ey * factor);
      (*sys_mat)(node_i*3,  node_i*3+2) += (ex*ez * factor);

      (*sys_mat)(node_i*3,  node_j*3  ) -= (ex*ex * factor);
      (*sys_mat)(node_i*3,  node_j*3+1) -= (ex*ey * factor);
      (*sys_mat)(node_i*3,  node_j*3+2) -= (ex*ez * factor);

      (*sys_mat)(node_i*3+1,  node_i*3  ) += (ex*ey * factor);
      (*sys_mat)(node_i*3+1,  node_i*3+1) += (ey*ey * factor);
      (*sys_mat)(node_i*3+1,  node_i*3+2) += (ey*ez * factor);

      (*sys_mat)(node_i*3+1,  node_j*3  ) -= (ex*ey * factor);
      (*sys_mat)(node_i*3+1,  node_j*3+1) -= (ey*ey * factor);
      (*sys_mat)(node_i*3+1,  node_j*3+2) -= (ey*ez * factor);

      (*sys_mat)(node_i*3+2,  node_i*3  ) += (ex*ez * factor);
      (*sys_mat)(node_i*3+2,  node_i*3+1) += (ey*ez * factor);
      (*sys_mat)(node_i*3+2,  node_i*3+2) += (ez*ez * factor);

      (*sys_mat)(node_i*3+2,  node_j*3  ) -= (ex*ez * factor);
      (*sys_mat)(node_i*3+2,  node_j*3+1) -= (ey*ez * factor);
      (*sys_mat)(node_i*3+2,  node_j*3+2) -= (ez*ez * factor);

      //rows for node_j
      (*sys_mat)(node_j*3,  node_i*3  ) -= (ex*ex * factor);
      (*sys_mat)(node_j*3,  node_i*3+1) -= (ex*ey * factor);
      (*sys_mat)(node_j*3,  node_i*3+2) -= (ex*ez * factor);

      (*sys_mat)(node_j*3,  node_j*3  ) += (ex*ex * factor);
      (*sys_mat)(node_j*3,  node_j*3+1) += (ex*ey * factor);
      (*sys_mat)(node_j*3,  node_j*3+2) += (ex*ez * factor);

      (*sys_mat)(node_j*3+1,  node_i*3  ) -= (ex*ey * factor);
      (*sys_mat)(node_j*3+1,  node_i*3+1) -= (ey*ey * factor);
      (*sys_mat)(node_j*3+1,  node_i*3+2) -= (ey*ez * factor);

      (*sys_mat)(node_j*3+1,  node_j*3  ) += (ex*ey * factor);
      (*sys_mat)(node_j*3+1,  node_j*3+1) += (ey*ey * factor);
      (*sys_mat)(node_j*3+1,  node_j*3+2) += (ey*ez * factor);

      (*sys_mat)(node_j*3+2,  node_i*3  ) -= (ex*ez * factor);
      (*sys_mat)(node_j*3+2,  node_i*3+1) -= (ey*ez * factor);
      (*sys_mat)(node_j*3+2,  node_i*3+2) -= (ez*ez * factor);

      (*sys_mat)(node_j*3+2,  node_j*3  ) += (ex*ez * factor);
      (*sys_mat)(node_j*3+2,  node_j*3+1) += (ey*ez * factor);
      (*sys_mat)(node_j*3+2,  node_j*3+2) += (ez*ez * factor);
    }
  }

  //build in torsional springs
  //and put edge nodes on the middle of the respective edge
  switch (distype)
  {
  case DRT::Element::tet10:

    (*sys_mat)(12,0)  = -0.5;
    (*sys_mat)(12,3)  = -0.5;
    (*sys_mat)(12,12) =  1.0;
    (*sys_mat)(13,1)  = -0.5;
    (*sys_mat)(13,4)  = -0.5;
    (*sys_mat)(13,13) =  1.0;
    (*sys_mat)(14,2)  = -0.5;
    (*sys_mat)(14,5)  = -0.5;
    (*sys_mat)(14,14) =  1.0;

    (*sys_mat)(15,3)  = -0.5;
    (*sys_mat)(15,6)  = -0.5;
    (*sys_mat)(15,15) =  1.0;
    (*sys_mat)(16,4)  = -0.5;
    (*sys_mat)(16,7)  = -0.5;
    (*sys_mat)(16,16) =  1.0;
    (*sys_mat)(17,5)  = -0.5;
    (*sys_mat)(17,8)  = -0.5;
    (*sys_mat)(17,17) =  1.0;

    (*sys_mat)(18,6)  = -0.5;
    (*sys_mat)(18,0)  = -0.5;
    (*sys_mat)(18,18) =  1.0;
    (*sys_mat)(19,7)  = -0.5;
    (*sys_mat)(19,1)  = -0.5;
    (*sys_mat)(19,19) =  1.0;
    (*sys_mat)(20,8)  = -0.5;
    (*sys_mat)(20,2)  = -0.5;
    (*sys_mat)(20,20) =  1.0;

    (*sys_mat)(21,0)  = -0.5;
    (*sys_mat)(21,9)  = -0.5;
    (*sys_mat)(21,21) =  1.0;
    (*sys_mat)(22,1)  = -0.5;
    (*sys_mat)(22,10) = -0.5;
    (*sys_mat)(22,22) =  1.0;
    (*sys_mat)(23,2)  = -0.5;
    (*sys_mat)(23,11) = -0.5;
    (*sys_mat)(23,23) =  1.0;

    (*sys_mat)(24,3)  = -0.5;
    (*sys_mat)(24,9)  = -0.5;
    (*sys_mat)(24,24) =  1.0;
    (*sys_mat)(25,4)  = -0.5;
    (*sys_mat)(25,10) = -0.5;
    (*sys_mat)(25,25) =  1.0;
    (*sys_mat)(26,5)  = -0.5;
    (*sys_mat)(26,11) = -0.5;
    (*sys_mat)(26,26) =  1.0;

    (*sys_mat)(27,6)  = -0.5;
    (*sys_mat)(27,9)  = -0.5;
    (*sys_mat)(27,27) =  1.0;
    (*sys_mat)(28,7)  = -0.5;
    (*sys_mat)(28,10) = -0.5;
    (*sys_mat)(28,28) =  1.0;
    (*sys_mat)(29,8)  = -0.5;
    (*sys_mat)(29,11) = -0.5;
    (*sys_mat)(29,29) =  1.0;

    ale3_tors_spring_tet4(sys_mat,xyze);
    break;

  case DRT::Element::tet4:
    ale3_tors_spring_tet4(sys_mat,xyze);
    break;

  case DRT::Element::pyramid5:
    ale3_tors_spring_pyramid5(sys_mat,xyze);
    break;

  case DRT::Element::wedge15:

    for (int k=0; k<3; k++)
    {
      (*sys_mat)(3*(6+k)   , 3*k    ) = -0.5;
      (*sys_mat)(3*(6+k) +1, 3*k +1 ) = -0.5;
      (*sys_mat)(3*(6+k) +2, 3*k +2 ) = -0.5;
      (*sys_mat)(3*(6+k)   , 3*((k+1)*((int)(k<2)))    ) = -0.5;
      (*sys_mat)(3*(6+k) +1, 3*((k+1)*((int)(k<2))) +1 ) = -0.5;
      (*sys_mat)(3*(6+k) +2, 3*((k+1)*((int)(k<2))) +2 ) = -0.5;
      (*sys_mat)(3*(6+k)   , 3*(6+k)   ) =  1.0;
      (*sys_mat)(3*(6+k) +1, 3*(6+k) +1) =  1.0;
      (*sys_mat)(3*(6+k) +2, 3*(6+k) +2) =  1.0;
    }

    for (int k=0; k<3; k++)
    {
      (*sys_mat)(3*(9+k)  , 3*k   ) = -0.5;
      (*sys_mat)(3*(9+k)+1, 3*k+1 ) = -0.5;
      (*sys_mat)(3*(9+k)+2, 3*k+2 ) = -0.5;
      (*sys_mat)(3*(9+k)  , 3*(3+k)   ) = -0.5;
      (*sys_mat)(3*(9+k)+1, 3*(3+k)+1 ) = -0.5;
      (*sys_mat)(3*(9+k)+2, 3*(3+k)+2 ) = -0.5;
      (*sys_mat)(3*(9+k)  , 3*(9+k)  ) =  1.0;
      (*sys_mat)(3*(9+k)+1, 3*(9+k)+1) =  1.0;
      (*sys_mat)(3*(9+k)+2, 3*(9+k)+2) =  1.0;
    }

    for (int k=0; k<3; k++)
    {
      (*sys_mat)(3*(12+k)   , 3*(3+k)    ) = -0.5;
      (*sys_mat)(3*(12+k) +1, 3*(3+k) +1 ) = -0.5;
      (*sys_mat)(3*(12+k) +2, 3*(3+k) +2 ) = -0.5;
      (*sys_mat)(3*(12+k)   , 3*(3+(k+1)*((int)(k<2)))    ) = -0.5;
      (*sys_mat)(3*(12+k) +1, 3*(3+(k+1)*((int)(k<2))) +1 ) = -0.5;
      (*sys_mat)(3*(12+k) +2, 3*(3+(k+1)*((int)(k<2))) +2 ) = -0.5;
      (*sys_mat)(3*(12+k)   , 3*(12+k)   ) =  1.0;
      (*sys_mat)(3*(12+k) +1, 3*(12+k) +1) =  1.0;
      (*sys_mat)(3*(12+k) +2, 3*(12+k) +2) =  1.0;
    }

    ale3_tors_spring_wedge6(sys_mat,xyze);
    break;

  case DRT::Element::wedge6:
    ale3_tors_spring_wedge6(sys_mat,xyze);
    break;


  case DRT::Element::hex20:

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(8+k)   , 3*k    ) = -0.5;
      (*sys_mat)(3*(8+k) +1, 3*k +1 ) = -0.5;
      (*sys_mat)(3*(8+k) +2, 3*k +2 ) = -0.5;
      (*sys_mat)(3*(8+k)   , 3*((k+1)*((int)(k<3)))    ) = -0.5;
      (*sys_mat)(3*(8+k) +1, 3*((k+1)*((int)(k<3))) +1 ) = -0.5;
      (*sys_mat)(3*(8+k) +2, 3*((k+1)*((int)(k<3))) +2 ) = -0.5;
      (*sys_mat)(3*(8+k)   , 3*(8+k)   ) =  1.0;
      (*sys_mat)(3*(8+k) +1, 3*(8+k) +1) =  1.0;
      (*sys_mat)(3*(8+k) +2, 3*(8+k) +2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(12+k)  , 3*k   ) = -0.5;
      (*sys_mat)(3*(12+k)+1, 3*k+1 ) = -0.5;
      (*sys_mat)(3*(12+k)+2, 3*k+2 ) = -0.5;
      (*sys_mat)(3*(12+k)  , 3*(4+k)   ) = -0.5;
      (*sys_mat)(3*(12+k)+1, 3*(4+k)+1 ) = -0.5;
      (*sys_mat)(3*(12+k)+2, 3*(4+k)+2 ) = -0.5;
      (*sys_mat)(3*(12+k)  , 3*(12+k)  ) =  1.0;
      (*sys_mat)(3*(12+k)+1, 3*(12+k)+1) =  1.0;
      (*sys_mat)(3*(12+k)+2, 3*(12+k)+2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(16+k)   , 3*(4+k)    ) = -0.5;
      (*sys_mat)(3*(16+k) +1, 3*(4+k) +1 ) = -0.5;
      (*sys_mat)(3*(16+k) +2, 3*(4+k) +2 ) = -0.5;
      (*sys_mat)(3*(16+k)   , 3*(4+(k+1)*((int)(k<3)))    ) = -0.5;
      (*sys_mat)(3*(16+k) +1, 3*(4+(k+1)*((int)(k<3))) +1 ) = -0.5;
      (*sys_mat)(3*(16+k) +2, 3*(4+(k+1)*((int)(k<3))) +2 ) = -0.5;
      (*sys_mat)(3*(16+k)   , 3*(16+k)   ) =  1.0;
      (*sys_mat)(3*(16+k) +1, 3*(16+k) +1) =  1.0;
      (*sys_mat)(3*(16+k) +2, 3*(16+k) +2) =  1.0;
    }

    ale3_tors_spring_hex8(sys_mat,xyze);
    break;

  case DRT::Element::hex27:

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(8+k)   , 3*k    ) = -0.5;
      (*sys_mat)(3*(8+k) +1, 3*k +1 ) = -0.5;
      (*sys_mat)(3*(8+k) +2, 3*k +2 ) = -0.5;
      (*sys_mat)(3*(8+k)   , 3*((k+1)*((int)(k<3)))    ) = -0.5;
      (*sys_mat)(3*(8+k) +1, 3*((k+1)*((int)(k<3))) +1 ) = -0.5;
      (*sys_mat)(3*(8+k) +2, 3*((k+1)*((int)(k<3))) +2 ) = -0.5;
      (*sys_mat)(3*(8+k)   , 3*(8+k)   ) =  1.0;
      (*sys_mat)(3*(8+k) +1, 3*(8+k) +1) =  1.0;
      (*sys_mat)(3*(8+k) +2, 3*(8+k) +2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(12+k)  , 3*k   ) = -0.5;
      (*sys_mat)(3*(12+k)+1, 3*k+1 ) = -0.5;
      (*sys_mat)(3*(12+k)+2, 3*k+2 ) = -0.5;
      (*sys_mat)(3*(12+k)  , 3*(4+k)   ) = -0.5;
      (*sys_mat)(3*(12+k)+1, 3*(4+k)+1 ) = -0.5;
      (*sys_mat)(3*(12+k)+2, 3*(4+k)+2 ) = -0.5;
      (*sys_mat)(3*(12+k)  , 3*(12+k)  ) =  1.0;
      (*sys_mat)(3*(12+k)+1, 3*(12+k)+1) =  1.0;
      (*sys_mat)(3*(12+k)+2, 3*(12+k)+2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(16+k)   , 3*(4+k)    ) = -0.5;
      (*sys_mat)(3*(16+k) +1, 3*(4+k) +1 ) = -0.5;
      (*sys_mat)(3*(16+k) +2, 3*(4+k) +2 ) = -0.5;
      (*sys_mat)(3*(16+k)   , 3*(4+(k+1)*((int)(k<3)))    ) = -0.5;
      (*sys_mat)(3*(16+k) +1, 3*(4+(k+1)*((int)(k<3))) +1 ) = -0.5;
      (*sys_mat)(3*(16+k) +2, 3*(4+(k+1)*((int)(k<3))) +2 ) = -0.5;
      (*sys_mat)(3*(16+k)   , 3*(16+k)   ) =  1.0;
      (*sys_mat)(3*(16+k) +1, 3*(16+k) +1) =  1.0;
      (*sys_mat)(3*(16+k) +2, 3*(16+k) +2) =  1.0;
    }

    (*sys_mat)(3*20  , 3*8   ) = -0.5;
    (*sys_mat)(3*20+1, 3*8+1 ) = -0.5;
    (*sys_mat)(3*20+2, 3*8+2 ) = -0.5;
    (*sys_mat)(3*20  , 3*10  ) = -0.5;
    (*sys_mat)(3*20+1, 3*10+1) = -0.5;
    (*sys_mat)(3*20+2, 3*10+2) = -0.5;
    (*sys_mat)(3*20  , 3*20  ) =  1.0;
    (*sys_mat)(3*20+1, 3*20+1) =  1.0;
    (*sys_mat)(3*20+2, 3*20+2) =  1.0;

    for (int k=0; k<4; k++)
    {
      (*sys_mat)(3*(21+k)   , 3*(8+k)    ) = -0.5;
      (*sys_mat)(3*(21+k) +1, 3*(8+k) + 1) = -0.5;
      (*sys_mat)(3*(21+k) +2, 3*(8+k) + 2) = -0.5;
      (*sys_mat)(3*(21+k)   , 3*(16+k)   ) = -0.5;
      (*sys_mat)(3*(21+k) +1, 3*(16+k)+ 1) = -0.5;
      (*sys_mat)(3*(21+k) +2, 3*(16+k)+ 2) = -0.5;
      (*sys_mat)(3*(21+k)   , 3*(21+k)   ) =  1.0;
      (*sys_mat)(3*(21+k) +1, 3*(21+k) +1) =  1.0;
      (*sys_mat)(3*(21+k) +2, 3*(21+k) +2) =  1.0;
    }

    (*sys_mat)(3*25  , 3*16  ) = -0.5;
    (*sys_mat)(3*25+1, 3*16+1) = -0.5;
    (*sys_mat)(3*25+2, 3*16+2) = -0.5;
    (*sys_mat)(3*25  , 3*18  ) = -0.5;
    (*sys_mat)(3*25+1, 3*18+1) = -0.5;
    (*sys_mat)(3*25+2, 3*18+2) = -0.5;
    (*sys_mat)(3*25  , 3*25  ) =  1.0;
    (*sys_mat)(3*25+1, 3*25+1) =  1.0;
    (*sys_mat)(3*25+2, 3*25+2) =  1.0;

    (*sys_mat)(3*26  , 3*21  ) = -0.5;
    (*sys_mat)(3*26+1, 3*21+1) = -0.5;
    (*sys_mat)(3*26+2, 3*21+2) = -0.5;
    (*sys_mat)(3*26  , 3*23  ) = -0.5;
    (*sys_mat)(3*26+1, 3*23+1) = -0.5;
    (*sys_mat)(3*26+2, 3*23+2) = -0.5;
    (*sys_mat)(3*26  , 3*26  ) =  1.0;
    (*sys_mat)(3*26+1, 3*26+1) =  1.0;
    (*sys_mat)(3*26+2, 3*26+2) =  1.0;

    ale3_tors_spring_hex8(sys_mat,xyze);
    break;

  case DRT::Element::hex8:
    ale3_tors_spring_hex8(sys_mat,xyze);
    break;

  default:
    dserror("unknown distype in ale spring dynamic");
    break;
  }


}

void DRT::ELEMENTS::Ale3::static_ke(vector<int>&              lm,
                                    Epetra_SerialDenseMatrix* sys_mat,
                                    Epetra_SerialDenseVector* residual,
				    RefCountPtr<MAT::Material> material,
                                    ParameterList&            params)
{
  const int iel = NumNode();
  const int nd  = 3 * iel;
  const DiscretizationType distype = this->Shape();


  //  get material using class StVenantKirchhoff
  if (material->MaterialType()!=m_stvenant)
    dserror("stvenant material expected but got type %d", material->MaterialType());
  MAT::StVenantKirchhoff* actmat = static_cast<MAT::StVenantKirchhoff*>(material.get());

  Epetra_SerialDenseMatrix xyze(3,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
    xyze(2,i)=Nodes()[i]->X()[2];
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix 	deriv(3,iel);
  Epetra_SerialDenseMatrix 	xjm(3,3);
  Epetra_SerialDenseMatrix 	bop(6,3*iel);
  Epetra_SerialDenseMatrix 	D(6,6);

  double                        vol=0.;

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
        const double e1 = intpoints.qxg[iquad][0];
        const double e2 = intpoints.qxg[iquad][1];
        const double e3 = intpoints.qxg[iquad][2];
        // shape functions and their derivatives
        DRT::UTILS::shape_function_3D(funct,e1,e2,e3,distype);
        DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

        // compute jacobian matrix

        // determine jacobian at point r,s,t
        for (int i=0; i<3; i++)
        {
          for (int j=0; j<3; j++)
          {
            double dum=0.;
            for (int l=0; l<iel; l++)
            {
              dum += deriv(i,l)*xyze(j,l);
            }
            xjm(i,j)=dum;
          }
        }

        // determinant of jacobian
        double det = xjm(0,0)*xjm(1,1)*xjm(2,2) +
                     xjm(0,1)*xjm(1,2)*xjm(2,0) +
                     xjm(0,2)*xjm(1,0)*xjm(2,1) -
                     xjm(0,2)*xjm(1,1)*xjm(2,0) -
                     xjm(0,0)*xjm(1,2)*xjm(2,1) -
                     xjm(0,1)*xjm(1,0)*xjm(2,2);


        // calculate element volume
        const double fac = intpoints.qwgt[iquad]*det;
        vol += fac;

        // calculate operator B

        // inverse of jacobian
        double x1r = xjm(0,0);
        double x2r = xjm(0,1);
        double x3r = xjm(0,2);
        double x1s = xjm(1,0);
        double x2s = xjm(1,1);
        double x3s = xjm(1,2);
        double x1t = xjm(2,0);
        double x2t = xjm(2,1);
        double x3t = xjm(2,2);

        const double dum=1.0/det;

        double xi11=dum*(x2s*x3t - x2t*x3s);
        double xi12=dum*(x3r*x2t - x2r*x3t);
        double xi13=dum*(x2r*x3s - x3r*x2s);
        double xi21=dum*(x3s*x1t - x3t*x1s);
        double xi22=dum*(x1r*x3t - x3r*x1t);
        double xi23=dum*(x3r*x1s - x1r*x3s);
        double xi31=dum*(x1s*x2t - x1t*x2s);
        double xi32=dum*(x2r*x1t - x1r*x2t);
        double xi33=dum*(x1r*x2s - x2r*x1s);

        // get operator b of global derivatives
        for (int i=0; i<iel; i++)
        {
          const int node_start = i*3;

          const double hr   = deriv(0,i);
          const double hs   = deriv(1,i);
          const double ht   = deriv(2,i);

          const double h1 = xi11*hr + xi12*hs + xi13*ht;
          const double h2 = xi21*hr + xi22*hs + xi23*ht;
          const double h3 = xi31*hr + xi32*hs + xi33*ht;

          bop(0,node_start+0) = h1 ;
          bop(0,node_start+1) = 0.0;
          bop(0,node_start+2) = 0.0;
          bop(1,node_start+0) = 0.0;
          bop(1,node_start+1) = h2 ;
          bop(1,node_start+2) = 0.0;
          bop(2,node_start+0) = 0.0;
          bop(2,node_start+1) = 0.0;
          bop(2,node_start+2) = h3 ;
          bop(3,node_start+0) = h2 ;
          bop(3,node_start+1) = h1 ;
          bop(3,node_start+2) = 0.0;
          bop(4,node_start+0) = 0.0;
          bop(4,node_start+1) = h3 ;
          bop(4,node_start+2) = h2 ;
          bop(5,node_start+0) = h3 ;
          bop(5,node_start+1) = 0.0;
          bop(5,node_start+2) = h1 ;
        }

        // call material law
	actmat->SetupCmat(&D);

        // elastic stiffness matrix ke
        //ale3_keku(estif,bop,D,fac,nd);

        for (int j=0; j<nd; j++)
        {
          double db[6];
          for (int k=0; k<6; k++)
          {
            db[k] = 0.0;
            for (int l=0; l<6; l++)
            {
              db[k] += D(k,l)*bop(l,j)*fac ;
            }
          }
          for (int i=0; i<nd; i++)
          {
            double dum = 0.0;
            for (int m=0; m<6; m++)
            {
              dum = dum + bop(m,i)*db[m] ;
            }
            (*sys_mat)(i,j) += dum ;
          }
        }

        // hourglass stabalization stiffness matrix ke
        /*
          see also:
          (1) T. Belytschko and L.P. Bindeman:
              Assumed strain stabilization of the 8-node hexahedral element
              Comp. Meth. Appl. Mech. Eng.: 105 (1993) p. 225-260.
          (2) D.P. Flanagan and T. Belytschko:
              A uniform strain hexahedron and quadrilateral with orthogonal
              hourglass control
              Int. J. Num. Meth. Ing.: Vol. 17 (1981) p. 679-706.
         */

	//Integration rule for hour-glass-stabilization. Not used in the moment. If needed,
	//it should be implemented within the getOptimalGaussrule-method
#if 0
        if (distype==hex8 and intpoints.nquad == 1)
        {
          const double ee = material->m.stvenant->youngs;
          const double nu = material->m.stvenant->possionratio;
          const double mu = ee / (2*(1+nu));

          // ASQBI
          const double c1 = 1.0/(1.0 - nu);
          const double c2 = (1.0 + nu)/3;
          const double c3 = 1.0/(1.0 - nu);

          double         xc[3][8];
          double         b[3][8];

          double         a[3][3];
          double         ba[3];
          double         r[3][3];

          double         h[4][8] = {{1,1,-1,-1,-1,-1,1,1},{1,-1,-1,1,-1,1,1,-1},
                                    {1,-1,1,-1,1,-1,1,-1},{-1,1,-1,1,1,-1,1,-1}};
          double         lam[3][8] = {{-1,1,1,-1,-1,1,1,-1},
                                      {-1,-1,1,1,-1,-1,1,1},
                                      {-1,-1,-1,-1,1,1,1,1}};
          double         gam[4][8];
          double         lx[3];
          double         hh[3][3];

          double         gg00[8][8], gg11[8][8], gg22[8][8], gg33[8][8];
          double         gg01[8][8], gg10[8][8], gg02[8][8], gg20[8][8];
          double         gg12[8][8], gg21[8][8];

          double         kstab[24][24];

          // corotational coordinate system: rotation tensor r[3][3]
          // accord. to (1) Appendix A, Table 9
          for (int i=0; i<2; i++)
          {
            for (int j=0; j<3; j++)
            {
              a[i][j] = 0.0;
              for (int k=0; k<8; k++)
              {
                a[i][j] += lam[i][k]*xyze(j,k);
              }
            }
          }

          const double dum =(a[0][0]*a[1][0]+a[0][1]*a[1][1]+a[0][2]*a[1][2])/
                            (a[0][0]*a[0][0]+a[0][1]*a[0][1]+a[0][2]*a[0][2]);

          a[1][0] = a[1][0] - dum * a[0][0];
          a[1][1] = a[1][1] - dum * a[0][1];
          a[1][2] = a[1][2] - dum * a[0][2];

          a[2][0] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
          a[2][1] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
          a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

          for(int i = 0; i<3; i++)
          {
            ba[i] = sqrt(a[i][0]*a[i][0]+a[i][1]*a[i][1]+a[i][2]*a[i][2]);
            r[i][0] = a[i][0] / ba[i];
            r[i][1] = a[i][1] / ba[i];
            r[i][2] = a[i][2] / ba[i];
          }

          // transforming nodal coordinates to the corotational system
          for (int i=0; i<8; i++)
          {
            for (int j=0; j<3; j++)
            {
              xc[j][i] = r[j][0]*xyze(0,i)+r[j][1]*xyze(1,i)+r[j][2]*xyze(2,i);
            }
          }

          int l0=0,l1=0,l2=0,l3=0,l4=0,l5=0,l6=0,l7=0;

          // B-matrix b[3][8] accord. to (2) Appendix I, eqn (79)
          for (int i=0; i<3; i++)
          {
            const int j = (i+1)%3;
            const int k = (j+1)%3;
            for (int l=0; l<8;l++)
            {
              switch (l)
              {
              case 0:
                l0=0;l1=1;l2=2;l3=3;l4=4;l5=5;l6=6;l7=7;
                break;
              case 1:
                l0=1;l1=2;l2=3;l3=0;l4=5;l5=6;l6=7;l7=4;
                break;
              case 2:
                l0=2;l1=3;l2=0;l3=1;l4=6;l5=7;l6=4;l7=5;
                break;
              case 3:
                l0=3;l1=0;l2=1;l3=2;l4=7;l5=4;l6=5;l7=6;
                break;
              case 4:
                l0=4;l1=7;l2=6;l3=5;l4=0;l5=3;l6=2;l7=1;
                break;
              case 5:
                l0=5;l1=4;l2=7;l3=6;l4=1;l5=0;l6=3;l7=2;
                break;
              case 6:
                l0=6;l1=5;l2=4;l3=7;l4=2;l5=1;l6=0;l7=3;
                break;
              case 7:
                l0=7;l1=6;l2=5;l3=4;l4=3;l5=2;l6=1;l7=0;
                break;
              }
              b[i][l0] =1/(12 * vol) * (xc[j][l1]*((xc[k][l5] - xc[k][l2]) -
                                                   (xc[k][l3] - xc[k][l4])) +
                                        xc[j][l2]* (xc[k][l1] - xc[k][l3]) +
                                        xc[j][l3]*((xc[k][l2] - xc[k][l7]) -
                                                   (xc[k][l4] - xc[k][l1])) +
                                        xc[j][l4]*((xc[k][l7] - xc[k][l5]) -
                                                   (xc[k][l1] - xc[k][l3])) +
                                        xc[j][l5]* (xc[k][l4] - xc[k][l1]) +
                                        xc[j][l7]* (xc[k][l3] - xc[k][l4]) );
            }
          }

          // gamma vectors, accord. to (1) eqn (2.12b)
          for (int i=0; i<4; i++)
          {
            for (int j=0; j<8; j++)
            {
              gam[i][j] = 0.125 * h[i][j];
              for (int k=0; k<3; k++)
              {
                const double dum = h[i][0]*xc[k][0]+h[i][1]*xc[k][1]+h[i][2]*xc[k][2]+
                                   h[i][3]*xc[k][3]+h[i][4]*xc[k][4]+h[i][5]*xc[k][5]+
                                   h[i][6]*xc[k][6]+h[i][7]*xc[k][7];
                gam[i][j] -= 0.125 * dum * b[k][j];
              }
            }
          }

          /* lambda * x (auxiliary vector) */
          lx[0] = lam[0][0]*xc[0][0]+lam[0][1]*xc[0][1]+lam[0][2]*xc[0][2]+
                  lam[0][3]*xc[0][3]+lam[0][4]*xc[0][4]+lam[0][5]*xc[0][5]+
                  lam[0][6]*xc[0][6]+lam[0][7]*xc[0][7];
          lx[1] = lam[1][0]*xc[1][0]+lam[1][1]*xc[1][1]+lam[1][2]*xc[1][2]+
                  lam[1][3]*xc[1][3]+lam[1][4]*xc[1][4]+lam[1][5]*xc[1][5]+
                  lam[1][6]*xc[1][6]+lam[1][7]*xc[1][7];
          lx[2] = lam[2][0]*xc[2][0]+lam[2][1]*xc[2][1]+lam[2][2]*xc[2][2]+
                  lam[2][3]*xc[2][3]+lam[2][4]*xc[2][4]+lam[2][5]*xc[2][5]+
                  lam[2][6]*xc[2][6]+lam[2][7]*xc[2][7];

          /* H_ij, accord. to (1) eqns. (3.15d) and (3.15e) */
          hh[0][0] = 1.0/3.0 * (lx[1]*lx[2])/lx[0];
          hh[1][1] = 1.0/3.0 * (lx[2]*lx[0])/lx[1];
          hh[2][2] = 1.0/3.0 * (lx[0]*lx[1])/lx[2];
          hh[0][1] = 1.0/3.0 * lx[2];
          hh[1][0] = 1.0/3.0 * lx[2];
          hh[0][2] = 1.0/3.0 * lx[1];
          hh[2][0] = 1.0/3.0 * lx[1];
          hh[1][2] = 1.0/3.0 * lx[0];
          hh[2][1] = 1.0/3.0 * lx[0];

          /* stabalization matrix with respect to the corotational ccord. system. */
          /* rearranging orders of dofs, accord. to (1) eqns. (3.15a) to (3.15c) */
          for (int i=0; i<8; i++)
          {
            for (int j=0; j<8; j++)
            {
              gg00[i][j] = gam[0][i] * gam[0][j];
              gg11[i][j] = gam[1][i] * gam[1][j];
              gg22[i][j] = gam[2][i] * gam[2][j];
              gg33[i][j] = gam[3][i] * gam[3][j];
              gg01[i][j] = gam[0][i] * gam[1][j];
              gg10[i][j] = gam[1][i] * gam[0][j];
              gg02[i][j] = gam[0][i] * gam[2][j];
              gg20[i][j] = gam[2][i] * gam[0][j];
              gg12[i][j] = gam[1][i] * gam[2][j];
              gg21[i][j] = gam[2][i] * gam[1][j];

              /* kstab 00 */
              kstab[i*3][j*3]     = 2*mu* (hh[0][0]*(c1*(gg11[i][j] + gg22[i][j])
                                                     + c2*gg33[i][j]) + 0.5 * (hh[1][1] + hh[2][2]) * gg00[i][j]);

              /* kstab 11 */
              kstab[i*3+1][j*3+1] = 2*mu* (hh[1][1]*(c1*(gg22[i][j] + gg00[i][j])
                                                     + c2*gg33[i][j]) + 0.5 * (hh[2][2] + hh[0][0]) * gg11[i][j]);

              /* kstab 22 */
              kstab[i*3+2][j*3+2] = 2*mu* (hh[2][2]*(c1*(gg00[i][j] + gg11[i][j])
                                                     + c2*gg33[i][j]) + 0.5 * (hh[0][0] + hh[1][1]) * gg22[i][j]);

              /* kstab 01 */
              kstab[i*3][j*3+1]   = 2*mu* (hh[0][1]*(c3*gg10[i][j]+0.5*gg01[i][j]));

              /* kstab 10 */
              kstab[i*3+1][j*3]   = 2*mu* (hh[1][0]*(c3*gg01[i][j]+0.5*gg10[i][j]));

              /* kstab 02 */
              kstab[i*3][j*3+2]   = 2*mu* (hh[0][2]*(c3*gg20[i][j]+0.5*gg02[i][j]));

              /* kstab 20 */
              kstab[i*3+2][j*3]   = 2*mu* (hh[2][0]*(c3*gg02[i][j]+0.5*gg20[i][j]));

              /* kstab 12 */
              kstab[i*3+1][j*3+2] = 2*mu* (hh[1][2]*(c3*gg21[i][j]+0.5*gg12[i][j]));

              /* kstab 21 */
              kstab[i*3+2][j*3+1] = 2*mu* (hh[2][1]*(c3*gg12[i][j]+0.5*gg21[i][j]));
            }
          }

          /* transforming kstab to the global coordinate system and */
          /* and adding to the one point quadrature matrix */
          for (int i=0; i<8; i++)
          {
            for (int j=0;j<8;j++)
            {
              for (int k=0;k<3;k++)
              {
                for (int l=0;l<3;l++)
                {
                  (*sys_mat)(i*3+k,j*3+l) += (r[0][k]*kstab[i*3+0][j*3+0] +
                                              r[1][k]*kstab[i*3+1][j*3+0] +
                                              r[2][k]*kstab[i*3+2][j*3+0]) * r[0][l] +
                                             (r[0][k]*kstab[i*3+0][j*3+1] +
                                              r[1][k]*kstab[i*3+1][j*3+1] +
                                              r[2][k]*kstab[i*3+2][j*3+1]) * r[1][l] +
                                             (r[0][k]*kstab[i*3+0][j*3+2] +
                                              r[1][k]*kstab[i*3+1][j*3+2] +
                                              r[2][k]*kstab[i*3+2][j*3+2]) * r[2][l];
                }
              }
            }
          }
        }
#endif
  }

}



// get optimal gaussrule for discretization type
GaussRule3D DRT::ELEMENTS::Ale3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule = intrule3D_undefined;
    switch (distype)
    {
    case hex8:
        rule = intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = intrule_hex_27point;
        break;
    case tet4:
        rule = intrule_tet_4point;
        break;
    case tet10:
        rule = intrule_tet_5point;
        break;
    case wedge6:
        rule = intrule_wedge_6point;
        break;
    case wedge15:
        rule = intrule_wedge_9point;
        break;
    case pyramid5:
      rule = intrule_pyramid_8point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
    }
    return rule;
}


//=======================================================================
//=======================================================================

int DRT::ELEMENTS::Ale3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif
#endif
