//-----------------------------------------------------------------------
/*!
\file ale3_evaluate.cpp

<pre>

</pre>
*/
//-----------------------------------------------------------------------
#ifdef D_ALE



#include "ale3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_mat/stvenantkirchhoff.H"

using namespace DRT::UTILS;


DRT::ELEMENTS::Ale3_Impl_Interface* DRT::ELEMENTS::Ale3_Impl_Interface::Impl(DRT::ELEMENTS::Ale3* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    return Ale3_Impl<DRT::Element::hex8>::Instance(true);
  }
  case DRT::Element::hex20:
  {
    return Ale3_Impl<DRT::Element::hex20>::Instance(true);
  }
  case DRT::Element::hex27:
  {
    return Ale3_Impl<DRT::Element::hex27>::Instance(true);
  }
  case DRT::Element::tet4:
  {
    return Ale3_Impl<DRT::Element::tet4>::Instance(true);
  }
  case DRT::Element::tet10:
  {
    return Ale3_Impl<DRT::Element::tet10>::Instance(true);
  }
  case DRT::Element::wedge6:
  {
    return Ale3_Impl<DRT::Element::wedge6>::Instance(true);
  }
/*  case DRT::Element::wedge15:
  {
    return Ale3_Impl<DRT::Element::wedge15>::Instance(true);
  }*/
  case DRT::Element::pyramid5:
  {
    return Ale3_Impl<DRT::Element::pyramid5>::Instance(true);
  }
  case DRT::Element::nurbs8:
  {
    return Ale3_Impl<DRT::Element::nurbs8>::Instance(true);
  }
  case DRT::Element::nurbs27:
  {
    return Ale3_Impl<DRT::Element::nurbs27>::Instance(true);
  }
  default:
    dserror("shape %d (%d nodes) not supported", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Ale3_Impl<distype> * DRT::ELEMENTS::Ale3_Impl<distype>::Instance( bool create )
{
  static Ale3_Impl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
      instance = new Ale3_Impl<distype>();
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Ale3_Impl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


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
  else if (action == "calc_ale_laplace")
      act = Ale3::calc_ale_laplace;
  else if (action == "calc_ale_lin_stiff")
    act = Ale3::calc_ale_lin_stiff;
  else if (action == "calc_ale_spring")
    act = Ale3::calc_ale_spring;
  else if (action == "calc_ale_spring_fixed_ref")
    act = Ale3::calc_ale_spring_fixed_ref;
  else if (action == "calc_ale_node_normal")
    act = Ale3::calc_ale_node_normal;
  else
    dserror("Unknown type of action for Ale3");

  // get the material
  RefCountPtr<MAT::Material> mat = Material();

  switch (act)
  {
  case calc_ale_laplace:
  {
    std::vector<double> my_dispnp;
    bool incremental = params.get<bool>("incremental");
    if (incremental)
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      my_dispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);
    }

    Ale3_Impl_Interface::Impl(this)->static_ke_laplace(
                           this,
                           discretization,
                           elemat1,
                           elevec1,
                           incremental,
                           my_dispnp,
                           mat,
                           params);

    break;
  }
  case calc_ale_lin_stiff:
  {
    std::vector<double> my_dispnp;
    bool incremental = params.get<bool>("incremental");
    if (incremental)
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      my_dispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);
    }

    Ale3_Impl_Interface::Impl(this)->static_ke(this,
					       discretization,
					       lm,
					       elemat1,
					       elevec1,
					       incremental,
					       my_dispnp,
					       mat,
					       params);

    break;
  }

  case calc_ale_spring:
  {
    RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    vector<double> my_dispnp(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

    Ale3_Impl_Interface::Impl(this)->static_ke_spring(this,elemat1,my_dispnp);

    break;
  }

  case calc_ale_spring_fixed_ref:
  {
    vector<double> my_dispnp(lm.size(),0.0);
    Ale3_Impl_Interface::Impl(this)->static_ke_spring(this,elemat1,my_dispnp);

    break;
  }

  case calc_ale_node_normal:
  {
    RefCountPtr<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
    vector<double> my_dispnp(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

    Ale3_Impl_Interface::Impl(this)->ElementNodeNormal(this,elevec1,my_dispnp);

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
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

//=================== ElementNodeNormal =============================
// Calculate node normals acc. to Wall (7.13)
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ElementNodeNormal(
  Ale3*                      ele,
  Epetra_SerialDenseVector& elevec1,
  std::vector<double>&       my_dispnp
  )
{
  if (distype == DRT::Element::nurbs8 or
      distype == DRT::Element::nurbs27)
  {
    dserror("not implemented!");
  }

  LINALG::Matrix<3,iel> xyze;

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for(int i=0;i<iel;i++)
  {
    const double* x = nodes[i]->X();
    xyze(0,i)=x[0];
    xyze(1,i)=x[1];
    xyze(2,i)=x[2];
  }

  for(int i=0;i<iel;i++)
  {
    xyze(0,i) += my_dispnp[3*i+0];
    xyze(1,i) += my_dispnp[3*i+1];
    xyze(2,i) += my_dispnp[3*i+2];
  }

  /*----------------------------------------- declaration of variables ---*/
  LINALG::Matrix<iel,1  > funct;
  LINALG::Matrix<3,  iel> deriv;
  LINALG::Matrix<3,  3  > xjm;
  LINALG::Matrix<3,  3  > xji;

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule();
  const IntegrationPoints3D  intpoints(gaussrule);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];

    // get values of shape functions and derivatives in the gausspoint
    DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

    // compute jacobian matrix
    // determine jacobian at point r,s,t
    xjm.MultiplyNT(deriv,xyze);

    // determinant and inverse of jacobian
    const double det = xji.Invert(xjm);

    // integrate shapefunction gradient over element
    const double fac = intpoints.qwgt[iquad]*det;

    for (int node=0; node<iel; ++node)
    {
      for (int dim=0; dim<3; ++dim)
      {
        int row = 3 * node + dim;
        elevec1(row) += (deriv(0,node)*xji(dim,0) + deriv(1,node)*xji(dim,1) + deriv(2,node)*xji(dim,2))
                        *fac;
      }
    }
  }
}

//////
// Calculate length of edge and differences in each dimension for two nodes.
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_edge_geometry(
  int i, int j,
  const LINALG::Matrix<3, iel>& xyze,
  double& length,
  double& dx,
  double& dy,
  double& dz)
{
  /*---------------------------------------------- x-, y- and z-difference ---*/
  dx = xyze(0,j)-xyze(0,i);
  dy = xyze(1,j)-xyze(1,i);
  dz = xyze(2,j)-xyze(2,i);
  /*------------------------------- determine distance between i and j ---*/
  length = sqrt(dx * dx + dy * dy + dz * dz);
#ifdef DEBUG
  if (length < (1.0E-14)) dserror("edge or diagonal of element has zero length");
#endif
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_add_tria_stiffness(
    int node_p, int node_q, int node_r, int node_s,
    const LINALG::Matrix<3, 1>& sq,
    const double len_sq,
    const LINALG::Matrix<3, 1>& rp,
    const double len_rp,
    const LINALG::Matrix<3, 1>& qp,
    const LINALG::Matrix<3, 1>& local_x,
    LINALG::Matrix<3*iel,3*iel>& sys_mat)
{
  //Positions for dynamic triangle (2D)
  //sequence: s,j,q
  LINALG::Matrix<2,3> xyze_dyn_tria;

  // Some matrices that can be found in the paper are not assembled
  // here, because it can be done with less memory. I use only one
  // transformation matrix to get from the triangle plane to the 12x12
  // stiffness component of this triangle.

  //transformation matrix from the plane of the triangle to the
  //three-dimensional global frame.
  //This corresponds to (R(sjq) x r(sjq) x S^T)^T from Farhat et al.
  LINALG::Matrix<12, 3> trans_matrix;

  //rotational stiffness matrix for tetrahedron with given dynamic triangle
  LINALG::Matrix<12,12> k_dyn_tet;

  //local x,y in the plane of the dynamic triangle
  // these are the 3d-vectors that span this plane
  LINALG::Matrix<3,1> local_y;

  // the point p relativ to the point (0,0) in the triangle plane
  // transformed into 3d space
  LINALG::Matrix<3,1> p;

  //local x-value of s xyze_dyn_tria(0,0) := (-1.0)*sq*local_x (local origin lies on plane pqr)
  //local y-value of s xyze_dyn_tria(1,0 := 0.0
  xyze_dyn_tria(0,0) = - sq.Dot(local_x);
  xyze_dyn_tria(1,0) = 0.0;

  //local_y = (sq + xyze_dyn_tria(0,0)*local_x)/|(sq + xyze_dyn_tria(0,0)*local_x)|
  //xyze_dyn_tria(1,2) = |sq + xyze_dyn_tria(0,0)*local_x|, xyze_dyn_tria(0,2) := 0.0
  local_y.Update(xyze_dyn_tria(0,0), local_x, 1.0, sq);  //just an intermediate step
  xyze_dyn_tria(1,2) = local_y.Norm2();
  xyze_dyn_tria(0,2) = 0.0;

  if (xyze_dyn_tria(1,2) != 0)  // == 0 will trigger the parallel check below
    local_y.Scale(1.0/xyze_dyn_tria(1,2));

  // check = (local_y x sq) * rp  If this is very small rp is parallel
  // to the plane spanned by local_y and sq.
  double check = (local_y(1)*sq(2) - local_y(2)*sq(1)) * rp(0)
                 + (local_y(2)*sq(0) - local_y(0)*sq(2)) * rp(1)
                 + (local_y(0)*sq(1) - local_y(1)*sq(0)) * rp(2);
  check /= sqrt((rp(0)*rp(0) + rp(1)*rp(1) + rp(2)*rp(2)) * (sq(0)*sq(0) + sq(1)*sq(1) + sq(2)*sq(2)));


  // if rp and local_y are parallel calculate stiffness of lineal spring s-q
  if (fabs(check) < 1e-2 or fabs(xyze_dyn_tria(1,2)) < 1e-9 )
  {
    const int ts   = 3*node_s;
    const int tsp  = ts+1;
    const int tspp = ts+2;
    const int tq   = 3*node_q;
    const int tqp  = tq+1;
    const int tqpp = tq+2;
    // we know the edge-information from above.
    const double factor = 1.0 / (len_sq*len_sq*len_sq);
    const double dxx_l3 = sq(0)*sq(0)*factor;
    const double dxy_l3 = sq(0)*sq(1)*factor;
    const double dxz_l3 = sq(0)*sq(2)*factor;
    const double dyy_l3 = sq(1)*sq(1)*factor;
    const double dyz_l3 = sq(1)*sq(2)*factor;
    const double dzz_l3 = sq(2)*sq(2)*factor;
    //put values in 'element stiffness'
    //rows for node_s
    sys_mat(ts,    ts  ) += dxx_l3;
    sys_mat(ts,    tsp ) += dxy_l3;
    sys_mat(ts,    tspp) += dxz_l3;

    sys_mat(ts,    tq  ) -= dxx_l3;
    sys_mat(ts,    tqp ) -= dxy_l3;
    sys_mat(ts,    tqpp) -= dxz_l3;

    sys_mat(tsp,   ts  ) += dxy_l3;
    sys_mat(tsp,   tsp ) += dyy_l3;
    sys_mat(tsp,   tspp) += dyz_l3;

    sys_mat(tsp,   tq  ) -= dxy_l3;
    sys_mat(tsp,   tqp ) -= dyy_l3;
    sys_mat(tsp,   tqpp) -= dyz_l3;

    sys_mat(tspp,  ts  ) += dxz_l3;
    sys_mat(tspp,  tsp ) += dyz_l3;
    sys_mat(tspp,  tspp) += dzz_l3;

    sys_mat(tspp,  tq  ) -= dxz_l3;
    sys_mat(tspp,  tqp ) -= dyz_l3;
    sys_mat(tspp,  tqpp) -= dzz_l3;

    //rows for node_q
    sys_mat(tq,    ts  ) -= dxx_l3;
    sys_mat(tq,    tsp ) -= dxy_l3;
    sys_mat(tq,    tspp) -= dxz_l3;

    sys_mat(tq,    tq  ) += dxx_l3;
    sys_mat(tq,    tqp ) += dxy_l3;
    sys_mat(tq,    tqpp) += dxz_l3;

    sys_mat(tqp,   ts  ) -= dxy_l3;
    sys_mat(tqp,   tsp ) -= dyy_l3;
    sys_mat(tqp,   tspp) -= dyz_l3;

    sys_mat(tqp,   tq  ) += dxy_l3;
    sys_mat(tqp,   tqp ) += dyy_l3;
    sys_mat(tqp,   tqpp) += dyz_l3;

    sys_mat(tqpp,  ts  ) -= dxz_l3;
    sys_mat(tqpp,  tsp ) -= dyz_l3;
    sys_mat(tqpp,  tspp) -= dzz_l3;

    sys_mat(tqpp,  tq  ) += dxz_l3;
    sys_mat(tqpp,  tqp ) += dyz_l3;
    sys_mat(tqpp,  tqpp) += dzz_l3;
  }
  else
  {
    //local x,y-values of j, using pO + Oj + jp = 0
    //(O is local origin on plane pqr)
    xyze_dyn_tria(0,1) = 0.0;
    p.Update(xyze_dyn_tria(1,2), local_y, 1, qp);
#if 0
    {
      FILE* errfile = DRT::Problem::Instance()->ErrorFile()->Handle();

      fprintf(errfile,"\n");
      fprintf(errfile,"rp=matrix([% e,% e,% e]).transpose()\n",rp(0),rp(1),rp(2));
      fprintf(errfile,"pq=matrix([% e,% e,% e]).transpose()\n",pq(0),pq(1),pq(2));
      fprintf(errfile,"sq=matrix([% e,% e,% e]).transpose()\n",sq(0),sq(1),sq(2));
      fprintf(errfile,"p=matrix([% e,% e,% e]).transpose()\n",p(0),p(1),p(2));
      fflush(errfile);
    }
#endif

    ///// solve xyze_dyn_tria(1,1) * local_y = p + lambda * (-rp)
    double d;
    double lambda;
    const double limit = 1e-4 * len_rp; // I think the limit can be even higher
    if (fabs(d = local_y(0)*rp(1) - local_y(1)*rp(0)) > limit) {
      const double fac = 1.0/d;
      xyze_dyn_tria(1,1) = (rp(1)*p(0) - rp(0)*p(1))*fac;
      lambda = (local_y(0)*p(1) - local_y(1)*p(0))*fac;
    } else if (fabs(d = local_y(0)*rp(2) - local_y(2)*rp(0)) > limit) {
      const double fac = 1.0/d;
      xyze_dyn_tria(1,1) = (rp(2)*p(0) - rp(0)*p(2))*fac;
      lambda = (local_y(0)*p(2) - local_y(2)*p(0))*fac;
    } else /* we know it has to work here, because the system is solvable  */ {
      const double fac = 1.0/(local_y(1)*rp(2) - local_y(2)*rp(1));
      xyze_dyn_tria(1,1) = (rp(2)*p(1) - rp(1)*p(2))*fac;
      lambda = (local_y(1)*p(2) - local_y(2)*p(1))*fac;
    }

    if (lambda < 0.0)
      lambda = 0.0;
    else if (lambda > 1.0)
      lambda = 1.0;
    const double one_minus_lambda = 1-lambda;

    ////// evaluate torsional stiffness of dynamic triangle
    const double& tmp = xyze_dyn_tria(0,0);
    const double y_jk = -xyze_dyn_tria(1,1) + xyze_dyn_tria(1,2);
    // squares of side lengths
    const double l_ij_sq = xyze_dyn_tria(1,1)*xyze_dyn_tria(1,1) + tmp*tmp;
    const double l_jk_sq = y_jk*y_jk;
    const double l_ki_sq = xyze_dyn_tria(1,2)*xyze_dyn_tria(1,2) + tmp*tmp;
    // auxiliary values same as in Farhat et al.
    const double a_ij = -tmp / (l_ij_sq);
    const double a_jk = 0.0;  // 0.0 / (l_jk_sq)
    const double a_ki = tmp / (l_ki_sq);
    const double b_ij = xyze_dyn_tria(1,1) / (l_ij_sq);
    const double b_jk = y_jk / (l_jk_sq);
    const double b_ki = -xyze_dyn_tria(1,2) / (l_ki_sq);

    const double a_ij_0 = a_ij * local_y(0);
    const double a_ij_1 = a_ij * local_y(1);
    const double a_ij_2 = a_ij * local_y(2);
    const double a_jk_0 = a_jk * local_y(0);
    const double a_jk_1 = a_jk * local_y(1);
    const double a_jk_2 = a_jk * local_y(2);
    const double a_ki_0 = a_ki * local_y(0);
    const double a_ki_1 = a_ki * local_y(1);
    const double a_ki_2 = a_ki * local_y(2);
    const double b_ij_0 = b_ij * local_x(0);
    const double b_ij_1 = b_ij * local_x(1);
    const double b_ij_2 = b_ij * local_x(2);
    const double b_jk_0 = b_jk * local_x(0);
    const double b_jk_1 = b_jk * local_x(1);
    const double b_jk_2 = b_jk * local_x(2);
    const double b_ki_0 = b_ki * local_x(0);
    const double b_ki_1 = b_ki * local_x(1);
    const double b_ki_2 = b_ki * local_x(2);

    // area of the triangle
    const double area_double = 0.5 * sqrt(2.0*l_ij_sq*l_jk_sq + 2.0*l_jk_sq*l_ki_sq + 2.0*l_ki_sq*l_ij_sq
                                          - l_ij_sq*l_ij_sq - l_jk_sq*l_jk_sq - l_ki_sq*l_ki_sq);
    const double area_double_sqare = area_double * area_double;


#ifdef DEBUG            /*---------------------------------- check edge lengths ---*/
    if (l_ij_sq < (1.0E-7)) dserror("edge or diagonal of element has zero length");
    if (l_jk_sq < (1.0E-7)) dserror("edge or diagonal of element has zero length");
    if (l_ki_sq < (1.0E-7)) dserror("edge or diagonal of element has zero length");
#endif

/*---------------------------------- determine torsional stiffnesses ---*/
    // instead of the whole matrix only the non-zero diagonal is stored
    const double C0 = l_ij_sq * l_ki_sq / area_double_sqare;
    const double C1 = l_ij_sq * l_jk_sq / area_double_sqare;
    const double C2 = l_ki_sq * l_jk_sq / area_double_sqare;

/*--------------------------------------- fill transformation matrix ---*/
    // This corresponds to (R(sjq) x r(sjq) x S^T)^T from Farhat et al.
    trans_matrix(0,0)  = one_minus_lambda*(b_ij_0 - a_ij_0);
    trans_matrix(1,0)  = one_minus_lambda*(b_ij_1 - a_ij_1);
    trans_matrix(2,0)  = one_minus_lambda*(b_ij_2 - a_ij_2);
    trans_matrix(3,0)  = b_ki_0 - a_ki_0;
    trans_matrix(4,0)  = b_ki_1 - a_ki_1;
    trans_matrix(5,0)  = b_ki_2 - a_ki_2;
    trans_matrix(6,0)  = lambda*(b_ij_0 - a_ij_0);
    trans_matrix(7,0)  = lambda*(b_ij_1 - a_ij_1);
    trans_matrix(8,0)  = lambda*(b_ij_2 - a_ij_2);
    trans_matrix(9,0)  = a_ij_0 + a_ki_0 - b_ki_0 - b_ij_0;
    trans_matrix(10,0) = a_ij_1 + a_ki_1 - b_ki_1 - b_ij_1;
    trans_matrix(11,0) = a_ij_2 + a_ki_2 - b_ki_2 - b_ij_2;

    trans_matrix(0,1)  = one_minus_lambda*(a_jk_0 + a_ij_0 - b_ij_0 - b_jk_0);
    trans_matrix(1,1)  = one_minus_lambda*(a_jk_1 + a_ij_1 - b_ij_1 - b_jk_1);
    trans_matrix(2,1)  = one_minus_lambda*(a_jk_2 + a_ij_2 - b_ij_2 - b_jk_2);
    trans_matrix(3,1)  = b_jk_0 - a_jk_0;
    trans_matrix(4,1)  = b_jk_1 - a_jk_1;
    trans_matrix(5,1)  = b_jk_2 - a_jk_2;
    trans_matrix(6,1)  = lambda*(a_jk_0 + a_ij_0 - b_ij_0 - b_jk_0);
    trans_matrix(7,1)  = lambda*(a_jk_1 + a_ij_1 - b_ij_1 - b_jk_1);
    trans_matrix(8,1)  = lambda*(a_jk_2 + a_ij_2 - b_ij_2 - b_jk_2);
    trans_matrix(9,1)  = b_ij_0 - a_ij_0;
    trans_matrix(10,1) = b_ij_1 - a_ij_1;
    trans_matrix(11,1) = b_ij_2 - a_ij_2;

    trans_matrix(0,2)  = one_minus_lambda*(b_jk_0 - a_jk_0);
    trans_matrix(1,2)  = one_minus_lambda*(b_jk_1 - a_jk_1);
    trans_matrix(2,2)  = one_minus_lambda*(b_jk_2 - a_jk_2);
    trans_matrix(3,2)  = a_ki_0 + a_jk_0 - b_jk_0 - b_ki_0;
    trans_matrix(4,2)  = a_ki_1 + a_jk_1 - b_jk_1 - b_ki_1;
    trans_matrix(5,2)  = a_ki_2 + a_jk_2 - b_jk_2 - b_ki_2;
    trans_matrix(6,2)  = lambda*(b_jk_0 - a_jk_0);
    trans_matrix(7,2)  = lambda*(b_jk_1 - a_jk_1);
    trans_matrix(8,2)  = lambda*(b_jk_2 - a_jk_2);
    trans_matrix(9,2)  = b_ki_0 - a_ki_0;
    trans_matrix(10,2) = b_ki_1 - a_ki_1;
    trans_matrix(11,2) = b_ki_2 - a_ki_2;

/*----------------------------------- perform matrix multiplications ---*/
    ////// Compute k_dyn_tet
    // trans_matrix x C x trans_matrix^(T)
    // This corresponds to S x r(sjq)^T x R(sjq)^T x C x R(sjq) x r(sjq) x S^T
    for (int i = 0; i < 12; ++i) {
        const double C0T0 = C0 * trans_matrix(i, 0);
        const double C1T1 = C1 * trans_matrix(i, 1);
        const double C2T2 = C2 * trans_matrix(i, 2);
        k_dyn_tet(i, i) = C0T0*trans_matrix(i, 0) +
                          C1T1*trans_matrix(i, 1) +
                          C2T2*trans_matrix(i, 2);
        for (int j = 0; j < i; ++j) {
            k_dyn_tet(i, j) = k_dyn_tet(j, i) = C0T0*trans_matrix(j, 0) +
                                                C1T1*trans_matrix(j, 1) +
                                                C2T2*trans_matrix(j, 2);
        }
    }

    // This makes writing the loop over sys_mat much easier
    std::vector<int> matrix_access(4);
    matrix_access[0] = 3*node_p;
    matrix_access[1] = 3*node_q;
    matrix_access[2] = 3*node_r;
    matrix_access[3] = 3*node_s;

    //Sort values in element's sys_mat
    for (int i=0; i<4; ++i) {
      const int ti   = 3*i;
      const int tip  = ti+1;
      const int tipp = ti+2;
      const int sys_mat_i = matrix_access[i];
      for (int j=0; j<4; ++j) {
        const int tj   = 3*j;
        const int tjp  = tj+1;
        const int tjpp = tj+2;
        const int sys_mat_j = matrix_access[j];

        sys_mat(sys_mat_i  , sys_mat_j  ) += k_dyn_tet(ti  , tj  );
        sys_mat(sys_mat_i  , sys_mat_j+1) += k_dyn_tet(ti  , tjp );
        sys_mat(sys_mat_i  , sys_mat_j+2) += k_dyn_tet(ti  , tjpp);

        sys_mat(sys_mat_i+1, sys_mat_j  ) += k_dyn_tet(tip , tj  );
        sys_mat(sys_mat_i+1, sys_mat_j+1) += k_dyn_tet(tip , tjp );
        sys_mat(sys_mat_i+1, sys_mat_j+2) += k_dyn_tet(tip , tjpp);

        sys_mat(sys_mat_i+2, sys_mat_j  ) += k_dyn_tet(tipp, tj  );
        sys_mat(sys_mat_i+2, sys_mat_j+1) += k_dyn_tet(tipp, tjp );
        sys_mat(sys_mat_i+2, sys_mat_j+2) += k_dyn_tet(tipp, tjpp);
      }
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_add_tetra_stiffness(
  int tet_0, int tet_1, int tet_2, int tet_3,
  LINALG::Matrix<3*iel,3*iel>& sys_mat,
  const LINALG::Matrix<3,iel>& xyze)
{
  //according to Farhat et al.
  //twelve-triangle configuration

  LINALG::Matrix<3, 1> e01, e02, e03, e10, e12, e13, e20, e21, e23, e30, e31, e32, local_x;
  double l01, l02, l03, l12, l13, l23;
  ale3_edge_geometry(tet_0, tet_1, xyze, l01, e01(0), e01(1), e01(2));
  ale3_edge_geometry(tet_0, tet_2, xyze, l02, e02(0), e02(1), e02(2));
  ale3_edge_geometry(tet_0, tet_3, xyze, l03, e03(0), e03(1), e03(2));
  ale3_edge_geometry(tet_1, tet_2, xyze, l12, e12(0), e12(1), e12(2));
  ale3_edge_geometry(tet_1, tet_3, xyze, l13, e13(0), e13(1), e13(2));
  ale3_edge_geometry(tet_2, tet_3, xyze, l23, e23(0), e23(1), e23(2));
  e10(0) = -e01(0);  e10(1) = -e01(1);  e10(2) = -e01(2);
  e20(0) = -e02(0);  e20(1) = -e02(1);  e20(2) = -e02(2);
  e21(0) = -e12(0);  e21(1) = -e12(1);  e21(2) = -e12(2);
  e30(0) = -e03(0);  e30(1) = -e03(1);  e30(2) = -e03(2);
  e31(0) = -e13(0);  e31(1) = -e13(1);  e31(2) = -e13(2);
  e32(0) = -e23(0);  e32(1) = -e23(1);  e32(2) = -e23(2);

  local_x(0) = e12(1)*e23(2)-e12(2)*e23(1);
  local_x(1) = e12(2)*e23(0)-e12(0)*e23(2);
  local_x(2) = e12(0)*e23(1)-e12(1)*e23(0);
  local_x.Scale(1.0/local_x.Norm2());
  ale3_add_tria_stiffness(tet_1, tet_2, tet_3, tet_0, e02, l02, e31, l13, e21, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_1, tet_3, tet_2, tet_0, e03, l03, e21, l12, e31, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_2, tet_1, tet_3, tet_0, e01, l01, e32, l23, e12, local_x, sys_mat);

  local_x(0) = e23(1)*e03(2)-e23(2)*e03(1);
  local_x(1) = e23(2)*e03(0)-e23(0)*e03(2);
  local_x(2) = e23(0)*e03(1)-e23(1)*e03(0);
  local_x.Scale(1.0/local_x.Norm2());
  ale3_add_tria_stiffness(tet_0, tet_2, tet_3, tet_1, e12, l12, e30, l03, e20, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_0, tet_3, tet_2, tet_1, e13, l13, e20, l02, e30, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_2, tet_0, tet_3, tet_1, e10, l01, e32, l23, e02, local_x, sys_mat);

  local_x(0) = e01(1)*e03(2)-e01(2)*e03(1);
  local_x(1) = e01(2)*e03(0)-e01(0)*e03(2);
  local_x(2) = e01(0)*e03(1)-e01(1)*e03(0);
  local_x.Scale(1.0/local_x.Norm2());
  ale3_add_tria_stiffness(tet_0, tet_1, tet_3, tet_2, e21, l12, e30, l03, e10, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_1, tet_0, tet_3, tet_2, e20, l02, e31, l13, e01, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_0, tet_3, tet_1, tet_2, e23, l23, e10, l01, e30, local_x, sys_mat);

  local_x(0) = e01(1)*e12(2)-e01(2)*e12(1);
  local_x(1) = e01(2)*e12(0)-e01(0)*e12(2);
  local_x(2) = e01(0)*e12(1)-e01(1)*e12(0);
  local_x.Scale(1.0/local_x.Norm2());
  ale3_add_tria_stiffness(tet_0, tet_1, tet_2, tet_3, e31, l13, e20, l02, e10, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_0, tet_2, tet_1, tet_3, e32, l23, e10, l01, e20, local_x, sys_mat);
  ale3_add_tria_stiffness(tet_1, tet_0, tet_2, tet_3, e30, l03, e21, l12, e01, local_x, sys_mat);
}

//////
// dummy function, "divide tetra into tetras"
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_tors_spring_tet4(
  LINALG::Matrix<3*iel,3*iel>& sys_mat,
  const LINALG::Matrix<3,iel>& xyze)
{
  ale3_add_tetra_stiffness(0,1,2,3,sys_mat,xyze);
}

//////
// divide pyramid into tetras
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_tors_spring_pyramid5(
  LINALG::Matrix<3*iel,3*iel>& sys_mat,
  const LINALG::Matrix<3,iel>& xyze)
{
  ale3_add_tetra_stiffness(0,1,3,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(0,1,2,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(1,2,3,4,sys_mat,xyze);
  ale3_add_tetra_stiffness(0,2,3,4,sys_mat,xyze);
}

//////
// divide wedge into tetras
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_tors_spring_wedge6(
  LINALG::Matrix<3*iel,3*iel>& sys_mat,
  const LINALG::Matrix<3,iel>& xyze)
{
  ale3_add_tetra_stiffness(2, 0, 1, 3, sys_mat, xyze);
  ale3_add_tetra_stiffness(2, 0, 1, 4, sys_mat, xyze);
  ale3_add_tetra_stiffness(2, 0, 1, 5, sys_mat, xyze);

  ale3_add_tetra_stiffness(3, 0, 1, 5, sys_mat, xyze);
  ale3_add_tetra_stiffness(3, 0, 2, 4, sys_mat, xyze);
  ale3_add_tetra_stiffness(3, 1, 2, 4, sys_mat, xyze);
  ale3_add_tetra_stiffness(3, 1, 2, 5, sys_mat, xyze);

  ale3_add_tetra_stiffness(4, 0, 1, 5, sys_mat, xyze);
  ale3_add_tetra_stiffness(4, 0, 2, 5, sys_mat, xyze);
  ale3_add_tetra_stiffness(4, 0, 3, 5, sys_mat, xyze);
  ale3_add_tetra_stiffness(4, 1, 3, 5, sys_mat, xyze);
  ale3_add_tetra_stiffness(4, 2, 3, 5, sys_mat, xyze);
}

//////
// divide hex into tetras
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_tors_spring_hex8(
  LINALG::Matrix<3*iel,3*iel>& sys_mat,
  const LINALG::Matrix<3,iel>& xyze)
{

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

//////
// divide nurbs into tetras
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::Ale3_Impl<distype>::ale3_tors_spring_nurbs27(
  LINALG::Matrix<3*iel,3*iel>& sys_mat,
  const LINALG::Matrix<3,iel>& xyze)
{

  //                          v
  //                         /
  //                        /
  //          X---------X---------X
  //         /|24      /|25      /|26
  //  w     / |       / |       / |
  //  ^    /  |      /  |      /  |
  //  |   X---------X---------X   |
  //  |  /|21 |    /|22 |    /|23 |
  //  | / |   X---/-|---X---/-|---X
  //   /  |  /|15/  |  /|16/  |  /|17                X---------X
  //  X---------X---------X   | / |                 /|7       /|6     0,1,3,4
  //  |18 |/  | |19 |/  | |20 |/  |                / |       / |      0,1,2,5
  //  |   X-----|---X-----|---X   |               /  |      /  |      1,2,3,6
  //  |  /|12 | |  /|13 | |  /|14 |              X---------X   |      0,2,3,7
  //  | / |   X-|-/-|---X-|-/-|---X              |4  |     |5  |
  //  |/  |  /6 |/  |  /7 |/  |  /8              |   X-----|---X      0,4,5,7
  //  X---------X---------X   | /                |  /3     |  /2      1,4,5,6
  //  |9  |/    |10 |/    |11 |/    	         | /       | /        2,5,6,7
  //  |   X-----|---X-----|---X     	         |/        |/         3,4,6,7
  //  |  /3     |  /4     |  /5		         X---------X
  //  | /       | /       | /		          0         1
  //  |/        |/        |/
  //  X---------X---------X ----->u
  //   0         1         2
  //

  //Use 8 tetrahedra to prevent node-face-penetration

  ale3_add_tetra_stiffness(0, 1, 3, 9,sys_mat,xyze);
  ale3_add_tetra_stiffness(0, 1, 4,10,sys_mat,xyze);
  ale3_add_tetra_stiffness(1, 4, 3,13,sys_mat,xyze);
  ale3_add_tetra_stiffness(0, 4, 3,12,sys_mat,xyze);

  ale3_add_tetra_stiffness(0, 9,10,12,sys_mat,xyze);
  ale3_add_tetra_stiffness(1, 9,10,13,sys_mat,xyze);
  ale3_add_tetra_stiffness(4,10,13,12,sys_mat,xyze);
  ale3_add_tetra_stiffness(3, 9,13,12,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness(1, 2, 4,10,sys_mat,xyze);
  ale3_add_tetra_stiffness(1, 2, 5,11,sys_mat,xyze);
  ale3_add_tetra_stiffness(2, 5, 4,14,sys_mat,xyze);
  ale3_add_tetra_stiffness(1, 5, 4,13,sys_mat,xyze);

  ale3_add_tetra_stiffness(1,10,11,13,sys_mat,xyze);
  ale3_add_tetra_stiffness(2,10,11,14,sys_mat,xyze);
  ale3_add_tetra_stiffness(5,11,14,13,sys_mat,xyze);
  ale3_add_tetra_stiffness(4,10,14,13,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness(4, 5, 7,13,sys_mat,xyze);
  ale3_add_tetra_stiffness(4, 5, 8,14,sys_mat,xyze);
  ale3_add_tetra_stiffness(5, 8, 7,17,sys_mat,xyze);
  ale3_add_tetra_stiffness(4, 8, 7,16,sys_mat,xyze);

  ale3_add_tetra_stiffness(4,13,14,16,sys_mat,xyze);
  ale3_add_tetra_stiffness(5,13,14,17,sys_mat,xyze);
  ale3_add_tetra_stiffness(8,14,17,16,sys_mat,xyze);
  ale3_add_tetra_stiffness(7,13,17,16,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness(3, 4, 6,12,sys_mat,xyze);
  ale3_add_tetra_stiffness(3, 4, 7,13,sys_mat,xyze);
  ale3_add_tetra_stiffness(4, 7, 6,16,sys_mat,xyze);
  ale3_add_tetra_stiffness(3, 7, 6,15,sys_mat,xyze);

  ale3_add_tetra_stiffness(3,12,13,15,sys_mat,xyze);
  ale3_add_tetra_stiffness(4,12,13,16,sys_mat,xyze);
  ale3_add_tetra_stiffness(7,13,16,15,sys_mat,xyze);
  ale3_add_tetra_stiffness(6,12,16,15,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness( 9,10,12,18,sys_mat,xyze);
  ale3_add_tetra_stiffness( 9,10,13,19,sys_mat,xyze);
  ale3_add_tetra_stiffness(10,13,12,22,sys_mat,xyze);
  ale3_add_tetra_stiffness( 9,13,12,21,sys_mat,xyze);

  ale3_add_tetra_stiffness( 9,18,19,21,sys_mat,xyze);
  ale3_add_tetra_stiffness(10,18,19,22,sys_mat,xyze);
  ale3_add_tetra_stiffness(13,19,22,21,sys_mat,xyze);
  ale3_add_tetra_stiffness(12,18,22,21,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness(10,11,13,19,sys_mat,xyze);
  ale3_add_tetra_stiffness(10,11,14,20,sys_mat,xyze);
  ale3_add_tetra_stiffness(11,14,13,23,sys_mat,xyze);
  ale3_add_tetra_stiffness(10,14,13,22,sys_mat,xyze);

  ale3_add_tetra_stiffness(10,19,20,22,sys_mat,xyze);
  ale3_add_tetra_stiffness(11,19,20,23,sys_mat,xyze);
  ale3_add_tetra_stiffness(14,20,23,22,sys_mat,xyze);
  ale3_add_tetra_stiffness(13,19,23,22,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness(13,14,16,22,sys_mat,xyze);
  ale3_add_tetra_stiffness(13,14,17,23,sys_mat,xyze);
  ale3_add_tetra_stiffness(14,17,16,26,sys_mat,xyze);
  ale3_add_tetra_stiffness(13,17,16,25,sys_mat,xyze);

  ale3_add_tetra_stiffness(13,22,23,25,sys_mat,xyze);
  ale3_add_tetra_stiffness(14,22,23,26,sys_mat,xyze);
  ale3_add_tetra_stiffness(17,23,26,25,sys_mat,xyze);
  ale3_add_tetra_stiffness(16,22,26,25,sys_mat,xyze);

  // -----------------------------------------------

  ale3_add_tetra_stiffness(12,13,15,21,sys_mat,xyze);
  ale3_add_tetra_stiffness(12,13,16,22,sys_mat,xyze);
  ale3_add_tetra_stiffness(13,16,15,25,sys_mat,xyze);
  ale3_add_tetra_stiffness(12,16,15,24,sys_mat,xyze);

  ale3_add_tetra_stiffness(12,21,22,24,sys_mat,xyze);
  ale3_add_tetra_stiffness(13,21,22,25,sys_mat,xyze);
  ale3_add_tetra_stiffness(16,22,25,24,sys_mat,xyze);
  ale3_add_tetra_stiffness(15,21,25,24,sys_mat,xyze);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Ale3_Impl<distype>::static_ke_spring(
  Ale3*                     ele,
  Epetra_SerialDenseMatrix& sys_mat_epetra,
  const vector<double>&           displacements)
{
  LINALG::Matrix<3*iel,3*iel> sys_mat(sys_mat_epetra.A(),true);
  int node_i, node_j;                                     // end nodes of spring
  double length;                                          // length of edge
  double dx, dy, dz;                                      // deltas in each direction
  double factor;

  // get node coordinates
  LINALG::Matrix<3,iel> xyze;
  DRT::Node** nodes = ele->Nodes();
  for(int i=0;i<iel;i++)
  {
    const double* x = nodes[i]->X();
    xyze(0,i) = x[0] + displacements[i*3];
    xyze(1,i) = x[1] + displacements[i*3+1];
    xyze(2,i) = x[2] + displacements[i*3+2];
  }

//lineal springs from all corner nodes to all corner nodes
//loop over all edges and diagonals of the element
  // We depend on the corner nodes to be the first nodes in xyze!
  for (node_i=0; node_i<(numcnd-1); ++node_i)
  {
    const int ti   = 3*node_i;
    const int tip  = ti+1;
    const int tipp = ti+2;

    double ii      = 0.0;
    double iip     = 0.0;
    double iipp    = 0.0;
    double ipi     = 0.0;
    double ipip    = 0.0;
    double ipipp   = 0.0;
    double ippi    = 0.0;
    double ippip   = 0.0;
    double ippipp  = 0.0;
    for (node_j=node_i+1; node_j<numcnd; ++node_j)
    {
      const int tj   = 3*node_j;
      const int tjp  = tj+1;
      const int tjpp = tj+2;
      ale3_edge_geometry(node_i,node_j,xyze,length,dx,dy,dz);
      factor = 1.0/(length*length*length);
      const double dxx_l3 = dx*dx*factor;
      const double dxy_l3 = dx*dy*factor;
      const double dxz_l3 = dx*dz*factor;
      const double dyy_l3 = dy*dy*factor;
      const double dyz_l3 = dy*dz*factor;
      const double dzz_l3 = dz*dz*factor;

      //put values in 'element stiffness'
      //rows for node_i
      ii     += dxx_l3;
      iip    += dxy_l3;
      iipp   += dxz_l3;

      sys_mat(ti,    tj  ) -= dxx_l3;
      sys_mat(ti,    tjp ) -= dxy_l3;
      sys_mat(ti,    tjpp) -= dxz_l3;

      ipi    += dxy_l3;
      ipip   += dyy_l3;
      ipipp  += dyz_l3;

      sys_mat(tip,   tj  ) -= dxy_l3;
      sys_mat(tip,   tjp ) -= dyy_l3;
      sys_mat(tip,   tjpp) -= dyz_l3;

      ippi   += dxz_l3;
      ippip  += dyz_l3;
      ippipp += dzz_l3;

      sys_mat(tipp,  tj  ) -= dxz_l3;
      sys_mat(tipp,  tjp ) -= dyz_l3;
      sys_mat(tipp,  tjpp) -= dzz_l3;

      //rows for node_j
      sys_mat(tj,    ti  ) -= dxx_l3;
      sys_mat(tj,    tip ) -= dxy_l3;
      sys_mat(tj,    tipp) -= dxz_l3;

      sys_mat(tj,    tj  ) += dxx_l3;
      sys_mat(tj,    tjp ) += dxy_l3;
      sys_mat(tj,    tjpp) += dxz_l3;

      sys_mat(tjp,   ti  ) -= dxy_l3;
      sys_mat(tjp,   tip ) -= dyy_l3;
      sys_mat(tjp,   tipp) -= dyz_l3;

      sys_mat(tjp,   tj  ) += dxy_l3;
      sys_mat(tjp,   tjp ) += dyy_l3;
      sys_mat(tjp,   tjpp) += dyz_l3;

      sys_mat(tjpp,  ti  ) -= dxz_l3;
      sys_mat(tjpp,  tip ) -= dyz_l3;
      sys_mat(tjpp,  tipp) -= dzz_l3;

      sys_mat(tjpp,  tj  ) += dxz_l3;
      sys_mat(tjpp,  tjp ) += dyz_l3;
      sys_mat(tjpp,  tjpp) += dzz_l3;
    }
    sys_mat(ti  ,ti  ) += ii;
    sys_mat(ti  ,tip ) += iip;
    sys_mat(ti  ,tipp) += iipp;
    sys_mat(tip ,ti  ) += ipi;
    sys_mat(tip ,tip ) += ipip;
    sys_mat(tip ,tipp) += ipipp;
    sys_mat(tipp,ti  ) += ippi;
    sys_mat(tipp,tip ) += ippip;
    sys_mat(tipp,tipp) += ippipp;
  }

  //build in torsional springs
  //and put edge nodes on the middle of the respective edge
  switch (distype)
  {
  case DRT::Element::tet10:

    sys_mat(12,0)  = -0.5;
    sys_mat(12,3)  = -0.5;
    sys_mat(12,12) =  1.0;
    sys_mat(13,1)  = -0.5;
    sys_mat(13,4)  = -0.5;
    sys_mat(13,13) =  1.0;
    sys_mat(14,2)  = -0.5;
    sys_mat(14,5)  = -0.5;
    sys_mat(14,14) =  1.0;

    sys_mat(15,3)  = -0.5;
    sys_mat(15,6)  = -0.5;
    sys_mat(15,15) =  1.0;
    sys_mat(16,4)  = -0.5;
    sys_mat(16,7)  = -0.5;
    sys_mat(16,16) =  1.0;
    sys_mat(17,5)  = -0.5;
    sys_mat(17,8)  = -0.5;
    sys_mat(17,17) =  1.0;

    sys_mat(18,6)  = -0.5;
    sys_mat(18,0)  = -0.5;
    sys_mat(18,18) =  1.0;
    sys_mat(19,7)  = -0.5;
    sys_mat(19,1)  = -0.5;
    sys_mat(19,19) =  1.0;
    sys_mat(20,8)  = -0.5;
    sys_mat(20,2)  = -0.5;
    sys_mat(20,20) =  1.0;

    sys_mat(21,0)  = -0.5;
    sys_mat(21,9)  = -0.5;
    sys_mat(21,21) =  1.0;
    sys_mat(22,1)  = -0.5;
    sys_mat(22,10) = -0.5;
    sys_mat(22,22) =  1.0;
    sys_mat(23,2)  = -0.5;
    sys_mat(23,11) = -0.5;
    sys_mat(23,23) =  1.0;

    sys_mat(24,3)  = -0.5;
    sys_mat(24,9)  = -0.5;
    sys_mat(24,24) =  1.0;
    sys_mat(25,4)  = -0.5;
    sys_mat(25,10) = -0.5;
    sys_mat(25,25) =  1.0;
    sys_mat(26,5)  = -0.5;
    sys_mat(26,11) = -0.5;
    sys_mat(26,26) =  1.0;

    sys_mat(27,6)  = -0.5;
    sys_mat(27,9)  = -0.5;
    sys_mat(27,27) =  1.0;
    sys_mat(28,7)  = -0.5;
    sys_mat(28,10) = -0.5;
    sys_mat(28,28) =  1.0;
    sys_mat(29,8)  = -0.5;
    sys_mat(29,11) = -0.5;
    sys_mat(29,29) =  1.0;

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
      const int tk = 3*k;
      const int t6k = tk + 3*6;
      const int a = 3 * ((k < 2) ? (k+1) : 0);
      sys_mat(t6k   , tk    ) = -0.5;
      sys_mat(t6k +1, tk +1 ) = -0.5;
      sys_mat(t6k +2, tk +2 ) = -0.5;
      sys_mat(t6k   , a    ) = -0.5;
      sys_mat(t6k +1, a +1 ) = -0.5;
      sys_mat(t6k +2, a +2 ) = -0.5;
      sys_mat(t6k   , t6k   ) =  1.0;
      sys_mat(t6k +1, t6k +1) =  1.0;
      sys_mat(t6k +2, t6k +2) =  1.0;
    }

    for (int k=0; k<3; k++)
    {
      const int tk = 3*k;
      const int t3k = tk + 3*3;
      const int t9k = 3*(9+k);
      sys_mat(t9k  , tk   ) = -0.5;
      sys_mat(t9k+1, tk+1 ) = -0.5;
      sys_mat(t9k+2, tk+2 ) = -0.5;
      sys_mat(t9k  , t3k  ) = -0.5;
      sys_mat(t9k+1, t3k+1) = -0.5;
      sys_mat(t9k+2, t3k+2) = -0.5;
      sys_mat(t9k  , t9k  ) =  1.0;
      sys_mat(t9k+1, t9k+1) =  1.0;
      sys_mat(t9k+2, t9k+2) =  1.0;
    }

    for (int k=0; k<3; k++)
    {
      const int t3k = 3*(3+k);
      const int t12k = 3*(12+k);
      const int a = 3 * (3 + ((k < 2) ? (k+1) : 0));
      sys_mat(t12k   , t3k    ) = -0.5;
      sys_mat(t12k +1, t3k +1 ) = -0.5;
      sys_mat(t12k +2, t3k +2 ) = -0.5;
      sys_mat(t12k   , a      ) = -0.5;
      sys_mat(t12k +1, a   +1 ) = -0.5;
      sys_mat(t12k +2, a   +2 ) = -0.5;
      sys_mat(t12k   , t12k   ) =  1.0;
      sys_mat(t12k +1, t12k +1) =  1.0;
      sys_mat(t12k +2, t12k +2) =  1.0;
    }

    ale3_tors_spring_wedge6(sys_mat,xyze);
    break;

  case DRT::Element::wedge6:
    ale3_tors_spring_wedge6(sys_mat,xyze);
    break;


  case DRT::Element::hex20:

    for (int k=0; k<4; k++)
    {
      const int tk = 3*k;
      const int t8k = tk + 3*8;
      const int a = 3 * ((k < 3) ? (k+1) : 0);
      sys_mat(t8k   , tk    ) = -0.5;
      sys_mat(t8k +1, tk +1 ) = -0.5;
      sys_mat(t8k +2, tk +2 ) = -0.5;
      sys_mat(t8k   , a     ) = -0.5;
      sys_mat(t8k +1, a +1  ) = -0.5;
      sys_mat(t8k +2, a +2  ) = -0.5;
      sys_mat(t8k   , t8k   ) =  1.0;
      sys_mat(t8k +1, t8k +1) =  1.0;
      sys_mat(t8k +2, t8k +2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      const int tk = 3*k;
      const int t4k = tk + 3*4;
      const int t12k = 3*(12+k);
      sys_mat(t12k  , tk    ) = -0.5;
      sys_mat(t12k+1, tk+1  ) = -0.5;
      sys_mat(t12k+2, tk+2  ) = -0.5;
      sys_mat(t12k  , t4k   ) = -0.5;
      sys_mat(t12k+1, t4k+1 ) = -0.5;
      sys_mat(t12k+2, t4k+2 ) = -0.5;
      sys_mat(t12k  , t12k  ) =  1.0;
      sys_mat(t12k+1, t12k+1) =  1.0;
      sys_mat(t12k+2, t12k+2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      const int t4k = 3*(4+k);
      const int t16k = 3*(16+k);
      const int a = 3 * ( 4 + ((k < 3) ? (k+1) : 0));
      sys_mat(t16k   , t4k    ) = -0.5;
      sys_mat(t16k +1, t4k +1 ) = -0.5;
      sys_mat(t16k +2, t4k +2 ) = -0.5;
      sys_mat(t16k   , a      ) = -0.5;
      sys_mat(t16k +1, a + 1  ) = -0.5;
      sys_mat(t16k +2, a + 2  ) = -0.5;
      sys_mat(t16k   , t16k   ) =  1.0;
      sys_mat(t16k +1, t16k +1) =  1.0;
      sys_mat(t16k +2, t16k +2) =  1.0;
    }

    ale3_tors_spring_hex8(sys_mat,xyze);
    break;

  case DRT::Element::hex27:

    for (int k=0; k<4; k++)
    {
      const int tk = 3*k;
      const int t8k = tk + 3*8;
      const int a = 3 * ((k < 3) ? (k+1) : 0);
      sys_mat(t8k   , tk    ) = -0.5;
      sys_mat(t8k +1, tk +1 ) = -0.5;
      sys_mat(t8k +2, tk +2 ) = -0.5;
      sys_mat(t8k   , a    ) = -0.5;
      sys_mat(t8k +1, a +1 ) = -0.5;
      sys_mat(t8k +2, a +2 ) = -0.5;
      sys_mat(t8k   , t8k   ) =  1.0;
      sys_mat(t8k +1, t8k +1) =  1.0;
      sys_mat(t8k +2, t8k +2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      const int tk = 3*k;
      const int t4k = tk + 3*4;
      const int t12k = tk + 3*12;
      sys_mat(t12k  , tk    ) = -0.5;
      sys_mat(t12k+1, tk+1  ) = -0.5;
      sys_mat(t12k+2, tk+2  ) = -0.5;
      sys_mat(t12k  , t4k   ) = -0.5;
      sys_mat(t12k+1, t4k+1 ) = -0.5;
      sys_mat(t12k+2, t4k+2 ) = -0.5;
      sys_mat(t12k  , t12k  ) =  1.0;
      sys_mat(t12k+1, t12k+1) =  1.0;
      sys_mat(t12k+2, t12k+2) =  1.0;
    }

    for (int k=0; k<4; k++)
    {
      const int t4k = 3*(4+k);
      const int t16k = 3*(16+k);
      const int a = 3 * (4 + ((k < 3) ? (k+1) : 0));
      sys_mat(t16k   , t4k    ) = -0.5;
      sys_mat(t16k +1, t4k +1 ) = -0.5;
      sys_mat(t16k +2, t4k +2 ) = -0.5;
      sys_mat(t16k   , a      ) = -0.5;
      sys_mat(t16k +1, a +1   ) = -0.5;
      sys_mat(t16k +2, a +2   ) = -0.5;
      sys_mat(t16k   , t16k   ) =  1.0;
      sys_mat(t16k +1, t16k +1) =  1.0;
      sys_mat(t16k +2, t16k +2) =  1.0;
    }

    sys_mat(3*20  , 3*8   ) = -0.5;
    sys_mat(3*20+1, 3*8+1 ) = -0.5;
    sys_mat(3*20+2, 3*8+2 ) = -0.5;
    sys_mat(3*20  , 3*10  ) = -0.5;
    sys_mat(3*20+1, 3*10+1) = -0.5;
    sys_mat(3*20+2, 3*10+2) = -0.5;
    sys_mat(3*20  , 3*20  ) =  1.0;
    sys_mat(3*20+1, 3*20+1) =  1.0;
    sys_mat(3*20+2, 3*20+2) =  1.0;

    for (int k=0; k<4; k++)
    {
      const int t8k = 3*(8+k);
      const int t16k = 3*(16+k);
      const int t21k = 3*(21+k);
      sys_mat(t21k   , t8k    ) = -0.5;
      sys_mat(t21k +1, t8k + 1) = -0.5;
      sys_mat(t21k +2, t8k + 2) = -0.5;
      sys_mat(t21k   , t16k   ) = -0.5;
      sys_mat(t21k +1, t16k+ 1) = -0.5;
      sys_mat(t21k +2, t16k+ 2) = -0.5;
      sys_mat(t21k   , t21k   ) =  1.0;
      sys_mat(t21k +1, t21k +1) =  1.0;
      sys_mat(t21k +2, t21k +2) =  1.0;
    }

    sys_mat(3*25  , 3*16  ) = -0.5;
    sys_mat(3*25+1, 3*16+1) = -0.5;
    sys_mat(3*25+2, 3*16+2) = -0.5;
    sys_mat(3*25  , 3*18  ) = -0.5;
    sys_mat(3*25+1, 3*18+1) = -0.5;
    sys_mat(3*25+2, 3*18+2) = -0.5;
    sys_mat(3*25  , 3*25  ) =  1.0;
    sys_mat(3*25+1, 3*25+1) =  1.0;
    sys_mat(3*25+2, 3*25+2) =  1.0;

    sys_mat(3*26  , 3*21  ) = -0.5;
    sys_mat(3*26+1, 3*21+1) = -0.5;
    sys_mat(3*26+2, 3*21+2) = -0.5;
    sys_mat(3*26  , 3*23  ) = -0.5;
    sys_mat(3*26+1, 3*23+1) = -0.5;
    sys_mat(3*26+2, 3*23+2) = -0.5;
    sys_mat(3*26  , 3*26  ) =  1.0;
    sys_mat(3*26+1, 3*26+1) =  1.0;
    sys_mat(3*26+2, 3*26+2) =  1.0;

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

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Ale3_Impl<distype>::static_ke(
  Ale3*                      ele,
  DRT::Discretization&       dis,
  vector<int>&               lm,
  Epetra_SerialDenseMatrix&  sys_mat_epetra,
  Epetra_SerialDenseVector&  /*residual*/,
  bool                       incremental,
  std::vector<double>&       my_dispnp,
  RefCountPtr<MAT::Material> material,
  ParameterList&             params)
{
  const int nd  = 3 * iel;
  // A view to sys_mat_epetra
  LINALG::Matrix<nd,nd> sys_mat(sys_mat_epetra.A(),true);

  //  get material using class StVenantKirchhoff
  if (material->MaterialType()!=INPAR::MAT::m_stvenant)
    dserror("stvenant material expected but got type %d", material->MaterialType());
  MAT::StVenantKirchhoff* actmat = static_cast<MAT::StVenantKirchhoff*>(material.get());

  LINALG::Matrix<3,iel> xyze;

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for(int i=0;i<iel;i++)
  {
    const double* x = nodes[i]->X();
    xyze(0,i)=x[0];
    xyze(1,i)=x[1];
    xyze(2,i)=x[2];
  }

  if (incremental)
  {
    for(int i=0;i<iel;i++)
    {
      xyze(0,i) += my_dispnp[3*i+0];
      xyze(1,i) += my_dispnp[3*i+1];
      xyze(2,i) += my_dispnp[3*i+2];
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);
  LINALG::Matrix<iel,1  >               weights(iel);

  if(distype==DRT::Element::nurbs8
     ||
     distype==DRT::Element::nurbs27)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

    bool zero_size=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());

    if(zero_size)
    {
      return;
    }

    for (int inode=0; inode<iel; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (ele->Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  /*----------------------------------------- declaration of variables ---*/
  LINALG::Matrix<iel,1  > funct;
  LINALG::Matrix<3,  iel> deriv;
  LINALG::Matrix<3,  3  > xjm;
  LINALG::Matrix<3,  3  > xji;
  LINALG::Matrix<6,  nd > bop;
  LINALG::Matrix<6,  6  > D(true);

  double vol=0.;

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule();
  const IntegrationPoints3D  intpoints(gaussrule);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];


    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs8
       &&
       distype != DRT::Element::nurbs27)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
      DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(3);
      gp(0)=e1;
      gp(1)=e2;
      gp(2)=e3;

      DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
	(funct  ,
	 deriv  ,
	 gp     ,
	 myknots,
	 weights,
	 distype);
    }

    // compute jacobian matrix

    // determine jacobian at point r,s,t
    xjm.MultiplyNT(deriv,xyze);

    // determinant and inverse of jacobian
    const double det = xji.Invert(xjm);

    // calculate element volume
    const double fac = intpoints.qwgt[iquad]*det;
    vol += fac;

    // get operator b of global derivatives
    for (int i=0; i<iel; i++)
    {
      const int node_start = i*3;

      // [h1,h2,h3] is the i-th column in (xji * deriv), but that
      // matrix is not calculated, because it is never needed.
      const double hr = deriv(0,i);
      const double hs = deriv(1,i);
      const double ht = deriv(2,i);

      const double h1 = xji(0,0)*hr + xji(0,1)*hs + xji(0,2)*ht;
      const double h2 = xji(1,0)*hr + xji(1,1)*hs + xji(1,2)*ht;
      const double h3 = xji(2,0)*hr + xji(2,1)*hs + xji(2,2)*ht;

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
    actmat->SetupCmat(D);

    // elastic stiffness matrix ke
    //ale3_keku(estif,bop,D,fac,nd);

    // Again what we really have here is a matrix multiplication. For
    // each j db is the j-th column of D * bop, so the whole
    // calculation is sys_mat += bop^T * D * bop.
    // (If you got the impression (or know) that this description is
    // wrong, don't believe me. I just read the code and tried to
    // figure out what it does.)
    for (int j=0; j<nd; j++) {
      double dum;
      double db[6];
      for (int k=0; k<6; k++) {
        db[k] = 0.0;
        for (int l=0; l<6; l++)
          db[k] += D(k,l)*bop(l,j)*fac ;
      }
      for (int i=0; i<nd; i++) {
        dum = 0.0;
        for (int m=0; m<6; m++)
          dum = dum + bop(m,i)*db[m] ;
        sys_mat(i,j) += dum;
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


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Ale3_Impl<distype>::static_ke_laplace(
    Ale3*                      ele,
    DRT::Discretization&       dis     ,
    Epetra_SerialDenseMatrix&  sys_mat_epetra ,
    Epetra_SerialDenseVector&  residual,
    bool                       incremental,
    std::vector<double>&        my_dispnp,
    RefCountPtr<MAT::Material>  material,
    ParameterList&              params  )
    {
  // ******************************************************
  // this method was copied from the ALE2 element and extended to 3D
  // ToDo: proper implementation of min_Jac for stiffness heuristics
  // ******************************************************

  const int nd  = 3 * iel;
  // A view to sys_mat_epetra
  LINALG::Matrix<nd,nd> sys_mat(sys_mat_epetra.A(),true);

  //  get material using class StVenantKirchhoff
  //  if (material->MaterialType()!=INPAR::MAT::m_stvenant)
  //    dserror("stvenant material expected but got type %d", material->MaterialType());
  //  MAT::StVenantKirchhoff* actmat = static_cast<MAT::StVenantKirchhoff*>(material.get());

  LINALG::Matrix<3,iel> xyze;

  // get node coordinates
  DRT::Node** nodes = ele->Nodes();
  for(int i=0;i<iel;i++)
  {
    const double* x = nodes[i]->X();
    xyze(0,i)=x[0];
    xyze(1,i)=x[1];
    xyze(2,i)=x[2];
  }

  if (incremental)  // Laplace with incremental????
  {
    for(int i=0;i<iel;i++)
    {
      xyze(0,i) += my_dispnp[3*i+0];
      xyze(1,i) += my_dispnp[3*i+1];
      xyze(2,i) += my_dispnp[3*i+2];
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);
  LINALG::Matrix<iel,1  >               weights(iel);

  if(distype==DRT::Element::nurbs8
      ||
      distype==DRT::Element::nurbs27)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
    =
        dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

    bool zero_size=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());

    if(zero_size)
    {
      return;
    }

    for (int inode=0; inode<iel; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
      =
          dynamic_cast<DRT::NURBS::ControlPoint* > (ele->Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  /*----------------------------------------- declaration of variables ---*/
  LINALG::Matrix<iel,1  > funct(true);
  LINALG::Matrix<3,  iel> deriv(true);
  LINALG::Matrix<3,  3  > xjm(true);
  LINALG::Matrix<3,  3  > xji(true);
  LINALG::Matrix<3,  iel> deriv_xy(true);
  LINALG::Matrix<iel,iel> tempmat(true);

  double vol=0.;

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule();
  const IntegrationPoints3D  intpoints(gaussrule);

  // This whole method was copied from the ALE2 element and extended to 3D.
  // ToDo: proper computation and usage of min_detF. Is there any detailed literature
  //       on this approach?

  // double min_detF(0.0);         /* minimal Jacobian determinant   */
  // ale2_min_jaco(Shape(),xyze,&min_detF);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];


    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs8
        &&
        distype != DRT::Element::nurbs27)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
      DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(3);
      gp(0)=e1;
      gp(1)=e2;
      gp(2)=e3;

      DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
      (funct  ,
          deriv  ,
          gp     ,
          myknots,
          weights,
          distype);
    }

    // determine jacobian matrix at point r,s,t
    xjm.MultiplyNT(deriv,xyze);

    // determinant and inverse of jacobian
    const double det = xji.Invert(xjm);

    // calculate element volume
    const double fac = intpoints.qwgt[iquad]*det;
    vol += fac;

    // compute global derivatives
    deriv_xy.Multiply(xji,deriv);

    /*------------------------- diffusivity depends on displacement ---*/
    //   This is how it is done in the 2d implementation ALE2:
    //   const double k_diff = 1.0/min_detF/min_detF;

    //   This is how we do it here for the time being due to lack of detailed knowledge
    //   on the underlying concept of min_detF (see comments above). We simply
    //   use here the Jacobi determinant evaluated at the Gauss points instead of
    //   min_detF determined at corner node positions.
    const double k_diff = 1.0/(det*det);
    //   Due to this heuristic, small elements are artificially made stiffer
    //   and large elements are made softer
    /*------------------------------- sort it into stiffness matrix ---*/

    tempmat.MultiplyTN(fac*k_diff,deriv_xy,deriv_xy,1.0);

  } // integration loop

  // insert finished temporary matrix
  for (int d=0; d < 3; d++)
  {
    for (int i=0; i < iel; i++)
    {
      for (int j=0; j < iel; j++)
      {
        sys_mat(i*3+d,j*3+d)+= tempmat(i,j);
      }
    }
  }

  return;
    }


// get optimal gaussrule for discretization type
template <DRT::Element::DiscretizationType distype>
inline GaussRule3D DRT::ELEMENTS::Ale3_Impl<distype>::getOptimalGaussrule()
{
    switch (distype)
    {
    case DRT::Element::hex8:
    case DRT::Element::nurbs8:
      return intrule_hex_8point;
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::nurbs27:
      return intrule_hex_27point;
    case DRT::Element::tet4:
      return intrule_tet_4point;
    case DRT::Element::tet10:
      return intrule_tet_5point;
    case DRT::Element::wedge6:
      return intrule_wedge_6point;
    case DRT::Element::wedge15:
      return intrule_wedge_9point;
    case DRT::Element::pyramid5:
      return intrule_pyramid_8point;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
      return intrule3D_undefined;
    }
}


#endif
