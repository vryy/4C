/*----------------------------------------------------------------------------*/
/*!
\file ale2_evaluate.cpp

\brief Evaluate routines of ALE element for 2D case

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "ale2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/elasthyper.H"


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale2::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Ale2::ActionType act = Ale2::none;

  // get the action required
  std::string action = params.get<std::string>("action","none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_ale_solid")
     act = Ale2::calc_ale_solid;
  else if (action == "calc_ale_solid_linear")
     act = Ale2::calc_ale_solid_linear;
  else if (action == "calc_ale_laplace_material")
    act = Ale2::calc_ale_laplace_material;
  else if (action == "calc_ale_laplace_spatial")
    act = Ale2::calc_ale_laplace_spatial;
  else if (action == "calc_ale_springs_material" )
    act = Ale2::calc_ale_springs_material;
  else if (action == "calc_ale_springs_spatial" )
    act = Ale2::calc_ale_springs_spatial;
  else if (action == "setup_material")
    act = Ale2::setup_material;
  else if (action == "calc_jacobian_determinant")
    act = Ale2::calc_det_jac;
  else
    dserror("%s is an unknown type of action for Ale2",action.c_str());

  bool spatialconfiguration = true;
  if (params.isParameter("use spatial configuration"))
    spatialconfiguration = params.get<bool>("use spatial configuration");


  // get the material
  Teuchos::RCP<MAT::Material> mat = Material();

  switch (act)
  {
    case calc_ale_solid:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke_nonlinear(lm, my_dispnp, &elemat1, &elevec1, params, true, false);

      break;
    }
    case calc_ale_solid_linear:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke_nonlinear(lm, my_dispnp, &elemat1, &elevec1, params, spatialconfiguration, true);

      break;
    }
    case calc_ale_laplace_material:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, my_dispnp, lm);
      static_ke_laplace(discretization, lm, &elemat1, elevec1, my_dispnp,
          spatialconfiguration);

      break;
    }
    case calc_ale_laplace_spatial:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp, my_dispnp, lm);
      static_ke_laplace(discretization, lm, &elemat1, elevec1, my_dispnp, true);

      break;
    }
    case calc_ale_springs_material:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp"); // get the displacements
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke_spring(&elemat1, elevec1, my_dispnp, spatialconfiguration);

      break;
    }
    case calc_ale_springs_spatial:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp"); // get the displacements
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke_spring(&elemat1, elevec1, my_dispnp, true);

      break;
    }
    case setup_material:
    {
      // get material
      Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<
          MAT::So3Material>(mat, true);

      if (so3mat->MaterialType() != INPAR::MAT::m_elasthyper
          and so3mat->MaterialType() != INPAR::MAT::m_stvenant) // ToDo (mayr): allow only materials without history
      {
        dserror("Illegal material type for ALE. Only materials allowed that do "
            "not store Gauss point data and do not need additional data from the "
            "element line definition.");
      }

      if (so3mat->MaterialType() == INPAR::MAT::m_elasthyper)
      {
        so3mat = Teuchos::rcp_dynamic_cast<MAT::ElastHyper>(mat, true);
        so3mat->Setup(0, NULL);
      }
      break; // no setup for St-Venant / classic_lin required
    }
    case calc_det_jac:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      compute_det_jac(elevec1, lm, my_dispnp);

      break;
    }
    default:
    {
      dserror("Unknown type of action for Ale2");
      break;
    }
  }

  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale2::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition,
    std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::edge_geometry(int i, int j,
    const Epetra_SerialDenseMatrix& xyze, double* length, double* sin_alpha,
    double* cos_alpha)
{
  double delta_x, delta_y;
  /*---------------------------------------------- x- and y-difference ---*/
  delta_x = xyze(0,j)-xyze(0,i);
  delta_y = xyze(1,j)-xyze(1,i);
  /*------------------------------- determine distance between i and j ---*/
  *length = sqrt( delta_x * delta_x
                + delta_y * delta_y );
  if (*length < (1.0E-14)) dserror("edge or diagonal of element has zero length");
  /*--------------------------------------- determine direction of i-j ---*/
  *sin_alpha = delta_y / *length;
  *cos_alpha = delta_x / *length;
  /*----------------------------------------------------------------------*/
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double DRT::ELEMENTS::Ale2::ale2_area_tria(const Epetra_SerialDenseMatrix& xyze,
    int i, int j, int k)
{
  double a, b, c;  /* geometrical values */
  double el_area;  /* element area */
  /*----------------------------------------------------------------------*/

  a = (xyze(0,i)-xyze(0,j))*(xyze(0,i)-xyze(0,j))
     +(xyze(1,i)-xyze(1,j))*(xyze(1,i)-xyze(1,j)); /* line i-j squared */
  b = (xyze(0,j)-xyze(0,k))*(xyze(0,j)-xyze(0,k))
     +(xyze(1,j)-xyze(1,k))*(xyze(1,j)-xyze(1,k)); /* line j-k squared */
  c = (xyze(0,k)-xyze(0,i))*(xyze(0,k)-xyze(0,i))
     +(xyze(1,k)-xyze(1,i))*(xyze(1,k)-xyze(1,i)); /* line k-i squared */
  el_area = 0.25 * sqrt(2.0*a*b + 2.0*b*c + 2.0*c*a - a*a - b*b - c*c);
  return el_area;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_torsional(int i, int j, int k,
    const Epetra_SerialDenseMatrix& xyze, Epetra_SerialDenseMatrix* k_torsion)
{
/*
                           k
                           *
                    / \
       y,v ^      l_ki /   \  l_jk
           |           /     \
      --->     i *-------* j
            x,u          l_ij
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
  area = ale2_area_tria(xyze,i,j,k);

/*---------------------------------- determine torsional stiffnesses ---*/
  C(0,0) = l_ij*l_ij * l_ki*l_ki / (4.0*area*area);
  C(1,1) = l_ij*l_ij * l_jk*l_jk / (4.0*area*area);
  C(2,2) = l_ki*l_ki * l_jk*l_jk / (4.0*area*area);

  C(0,1) = C(0,2) = C(1,0) = C(1,2) = C(2,0) = C(2,1) = 0;

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


  int err = A.Multiply('T','N',1,R,C,0);  // A = R^t * C
  if (err!=0)
    dserror("Multiply failed");
  err = k_torsion->Multiply('N','N',1,A,R,0);  // stiff = A * R
  if (err!=0)
    dserror("Multiply failed");

}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_tors_spring_tri3(
    Epetra_SerialDenseMatrix* sys_mat, const Epetra_SerialDenseMatrix& xyze)
{
  int i, j;      /* counters */

  Epetra_SerialDenseMatrix k_tria(6,6);  // rotational stiffness matrix of one triangle

  /*-------------------------------- evaluate torsional stiffness part ---*/
  ale2_torsional(0,1,2,xyze,&k_tria);

  /*-------------------------- add everything to the element stiffness ---*/
  for (i=0; i<6; i++)
    for (j=0; j<6; j++)
      (*sys_mat)(i,j) += k_tria(i,j);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_tors_spring_quad4(
    Epetra_SerialDenseMatrix* sys_mat, const Epetra_SerialDenseMatrix& xyze)
{
  int i, j;                 /* counters */

  Epetra_SerialDenseMatrix k_tria(6,6);  // rotational stiffness matrix of one triangle

  /*--- pass all nodes and determine the triangle defined by the node i and
  adjunct nodes... ---*/
  ale2_torsional(0,1,2,xyze,&k_tria);

  /*---------- ...sort everything into the element stiffness matrix... ---*/
  for (i=0; i<6; i++)
    for (j=0; j<6; j++)
      (*sys_mat)(i,j) += k_tria(i,j);

  /*--------------------------------- ...repeat for second triangle... ---*/
  ale2_torsional(1,2,3,xyze,&k_tria);
  for (i=2; i<8; i++)
    for (j=2; j<8; j++)
      (*sys_mat)(i,j) += k_tria(i-2,j-2);

  /*------------------------------------------ ...and for the third... ---*/
  ale2_torsional(2,3,0,xyze,&k_tria);
  for (i=4; i<8; i++)
    for (j=4; j<8; j++)
      (*sys_mat)(i,j) += k_tria(i-4,j-4);
  for (i=0; i<2; i++)
    for (j=0; j<2; j++)
      (*sys_mat)(i,j) += k_tria(i+4,j+4);
  for (i=4; i<8; i++)
  {
    for (j=0; j<2; j++)
    {
      (*sys_mat)(i,j) += k_tria(i-4,j+4);
      (*sys_mat)(j,i) += k_tria(i-4,j+4);
    }
  }

  /*------------------------------- ...and eventually for a forth time ---*/
  ale2_torsional(3,0,1,xyze,&k_tria);
  for (i=0; i<4; i++)
    for (j=0; j<4; j++)
      (*sys_mat)(i,j) += k_tria(i+2,j+2);
  for (i=6; i<8; i++)
    for (j=6; j<8; j++)
      (*sys_mat)(i,j) += k_tria(i-6,j-6);
  for (i=6; i<8; i++)
  {
    for (j=0; j<4; j++)
    {
      (*sys_mat)(i,j) += k_tria(i-6,j+2);
      (*sys_mat)(j,i) += k_tria(i-6,j+2);
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::static_ke_spring(Epetra_SerialDenseMatrix* sys_mat,
    Epetra_SerialDenseVector& residual, std::vector<double>& displacements,
    const bool spatialconfiguration)
{
  const int iel = NumNode(); // numnode of this element
  const DiscretizationType distype = this->Shape();
  int numcnd; // number of corner nodes
  int node_i, node_j; // end nodes of actual spring
  double length; // length of actual edge
  double sin, cos; // direction of actual edge
  double factor;


  //number of corner nodes
  switch (distype)
  {
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
    numcnd = 4;
    break;
  case DRT::Element::tri3:
  case DRT::Element::tri6:
    numcnd = 3;
    break;
  default:
    numcnd = 0;
    dserror("distype unkown");
    break;
  }

  //Actual element coordinates
  Epetra_SerialDenseMatrix xyze(2,iel);

  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  // compute spatial configuration (if necessary)
  if (spatialconfiguration)
  {
    for(int i=0;i<iel;i++)
    {
      xyze(0,i)+= displacements[i*2];
      xyze(1,i)+= displacements[i*2+1];
    }
  }

  //Linear springs from all corner nodes to all corner nodes
  //loop over all edges and diagonals of the element
  for (node_i=0; node_i<numcnd; node_i++)
  {
    for (node_j=node_i+1; node_j<numcnd; node_j++)
    {
      edge_geometry(node_i,node_j,xyze,&length,&sin,&cos);
      factor = 1.0 / length;
      //put values in 'element stiffness'
      (*sys_mat)(node_i*2,  node_i*2  ) += cos*cos * factor;
      (*sys_mat)(node_i*2+1,node_i*2+1) += sin*sin * factor;
      (*sys_mat)(node_i*2,  node_i*2+1) += sin*cos * factor;
      (*sys_mat)(node_i*2+1,node_i*2  ) += sin*cos * factor;

      (*sys_mat)(node_j*2,  node_j*2)   += cos*cos * factor;
      (*sys_mat)(node_j*2+1,node_j*2+1) += sin*sin * factor;
      (*sys_mat)(node_j*2,  node_j*2+1) += sin*cos * factor;
      (*sys_mat)(node_j*2+1,node_j*2)   += sin*cos * factor;

      (*sys_mat)(node_i*2,  node_j*2)   -= cos*cos * factor;
      (*sys_mat)(node_i*2+1,node_j*2+1) -= sin*sin * factor;
      (*sys_mat)(node_i*2,  node_j*2+1) -= sin*cos * factor;
      (*sys_mat)(node_i*2+1,node_j*2)   -= sin*cos * factor;

      (*sys_mat)(node_j*2,  node_i*2)   -= cos*cos * factor;
      (*sys_mat)(node_j*2+1,node_i*2+1) -= sin*sin * factor;
      (*sys_mat)(node_j*2,  node_i*2+1) -= sin*cos * factor;
      (*sys_mat)(node_j*2+1,node_i*2)   -= sin*cos * factor;
    }
  }

  //build in torsional springs
  //and put edge nodes on the middle of the respective edge
  switch (distype)
  {
    case DRT::Element::quad8:
      (*sys_mat)(8,8) =  1.0;
      (*sys_mat)(8,0) = -0.5;
      (*sys_mat)(8,2) = -0.5;
      (*sys_mat)(9,9) =  1.0;
      (*sys_mat)(9,1) = -0.5;
      (*sys_mat)(9,3) = -0.5;
      (*sys_mat)(10,10) =  1.0;
      (*sys_mat)(10,2)  = -0.5;
      (*sys_mat)(10,4)  = -0.5;
      (*sys_mat)(11,11) =  1.0;
      (*sys_mat)(11,3)  = -0.5;
      (*sys_mat)(11,5)  = -0.5;
      (*sys_mat)(12,12) =  1.0;
      (*sys_mat)(12,4)  = -0.5;
      (*sys_mat)(12,6)  = -0.5;
      (*sys_mat)(13,13) =  1.0;
      (*sys_mat)(13,5)  = -0.5;
      (*sys_mat)(13,7)  = -0.5;
      (*sys_mat)(14,14) =  1.0;
      (*sys_mat)(14,6)  = -0.5;
      (*sys_mat)(14,0)  = -0.5;
      (*sys_mat)(15,15) =  1.0;
      (*sys_mat)(15,7)  = -0.5;
      (*sys_mat)(15,1)  = -0.5;
      ale2_tors_spring_quad4(sys_mat,xyze);
      break;
    case DRT::Element::quad9:
      (*sys_mat)(8,8) =  1.0;
      (*sys_mat)(8,0) = -0.5;
      (*sys_mat)(8,2) = -0.5;
      (*sys_mat)(9,9) =  1.0;
      (*sys_mat)(9,1) = -0.5;
      (*sys_mat)(9,3) = -0.5;
      (*sys_mat)(10,10) =  1.0;
      (*sys_mat)(10,2)  = -0.5;
      (*sys_mat)(10,4)  = -0.5;
      (*sys_mat)(11,11) =  1.0;
      (*sys_mat)(11,3)  = -0.5;
      (*sys_mat)(11,5)  = -0.5;
      (*sys_mat)(12,12) =  1.0;
      (*sys_mat)(12,4)  = -0.5;
      (*sys_mat)(12,6)  = -0.5;
      (*sys_mat)(13,13) =  1.0;
      (*sys_mat)(13,5)  = -0.5;
      (*sys_mat)(13,7)  = -0.5;
      (*sys_mat)(14,14) =  1.0;
      (*sys_mat)(14,6)  = -0.5;
      (*sys_mat)(14,0)  = -0.5;
      (*sys_mat)(15,15) =  1.0;
      (*sys_mat)(15,7)  = -0.5;
      (*sys_mat)(15,1)  = -0.5;
      (*sys_mat)(16,16) =  1.0;
      (*sys_mat)(16,8)  = -0.5;
      (*sys_mat)(16,12) = -0.5;
      (*sys_mat)(17,17) =  1.0;
      (*sys_mat)(17,9)  = -0.5;
      (*sys_mat)(17,13) = -0.5;
      ale2_tors_spring_quad4(sys_mat,xyze);
      break;
    case DRT::Element::quad4:
      ale2_tors_spring_quad4(sys_mat,xyze);
      break;
    case DRT::Element::tri3:
      ale2_tors_spring_tri3(sys_mat,xyze);
      break;
    case DRT::Element::tri6:
      (*sys_mat)(6,6) =  1.0;
      (*sys_mat)(6,0) = -0.5;
      (*sys_mat)(6,2) = -0.5;
      (*sys_mat)(7,7) =  1.0;
      (*sys_mat)(7,1) = -0.5;
      (*sys_mat)(7,3) = -0.5;
      (*sys_mat)(8,8) =  1.0;
      (*sys_mat)(8,2) = -0.5;
      (*sys_mat)(8,4) = -0.5;
      (*sys_mat)(9,9) =  1.0;
      (*sys_mat)(9,3) = -0.5;
      (*sys_mat)(9,5) = -0.5;
      (*sys_mat)(10,10) =  1.0;
      (*sys_mat)(10,4)  = -0.5;
      (*sys_mat)(10,0)  = -0.5;
      (*sys_mat)(11,11) =  1.0;
      (*sys_mat)(11,5)  = -0.5;
      (*sys_mat)(11,1)  = -0.5;
      ale2_tors_spring_tri3(sys_mat,xyze);
      break;
    default:
      dserror("unknown distype in ale spring dynamic");
      break;
  }

  // compute residual vector
  residual.Scale(0.0);
  for(int i =0; i< 2*iel; ++i)
    for(int j=0; j<2*iel; ++j)
      residual[i]+=(*sys_mat)(i,j)*displacements[j];

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::static_ke_nonlinear(const std::vector<int> & lm,
    const std::vector<double> & disp, Epetra_SerialDenseMatrix * stiffmatrix,
    Epetra_SerialDenseVector * force, Teuchos::ParameterList& params,
    const bool spatialconfiguration, const bool pseudolinear)
{
  const int numnode = NumNode();
  const int numdf = 2;
  const int nd = numnode * numdf;

  // general arrays
  Epetra_SerialDenseVector funct(numnode);
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  Epetra_SerialDenseMatrix boplin;
  boplin.Shape(4,2*numnode);
  Epetra_SerialDenseVector F;
  F.Size(4);
  Epetra_SerialDenseVector strain;
  strain.Size(4);
  double det;
  Epetra_SerialDenseMatrix xrefe(2,numnode);
  Epetra_SerialDenseMatrix xcure(2,numnode);
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure;
  b_cure.Shape(numeps,nd);
  Epetra_SerialDenseMatrix stress;
  stress.Shape(4,4);
  Epetra_SerialDenseMatrix C;
  C.Shape(4,4);


  /*------- get integration data ---------------------------------------- */
  const DiscretizationType distype = Shape();

  // gaussian points
  const DRT::UTILS::GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<numnode; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];

    xcure(0,k) = xrefe(0,k);
    xcure(1,k) = xrefe(1,k);

    if (spatialconfiguration)
    {
      xcure(0,k) += disp[k*numdf+0];
      xcure(1,k) += disp[k*numdf+1];
    }
  }

  /*--------------------------------- get node weights for nurbs elements */
  Epetra_SerialDenseVector weights(numnode);
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnode; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }


  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    /*================================== gaussian point and weight at it */
    const double e1 = intpoints.qxg[ip][0];
    const double e2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // get values of shape functions and derivatives in the gausspoint
    if (distype != DRT::Element::nurbs4 && distype != DRT::Element::nurbs9)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      dserror("Not implemented yet!");
    }

    /*--------------------------------------- compute jacobian Matrix */
    JacobianMatrix(xrefe,deriv,xjm,&det,numnode);

    /*------------------------------------ integration factor  -------*/
    const double fac = wgt * det;
    /*----------------------------------- calculate operator Blin  ---*/
    CalcBOpLin(boplin,deriv,xjm,det,numnode);
    //cout.precision(16);
    /*------------ calculate defgrad F^u, Green-Lagrange-strain E^u --*/
    DefGrad(F,strain,xrefe,xcure,boplin,numnode);


    /*-calculate defgrad F in matrix notation and Blin in current conf.*/
    BOpLinCure(b_cure,boplin,F,numeps,nd);


    CallMatGeoNonl(strain,stress,C,numeps,Material(),params);



    /*---------------------- geometric part of stiffness matrix kg ---*/
    if (stiffmatrix) Kg(*stiffmatrix,boplin,stress,fac,nd,numeps);

    /*------------------ elastic+displacement stiffness matrix keu ---*/
    if (stiffmatrix) Keu(*stiffmatrix,b_cure,C,fac,nd,numeps);

    /*--------------- nodal forces fi from integration of stresses ---*/
    if (not pseudolinear and force) Fint(stress,b_cure,*force,fac,nd);


  } // for (int ip=0; ip<totngp; ++ip)

  if (pseudolinear and force)
  {
    Epetra_SerialDenseVector displacements;
    displacements.Resize(nd);
    for (int i = 0; i < nd; ++i)
      displacements(i) = disp[i];

    stiffmatrix->Multiply(false, displacements, *force);
  }

  return;
}

///*----------------------------------------------------------------------------*/
///*----------------------------------------------------------------------------*/
//static void ale2_min_jaco(DRT::Element::DiscretizationType distyp,
//    Epetra_SerialDenseMatrix xyz, double *min_detF)
//{
//  double  detF[4]; // Jacobian determinant at nodes
//
//  switch (distyp)
//  {
//    case DRT::Element::quad4:
//    case DRT::Element::quad8:
//    case DRT::Element::quad9:
//      /*--------------------- evaluate Jacobian determinant at nodes ---*/
//      detF[0] = 0.25 * ( (xyz[0][0]-xyz[0][1]) * (xyz[1][0]-xyz[1][3])
//                       - (xyz[1][0]-xyz[1][1]) * (xyz[0][0]-xyz[0][3]) );
//      detF[1] = 0.25 * ( (xyz[0][0]-xyz[0][1]) * (xyz[1][1]-xyz[1][2])
//                       - (xyz[1][0]-xyz[1][1]) * (xyz[0][1]-xyz[0][2]) );
//      detF[2] = 0.25 * ( (xyz[0][3]-xyz[0][2]) * (xyz[1][1]-xyz[1][2])
//                       - (xyz[1][3]-xyz[1][2]) * (xyz[0][1]-xyz[0][2]) );
//      detF[3] = 0.25 * ( (xyz[0][3]-xyz[0][2]) * (xyz[1][0]-xyz[1][3])
//                       - (xyz[1][3]-xyz[1][2]) * (xyz[0][0]-xyz[0][3]) );
//
//      std::cout << "detF[0] = " << detF[0] << std::endl;
//      std::cout << "detF[1] = " << detF[1] << std::endl;
//      std::cout << "detF[2] = " << detF[2] << std::endl;
//      std::cout << "detF[3] = " << detF[3] << std::endl;
//
//      /*------------------------------------------------- check sign ---*/
//      if (detF[0] <= 0.0) dserror("Negative JACOBIAN ");
//      if (detF[1] <= 0.0) dserror("Negative JACOBIAN ");
//      if (detF[2] <= 0.0) dserror("Negative JACOBIAN ");
//      if (detF[3] <= 0.0) dserror("Negative JACOBIAN ");
//      /*-------------------------------------- look for the smallest ---*/
//      *min_detF = ( detF[0]  < detF[1]) ?  detF[0]  : detF[1];
//      *min_detF = (*min_detF < detF[2]) ? *min_detF : detF[2];
//      *min_detF = (*min_detF < detF[3]) ? *min_detF : detF[3];
//      /*----------------------------------------------------------------*/
//      break;
//    case DRT::Element::tri3:
//      *min_detF = (-xyz[0][0]+xyz[0][1]) * (-xyz[1][0]+xyz[1][2])
//                - (-xyz[0][0]+xyz[0][2]) * (-xyz[1][0]+xyz[1][1]);
//      if (*min_detF <= 0.0) dserror("Negative JACOBIAN ");
//      break;
//    default:
//      dserror("minimal Jacobian determinant for this distyp not implemented");
//      break;
//  }
//  return;
//}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::static_ke_laplace(DRT::Discretization& dis,
    std::vector<int> &lm, Epetra_SerialDenseMatrix *sys_mat,
    Epetra_SerialDenseVector& residual, std::vector<double>& displacements,
    const bool spatialconfiguration)
{
//  dserror("We don't know what is really done in the element evaluation"
//      "of the Laplace smoothing strategy. Check this CAREFULLY before"
//      "using it.");

  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  // update spatial configuration if necessary
  if (spatialconfiguration)
  {
    for(int i=0;i<iel;i++)
    {
      xyze(0,i) += displacements[2*i+0];
      xyze(1,i) += displacements[2*i+1];
    }
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(2);
  Epetra_SerialDenseVector              weights(iel);

  if(distype==DRT::Element::nurbs4
     ||
     distype==DRT::Element::nurbs9)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

    (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());

    for (int inode=0; inode<iel; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        = dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector funct(iel);
  Epetra_SerialDenseMatrix deriv(2,iel);
  Epetra_SerialDenseMatrix deriv_xy(2,iel);
  Epetra_SerialDenseMatrix xjm(2,2);
  Epetra_SerialDenseMatrix xji(2,2);

  // Gauss quadrature points
  const DRT::UTILS::GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);
//  double min_detF = 0.0;         /* minimal Jacobian determinant   */
//  ale2_min_jaco(Shape(),xyze,&min_detF);

  // integration loop
  for (int iquad = 0; iquad < intpoints.nquad; ++iquad)
  {
    // parameter coordinates in 1- and 2-direction at quadrature point 'iquad'
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    // get values of shape functions and derivatives at the gausspoint
    if(distype != DRT::Element::nurbs4 && distype != DRT::Element::nurbs9)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D(funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv(funct, deriv, gp, myknots,
          weights, distype);
    }

    // compute jacobian matrix

    // determine jacobian at point r,s,t
    for (int i=0; i<2; ++i)
    {
      for (int j=0; j<2; ++j)
      {
        double dum=0.;
        for (int l=0; l<iel; ++l)
        {
          dum += deriv(i,l)*xyze(j,l);
        }
        xjm(i,j)=dum;
      }
    }

    // determinant of jacobian
    const double det = xjm(0,0)*xjm(1,1) - xjm(0,1)*xjm(1,0);
    const double fac = intpoints.qwgt[iquad]*det;

    // inverse of jacobian
    const double dum=1.0/det;
    xji(0,0) = xjm(1,1)* dum;
    xji(0,1) =-xjm(0,1)* dum;
    xji(1,0) =-xjm(1,0)* dum;
    xji(1,1) = xjm(0,0)* dum;

    for (int isd=0; isd<2; isd++)
      for (int jsd=0; jsd<2; jsd++)
        for (int inode=0; inode<iel; inode++)
          deriv_xy(isd,inode) = xji(isd,jsd) * deriv(jsd,inode);

    /*------------------------- diffusivity depends on displacement ---*/
    const double k_diff = 1.0; // 1.0/min_detF/min_detF;

    /*------------------------------- sort it into stiffness matrix ---*/
    for (int i=0; i<iel; ++i)
    {
       for (int j=0; j<iel; ++j)
       {
         (*sys_mat)(i*2,j*2)     += ( deriv_xy(0,i) * deriv_xy(0,j)
                                    + deriv_xy(1,i) * deriv_xy(1,j) )*fac*k_diff;
         (*sys_mat)(i*2+1,j*2+1) += ( deriv_xy(0,i) * deriv_xy(0,j)
                                    + deriv_xy(1,i) * deriv_xy(1,j) )*fac*k_diff;
       }
    }
  }

  residual.Scale(0.0);
  for (int i = 0; i < 2 * iel; ++i)
    for (int j = 0; j < 2 * iel; ++j)
      residual[i] += (*sys_mat)(i, j) * displacements[j];

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::CalcBOpLin(Epetra_SerialDenseMatrix& boplin,
                                     Epetra_SerialDenseMatrix& deriv,
                                     Epetra_SerialDenseMatrix& xjm,
                                     double& det,
                                     const int iel)
{

  double dum;
  double xji[2][2];
  /*---------------------------------------------- inverse of jacobian ---*/
  dum = 1.0/det;
  xji[0][0] = xjm(1,1)* dum;
  xji[0][1] =-xjm(0,1)* dum;
  xji[1][0] =-xjm(1,0)* dum;
  xji[1][1] = xjm(0,0)* dum;
  /*----------------------------- get operator boplin of global derivatives -*/
  /*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,x    0   |
       |   0    Nk,y |
       | Nk,y    0   |
       |  0     Nk,x |
  */
  for (int inode=0; inode<iel; inode++)
  {
    int dnode = inode*2;

    boplin(0,dnode+0) = deriv(0,inode)*xji[0][0] + deriv(1,inode)*xji[0][1];
    boplin(1,dnode+1) = deriv(0,inode)*xji[1][0] + deriv(1,inode)*xji[1][1];
    boplin(2,dnode+0) = boplin(1,dnode+1);
    boplin(3,dnode+1) = boplin(0,dnode+0);
  } /* end of loop over nodes */
  return;
}


/*-----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::JacobianMatrix(
  const Epetra_SerialDenseMatrix& xrefe,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix& xjm,
  double* det,
  const int iel
)
{
  memset(xjm.A(),0,xjm.N()*xjm.M()*sizeof(double));

  for (int k=0; k<iel; k++)
  {
    xjm(0,0) += deriv(0,k) * xrefe(0,k);
    xjm(0,1) += deriv(0,k) * xrefe(1,k);
    xjm(1,0) += deriv(1,k) * xrefe(0,k);
    xjm(1,1) += deriv(1,k) * xrefe(1,k);
  }

/*------------------------------------------ determinant of jacobian ---*/
  *det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];

  if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT %8.5f in ELEMENT %d\n",*det,Id());
/*----------------------------------------------------------------------*/

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::DefGrad(Epetra_SerialDenseVector& F,
                           Epetra_SerialDenseVector& strain,
                           const Epetra_SerialDenseMatrix& xrefe,
                           const Epetra_SerialDenseMatrix& xcure,
                           Epetra_SerialDenseMatrix& boplin,
                           const int iel)
{
  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,X  |
        |  1 + Uy,Y  |
        |      Ux,Y  |
        |      Uy,X  |
  */

  memset(F.A(),0,F.N()*F.M()*sizeof(double));

  F[0] = 1;
  F[1] = 1;
  for (int inode=0; inode<iel; inode++)
  {
     F[0] += boplin(0,2*inode)   * (xcure(0,inode) - xrefe(0,inode));  // F_11
     F[1] += boplin(1,2*inode+1) * (xcure(1,inode) - xrefe(1,inode));  // F_22
     F[2] += boplin(2,2*inode)   * (xcure(0,inode) - xrefe(0,inode));  // F_12
     F[3] += boplin(3,2*inode+1) * (xcure(1,inode) - xrefe(1,inode));  // F_21
  } /* end of loop over nodes */

  /*-----------------------calculate Green-Lagrange strain E -------------*/
  strain[0] = 0.5 * (F[0] * F[0] + F[3] * F[3] - 1.0);  // E_11
  strain[1] = 0.5 * (F[2] * F[2] + F[1] * F[1] - 1.0);  // E_22
  strain[2] = 0.5 * (F[0] * F[2] + F[3] * F[1]);        // E_12
  strain[3] = strain[2];                                // E_21

  return;
}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::Kg(Epetra_SerialDenseMatrix& estif,
                                 const Epetra_SerialDenseMatrix& boplin,
                                 const Epetra_SerialDenseMatrix& stress,
                                 const double fac,
                                 const int nd,
                                 const int numeps)
{
  /*---------------------------------------------- perform B^T * SIGMA * B*/
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
      for(int r=0; r<numeps; r++)
         for(int m=0; m<numeps; m++)
            estif(i,j) += boplin(r,i)*stress(r,m)*boplin(m,j)*fac;

  return;

}

/*----------------------------------------------------------------------*
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::Keu(Epetra_SerialDenseMatrix& estif,
                                  const Epetra_SerialDenseMatrix& b_cure,
                                  const Epetra_SerialDenseMatrix& C,
                                  const double fac,
                                  const int nd,
                                  const int numeps)
{

  /*------------- perform B_cure^T * D * B_cure, whereas B_cure = F^T * B */
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
        for(int k=0; k<numeps; k++)
           for(int m=0; m<numeps; m++)
             estif(i,j) +=  b_cure(k,i)*C(k,m)*b_cure(m,j)*fac;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::Fint(const Epetra_SerialDenseMatrix& stress,
                                   const Epetra_SerialDenseMatrix& b_cure,
                                   Epetra_SerialDenseVector& intforce,
                                   const double fac,
                                   const int nd)

{
  Epetra_SerialDenseVector st;
  st.Size(4);

  st[0] = fac * stress(0,0);
  st[1] = fac * stress(1,1);
  st[2] = fac * stress(0,2);
  st[3] = fac * stress(0,2);

  for(int i=0; i<nd; i++)
    for(int j=0; j<4; j++)
      intforce[i] += b_cure(j,i)*st[j];

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::CallMatGeoNonl(
  const Epetra_SerialDenseVector& strain,  ///< Green-Lagrange strain vector
  Epetra_SerialDenseMatrix& stress,  ///< stress vector
  Epetra_SerialDenseMatrix& C,  ///< elasticity matrix
  const int numeps,  ///< number of strains
  Teuchos::RCP<const MAT::Material> material,  ///< the material data
  Teuchos::ParameterList& params ///< element parameter list
)
{

  /*--------------------------- call material law -> get tangent modulus--*/
  switch(material->MaterialType())
  {
    case INPAR::MAT::m_stvenant:/*----------------------- linear elastic ---*/
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(material.get());
      double ym = actmat->Youngs();
      double pv = actmat->PoissonRatio();

  /*----------- material-tangente - plane strain, rotational symmetry ---*/

      const double c1=ym/(1.0+pv);
      const double b1=c1*pv/(1.0-2.0*pv);
      const double a1=b1+c1;

      C(0,0)=a1;
      C(0,1)=b1;
      C(0,2)=0.;
      C(0,3)=0.;

      C(1,0)=b1;
      C(1,1)=a1;
      C(1,2)=0.;
      C(1,3)=0.;

      C(2,0)=0.;
      C(2,1)=0.;
      C(2,2)=c1/2.;
      C(2,3)=c1/2.;

      C(3,0)=0.;
      C(3,1)=0.;
      C(3,2)=c1/2;
      C(3,3)=c1/2;


      /*-------------------------- evaluate 2.PK-stresses -------------------*/
      /*------------------ Summenschleife -> += (2.PK stored as vector) ------*/

      Epetra_SerialDenseVector svector;
      svector.Size(3);

      for (int k=0; k<3; k++)
      {
        for (int i=0; i<numeps; i++)
        {
          svector(k) += C(k,i) * strain(i);
        }
      }
      /*------------------ 2.PK stored as matrix -----------------------------*/
      stress(0,0)=svector(0);
      stress(0,2)=svector(2);
      stress(1,1)=svector(1);
      stress(1,3)=svector(2);
      stress(2,0)=svector(2);
      stress(2,2)=svector(1);
      stress(3,1)=svector(2);
      stress(3,3)=svector(0);


      break;
    }
    case INPAR::MAT::m_elasthyper: // general hyperelastic matrial (bborn, 06/09)
    {
      MaterialResponse3dPlane(stress,C,strain,params);
      break;
    }
    default:
    {
      dserror("Invalid type of material law for wall element");
      break;
    }
  } // switch(material->MaterialType())

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::MaterialResponse3dPlane(
  Epetra_SerialDenseMatrix& stress,
  Epetra_SerialDenseMatrix& C,
  const Epetra_SerialDenseVector& strain,
  Teuchos::ParameterList& params )
{
  // make 3d equivalent of Green-Lagrange strain
  LINALG::Matrix<6,1> gl(false);
  GreenLagrangePlane3d(strain,gl);

  // call 3d stress response
  LINALG::Matrix<6,1> pk2(true);  // must be zerofied!!!
  LINALG::Matrix<6,6> cmat(true);  // must be zerofied!!!
  MaterialResponse3d(&pk2, &cmat, &gl,params);

// we have plain strain

  // transform 2nd Piola--Kirchhoff stress back to 2d stress matrix
  memset(stress.A(), 0, stress.M()*stress.N()*sizeof(double));  // zerofy
  stress(0,0) = stress(3,3) = pk2(0);  // S_{11}
  stress(1,1) = stress(2,2) = pk2(1);  // S_{22}
  stress(0,2) = stress(1,3) = stress(3,1) = stress(2,0) = pk2(3); // S_{12}

  // transform elasticity matrix  back to 2d matrix
  C(0,0) = cmat(0,0);  // C_{1111}
  C(0,1) = cmat(0,1);  // C_{1122}
  C(0,2) = cmat(0,3);  // C_{1112}
  C(0,3) = cmat(0,3);  // C_{1112} = C_{1121}

  C(1,0) = cmat(1,0);  // C_{2211}
  C(1,1) = cmat(1,1);  // C_{2222}
  C(1,2) = cmat(1,3);  // C_{2212}
  C(1,3) = cmat(1,3);  // C_{2221} = C_{2212}

  C(2,0) = cmat(3,0);  // C_{1211}
  C(2,1) = cmat(3,1);  // C_{1222}
  C(2,2) = cmat(3,3);  // C_{1212}
  C(2,3) = cmat(3,3);  // C_{1221} = C_{1212}

  C(3,0) = cmat(3,0);  // C_{2111} = C_{1211}
  C(3,1) = cmat(3,1);  // C_{2122} = C_{1222}
  C(3,2) = cmat(3,3);  // C_{2112} = C_{1212}
  C(3,3) = cmat(3,3);  // C_{2121} = C_{1212}

  // leave this dump
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::MaterialResponse3d(
  LINALG::Matrix<6,1>* stress,
  LINALG::Matrix<6,6>* cmat,
  const LINALG::Matrix<6,1>* glstrain,
  Teuchos::ParameterList& params
  )
{

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  if(so3mat == Teuchos::null)
    dserror("cast to So3Material failed!");

  so3mat->Evaluate(NULL,glstrain,params,stress,cmat,Id());

  return;
}

/*-----------------------------------------------------------------------------*
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::GreenLagrangePlane3d(
  const Epetra_SerialDenseVector& glplane, LINALG::Matrix<6, 1>& gl3d)
{
  gl3d(0) = glplane(0);             // E_{11}
  gl3d(1) = glplane(1);             // E_{22}
  gl3d(2) = 0.0;                    // E_{33}
  gl3d(3) = glplane(2)+glplane(3);  // 2*E_{12}=E_{12}+E_{21}
  gl3d(4) = 0.0;                    // 2*E_{23}
  gl3d(5) = 0.0;                    // 2*E_{31}

  return;
}

/*-----------------------------------------------------------------------------*
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::BOpLinCure(
    Epetra_SerialDenseMatrix& b_cure,
    const Epetra_SerialDenseMatrix& boplin, const Epetra_SerialDenseVector& F,
    const int numeps, const int nd)
{


     Epetra_SerialDenseMatrix Fmatrix;
     Fmatrix.Shape(4,4);


  /*---------------------------write Vector F as a matrix Fmatrix*/

     Fmatrix(0,0) = F[0];
     Fmatrix(0,2) = 0.5 * F[2];
     Fmatrix(0,3) = 0.5 * F[2];
     Fmatrix(1,1) = F[1];
     Fmatrix(1,2) = 0.5 * F[3];
     Fmatrix(1,3) = 0.5 * F[3];
     Fmatrix(2,1) = F[2];
     Fmatrix(2,2) = 0.5 * F[0];
     Fmatrix(2,3) = 0.5 * F[0];
     Fmatrix(3,0) = F[3];
     Fmatrix(3,2) = 0.5 * F[1];
     Fmatrix(3,3) = 0.5 * F[1];

    /*-------------------------------------------------int_b_cure operator*/
      memset(b_cure.A(),0,b_cure.N()*b_cure.M()*sizeof(double));
      for(int i=0; i<numeps; i++)
        for(int j=0; j<nd; j++)
          for(int k=0; k<numeps; k++)
            b_cure(i,j) += Fmatrix(k,i)*boplin(k,j);
    /*----------------------------------------------------------------*/

  return;
}

/*-----------------------------------------------------------------------------*
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::compute_det_jac(Epetra_SerialDenseVector& elevec1,
    const std::vector<int> & lm, const std::vector<double> & disp)
{
  const int numnode = NumNode();
  const int numdf = 2;

  // general arrays
  Epetra_SerialDenseVector funct(numnode);
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  double det;
  double qm = 0.0;
  Epetra_SerialDenseMatrix xrefe(2,numnode);
  Epetra_SerialDenseMatrix xcure(2,numnode);

  Epetra_SerialDenseVector detjac(4);
  Epetra_SerialDenseVector quality(4);

  /*------- get integration data ---------------------------------------- */
  const DiscretizationType distype = Shape();
  if (distype != quad4)
    dserror("Quality metric is currently implemented for Quad4 elements, only.");

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<numnode; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];

    // We always evaluate the current configuration
    xcure(0,k) = xrefe(0,k) += disp[k*numdf+0];
    xcure(1,k) = xrefe(1,k) += disp[k*numdf+1];
  }

  // array with x- and y-coordinates of nodes in parameter space
  double nodepos[4][2];
  nodepos[0][0] = -1.0;
  nodepos[0][1] = -1.0;
  nodepos[1][0] =  1.0;
  nodepos[1][1] = -1.0;
  nodepos[2][0] =  1.0;
  nodepos[2][1] =  1.0;
  nodepos[3][0] = -1.0;
  nodepos[3][1] =  1.0;

  /*------------- Loop over all nodes -----------------------------------*/
  for (int node = 0; node < 4; ++node)
  {
    /*================================== gaussian point and weight at it */
    const double e1 = nodepos[node][0];
    const double e2 = nodepos[node][1];

    // get values of shape functions and derivatives in the gausspoint
    // shape functions and their derivatives for polynomials
    DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);

    /*--------------------------------------- compute jacobian Matrix */
    JacobianMatrix(xrefe,deriv,xjm,&det,numnode);

    /*---------------------------- Evaluate quality measure */
    EvaluateOddy(xjm,det,qm);


    detjac[node] = det;
    quality[node] = qm;
  } // loop over nodes

  // assign results
  double mindetjac = detjac[0];
  double minqm = quality[0];
  for (int i = 1; i < 4; ++i)
  {
    if (detjac[i] < mindetjac)
      mindetjac = detjac[i];
    if (quality[i] < minqm)
      minqm = quality[i];
  }

  elevec1[0] = mindetjac;
  elevec1[1] = minqm;

  return;
}

/*-----------------------------------------------------------------------------*
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::EvaluateOddy(const Epetra_SerialDenseMatrix& xjm,
    double det, double& qm)
{
  // compute C
  Epetra_SerialDenseMatrix c(2,2,true);
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      for (int k = 0; k < 2; ++k)
        c(i,j) += xjm[k][i] * xjm[k][j];
      c(i,j) /= det;
    }
  }

  // compute D
  double d = 0.0;
  for (int i = 0; i < 2; ++i)
  {
    for (int j = 0; j < 2; ++j)
    {
      d += c(i,j)*c(i,j);
    }
  }

  for (int k = 0; k < 2; ++k)
    d -= 0.5 * c(k,k);

  qm = d;
}
