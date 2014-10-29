/*----------------------------------------------------------------------------*/
/*!
\file ale2_evaluate.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289 15362
</pre>
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
  else if (action == "calc_ale_lin")
     act = Ale2::calc_ale_lin;
  else if (action == "calc_ale_lin_stiff")
    act = Ale2::calc_ale_lin_fixed_ref;
  else if (action == "calc_ale_solid")
     act = Ale2::calc_ale_solid;
  else if (action == "calc_ale_laplace")
      act = Ale2::calc_ale_laplace;
  else if (action == "calc_ale_spring" || action == "calc_ale_springs" )
    act = Ale2::calc_ale_springs;
  else if (action == "calc_ale_spring_fixed_ref")
    act = Ale2::calc_ale_springs_fixed_ref;
  else if (action == "setup_material")
    act = Ale2::setup_material;
  else
    dserror("%s is an unknown type of action for Ale2",action.c_str());

  // get the material
  Teuchos::RCP<MAT::Material> mat = Material();

  switch (act)
  {
    case calc_ale_lin:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke(discretization,lm,&elemat1,elevec1,true,my_dispnp,params);

      break;
    }
    case calc_ale_lin_fixed_ref:
    {
      std::vector<double> my_dispnp(lm.size(),0.0);
      static_ke(discretization,lm,&elemat1,elevec1,false,my_dispnp,params);

      break;
    }
    case calc_ale_solid:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke_nonlinear(discretization,lm,&elemat1,elevec1,my_dispnp,params);

      break;
    }
    case calc_ale_laplace:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);
      static_ke_laplace(discretization,lm,&elemat1,elevec1,my_dispnp,params);

      break;
    }
    case calc_ale_springs:
    {
      Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp"); // get the displacements
      std::vector<double> my_dispnp(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,my_dispnp,lm);

      static_ke_spring(&elemat1,elevec1,true,my_dispnp);

      break;
    }
    case calc_ale_springs_fixed_ref:
    {
      // same as calc_ale_spring, however, no displ. and hence initial/reference configuration
      // is used for stiffness matrix computation
      std::vector<double> my_dispnp(lm.size(),0.0);
      static_ke_spring(&elemat1,elevec1,false,my_dispnp);

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
    Epetra_SerialDenseVector& residual, bool incremental,
    std::vector<double>& displacements)
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

  if(incremental)
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
  residual.Scale(0.0);
  for(int i =0; i< 2*iel; ++i)
  {
    for(int j=0; j<2*iel; ++j)
    {
      residual[i]+=(*sys_mat)(i,j)*displacements[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::static_ke(DRT::Discretization& dis,
    std::vector<int> &lm, Epetra_SerialDenseMatrix *sys_mat,
    Epetra_SerialDenseVector& residual, bool incremental,
    std::vector<double>& my_dispnp, Teuchos::ParameterList &params)
{
  const int iel = NumNode();
  const int nd  = 2 * iel;
  const DiscretizationType distype = this->Shape();

  // get material using class StVenantKirchhoff
  if (Material()->MaterialType()!=INPAR::MAT::m_stvenant)
    dserror("stvenant material expected but got type %d", Material()->MaterialType());
  MAT::StVenantKirchhoff* actmat = static_cast<MAT::StVenantKirchhoff*>(Material().get());

  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  if (incremental)
  {
    for(int i=0;i<iel;i++)
    {
      xyze(0,i) += my_dispnp[2*i+0];
      xyze(1,i) += my_dispnp[2*i+1];
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
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

    bool zero_sized=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());

    if(zero_sized)
    {
      return;
    }

    for (int inode=0; inode<iel; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  /*----------------------------------------- declaration of variables ---*/
  Epetra_SerialDenseVector funct(iel);
  Epetra_SerialDenseMatrix deriv(2,iel);
  Epetra_SerialDenseMatrix xjm(2,2);
  Epetra_SerialDenseMatrix xji(2,2);
  Epetra_SerialDenseMatrix bop(3,2*iel);
  Epetra_SerialDenseMatrix d(4,4);

  // gaussian points
  const DRT::UTILS::GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs4
       &&
       distype != DRT::Element::nurbs9)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct  ,
         deriv  ,
         gp     ,
         myknots,
         weights,
         distype);
    }

    // compute jacobian matrix

    // determine jacobian at point r,s,t
    for (int i=0; i<2; i++)
    {
      for (int j=0; j<2; j++)
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
    const double det = xjm(0,0)*xjm(1,1) - xjm(0,1)*xjm(1,0);
    if(abs(det)<1.0e-6)
    {
      dserror("det %12.5e in ale element %d\n",det,Id());
    }
    const double fac = intpoints.qwgt[iquad]*det;

    // calculate operator B

    // inverse of jacobian
    const double dum=1.0/det;
    xji(0,0) = xjm(1,1)* dum;
    xji(0,1) =-xjm(0,1)* dum;
    xji(1,0) =-xjm(1,0)* dum;
    xji(1,1) = xjm(0,0)* dum;

    // get operator b of global derivatives
    for (int i=0; i<iel; i++)
    {
      const int node_start = i*2;

      const double hr   = deriv(0,i);
      const double hs   = deriv(1,i);

      const double h1 = xji(0,0)*hr + xji(0,1)*hs;
      const double h2 = xji(1,0)*hr + xji(1,1)*hs;
/*
           | Nk,x    0   |
           |   0    Nk,y |, k=0...iel-1
           | Nk,y   Nk,x |
*/

      bop(0,node_start+0) = h1 ;
      bop(0,node_start+1) = 0.0;
      bop(1,node_start+0) = 0.0;
      bop(1,node_start+1) = h2 ;
      bop(2,node_start+0) = h2;
      bop(2,node_start+1) = h1;
    }

    // call material law
    actmat->SetupCmat2d(&d);

    for (int j=0; j<nd; j++)
    {
      double db[3];
      for (int k=0; k<3; k++)
      {
        db[k] = 0.0;
        for (int l=0; l<3; l++)
        {
          db[k] += d(k,l)*bop(l,j)*fac ;
        }
      }
      for (int i=0; i<nd; i++)
      {
        double dum = 0.0;
        for (int m=0; m<3; m++)
        {
          dum = dum + bop(m,i)*db[m] ;
        }
        (*sys_mat)(i,j) += dum ;
      }
    }
  }
  residual.Scale(0.0);
  for(int i =0; i< 2*iel; ++i)
  {
    for(int j=0; j<2*iel; ++j)
    {
      residual[i]+=(*sys_mat)(i,j)*my_dispnp[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::static_ke_nonlinear(
    DRT::Discretization& dis,
    const std::vector<int> &lm,
    Epetra_SerialDenseMatrix *sys_mat,
    Epetra_SerialDenseVector& residual,
    const std::vector<double>& my_dispnp,
    Teuchos::ParameterList &params)
{
  // declaration of variables
  const DiscretizationType distype = this->Shape();
  const int numnode = NumNode(); // number of nodes
  const int numdf   = 2; // number of degrees of freedom
  const int nd      = numnode * numdf; // number of nodal degrees of freedom
  double det        = 0.0; // jacobian determinant wrt. xrefe
  Epetra_SerialDenseVector funct(numnode); // shape functions
  Epetra_SerialDenseMatrix deriv(2, numnode); // shape function derivatives
  Epetra_SerialDenseMatrix xjm(2, 2); // jacobian matrix wrt. xrefe
  Epetra_SerialDenseMatrix xji(2, 2); // inverse of jacobian matrix wrt. xrefe
  Epetra_SerialDenseMatrix boplin(4, nd); // blin operator
  Epetra_SerialDenseMatrix b_cure(4, nd); // blin operator in current configuration
  Epetra_SerialDenseMatrix F(2, 2); // deformation gradient
  Epetra_SerialDenseMatrix strain(2, 2); // strains
  Epetra_SerialDenseMatrix stress(4, 4); // stresses
  Epetra_SerialDenseMatrix C(4, 4); // constitutive matrix
  Epetra_SerialDenseMatrix stiffmatrix(nd, nd); // stiffness matrix
  Epetra_SerialDenseVector force(nd); // force vector

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(2);
  Epetra_SerialDenseVector              weights(numnode);

  if(distype==DRT::Element::nurbs4
     ||
     distype==DRT::Element::nurbs9)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

    bool zero_sized=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());

    if(zero_sized)
    {
      return;
    }

    for (int inode=0; inode<numnode; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  // get element node coordinates and calculate current configuration
  Epetra_SerialDenseMatrix xrefe(2, numnode); // reference configuration
  for (int k = 0; k < numnode; ++k)
  {
    xrefe(0, k) = Nodes()[k]->X()[0];
    xrefe(1, k) = Nodes()[k]->X()[1];
  }

  // gaussian points
  const DRT::UTILS::GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  // integration loop
  for (int ip = 0; ip < intpoints.nquad; ++ip)
  {
    // get gausspoints and integration weight
    const double e1 = intpoints.qxg[ip][0];
    const double e2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs4
       &&
       distype != DRT::Element::nurbs9)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct  ,
         deriv  ,
         gp     ,
         myknots,
         weights,
         distype);
    }

    // compute jacobian matrix wrt. xrefe
    ale2_jacobianmatrix(xrefe, deriv, xjm, &det, numnode);
    double fac = wgt * det;

    // calculate shape function derivatives wrt. xrefe
    Epetra_SerialDenseMatrix bop_temp(2, numnode);
    ale2_bop(bop_temp, deriv, xjm, xji, numnode);

    // fill blin operator
    for (int j = 0; j < numnode; j++)
    {
      boplin(0, 2 * j)     = bop_temp(0, j);
      boplin(1, 2 * j + 1) = bop_temp(1, j);
      boplin(2, 2 * j)     = bop_temp(1, j);
      boplin(3, 2 * j + 1) = bop_temp(0, j);
    }

    // calculate deformation gradient F and Green-Lagrange deformations
    ale2_greenlag(F, strain, deriv, xji, my_dispnp);

    // write the deformation gradient F in matrix notation
    Epetra_SerialDenseMatrix Fmatrix(4, 4);
    Fmatrix(0, 0) = F(0, 0);
    Fmatrix(0, 2) = 0.5 * F(0, 1);
    Fmatrix(0, 3) = 0.5 * F(0, 1);
    Fmatrix(1, 1) = F(1, 1);
    Fmatrix(1, 2) = 0.5 * F(1, 0);
    Fmatrix(1, 3) = 0.5 * F(1, 0);
    Fmatrix(2, 1) = F(0, 1);
    Fmatrix(2, 2) = 0.5 * F(0, 0);
    Fmatrix(2, 3) = 0.5 * F(0, 0);
    Fmatrix(3, 0) = F(1, 0);
    Fmatrix(3, 2) = 0.5 * F(1, 1);
    Fmatrix(3, 3) = 0.5 * F(1, 1);

    // calculate b_cure operator
    int err = b_cure.Multiply('T', 'N', 1.0, Fmatrix, boplin, 0.0);
    if (err != 0)
      dserror("Multiply failed");

    // make 3D equivalent of Green-Lagrange strain
    LINALG::Matrix<6, 1> gl(false);
    gl(0) = strain(0, 0); // E_{11}
    gl(1) = strain(1, 1); // E_{22}
    gl(2) = 0.0; // E_{33}
    gl(3) = strain(0, 1) + strain(1, 0); // 2*E_{12}=E_{12}+E_{21}
    gl(4) = 0.0; // 2*E_{23}
    gl(5) = 0.0; // 2*E_{31}

    // call 3D stress response under plain strain
    LINALG::Matrix<6, 1> pk2(true); // Piola-Kirchhoff 2 stresses (must be zerofied!)
    LINALG::Matrix<6, 6> cmat(true); // constitutive matrix (must be zerofied!)
    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<
        MAT::So3Material>(Material(), true);
    so3mat->Evaluate(NULL, &gl, params, &pk2, &cmat, Id());

    // transform 2nd Piola-Kirchhoff stress back to 2d stress matrix
    memset(stress.A(), 0, stress.M() * stress.N() * sizeof(double)); // set all values to zero
    stress(0, 0) = stress(3, 3) = pk2(0); // S_{11}
    stress(1, 1) = stress(2, 2) = pk2(1); // S_{22}
    stress(0, 2) = stress(1, 3) = stress(3, 1) = stress(2, 0) = pk2(3); // S_{12}

    // transform elasticity matrix  back to 2d matrix
    C(0, 0) = cmat(0, 0); // C_{1111}
    C(0, 1) = cmat(0, 1); // C_{1122}
    C(0, 2) = cmat(0, 3); // C_{1112}
    C(0, 3) = cmat(0, 3); // C_{1112} = C_{1121}

    C(1, 0) = cmat(1, 0); // C_{2211}
    C(1, 1) = cmat(1, 1); // C_{2222}
    C(1, 2) = cmat(1, 3); // C_{2212}
    C(1, 3) = cmat(1, 3); // C_{2221} = C_{2212}

    C(2, 0) = cmat(3, 0); // C_{1211}
    C(2, 1) = cmat(3, 1); // C_{1222}
    C(2, 2) = cmat(3, 3); // C_{1212}
    C(2, 3) = cmat(3, 3); // C_{1221} = C_{1212}

    C(3, 0) = cmat(3, 0); // C_{2111} = C_{1211}
    C(3, 1) = cmat(3, 1); // C_{2122} = C_{1222}
    C(3, 2) = cmat(3, 3); // C_{2112} = C_{1212}
    C(3, 3) = cmat(3, 3); // C_{2121} = C_{1212}

    // geometric part of stiffness matrix k_g perform B^T * SIGMA * B * fac
    Epetra_SerialDenseMatrix dumdum(4, nd); // temporary variable
    err = dumdum.Multiply('N', 'N', 1.0, stress, boplin, 0.0);
    if (err != 0)
      dserror("Multiply failed");
    err = (*sys_mat).Multiply('T', 'N', fac, boplin, dumdum, 1.0);
    if (err != 0)
      dserror("Multiply failed");

    // linear stiffness matrix k_eu perform B_cure^T * C * B_cure * fac
    err = dumdum.Multiply('N', 'N', 1.0, C, b_cure, 0.0);
    if (err != 0)
      dserror("Multiply failed");
    err = (*sys_mat).Multiply('T', 'N', fac, b_cure, dumdum, 1.0);
    if (err != 0)
      dserror("Multiply failed");

    // compute residual forces
    Epetra_SerialDenseMatrix stress_vec(4, 1); // stress vector
    stress_vec(0, 0) = stress(0, 0);
    stress_vec(1, 0) = stress(1, 1);
    stress_vec(2, 0) = stress(0, 2);
    stress_vec(3, 0) = stress(0, 2);

    // perform B_cure^T * S * fac
    err = (residual).Multiply('T', 'N', fac, b_cure, stress_vec, 1.0);
    if (err != 0)
      dserror("Multiply failed");
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
static void ale2_min_jaco(DRT::Element::DiscretizationType distyp,
    Epetra_SerialDenseMatrix xyz, double *min_detF)
{
  double  detF[4]; // Jacobian determinant at nodes

  switch (distyp)
  {
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      /*--------------------- evaluate Jacobian determinant at nodes ---*/
      detF[0] = 0.25 * ( (xyz[0][0]-xyz[0][1]) * (xyz[1][0]-xyz[1][3])
                       - (xyz[1][0]-xyz[1][1]) * (xyz[0][0]-xyz[0][3]) );
      detF[1] = 0.25 * ( (xyz[0][0]-xyz[0][1]) * (xyz[1][1]-xyz[1][2])
                       - (xyz[1][0]-xyz[1][1]) * (xyz[0][1]-xyz[0][2]) );
      detF[2] = 0.25 * ( (xyz[0][3]-xyz[0][2]) * (xyz[1][1]-xyz[1][2])
                       - (xyz[1][3]-xyz[1][2]) * (xyz[0][1]-xyz[0][2]) );
      detF[3] = 0.25 * ( (xyz[0][3]-xyz[0][2]) * (xyz[1][0]-xyz[1][3])
                       - (xyz[1][3]-xyz[1][2]) * (xyz[0][0]-xyz[0][3]) );

      std::cout << "detF[0] = " << detF[0] << std::endl;
      std::cout << "detF[1] = " << detF[1] << std::endl;
      std::cout << "detF[2] = " << detF[2] << std::endl;
      std::cout << "detF[3] = " << detF[3] << std::endl;

      /*------------------------------------------------- check sign ---*/
      if (detF[0] <= 0.0) dserror("Negative JACOBIAN ");
      if (detF[1] <= 0.0) dserror("Negative JACOBIAN ");
      if (detF[2] <= 0.0) dserror("Negative JACOBIAN ");
      if (detF[3] <= 0.0) dserror("Negative JACOBIAN ");
      /*-------------------------------------- look for the smallest ---*/
      *min_detF = ( detF[0]  < detF[1]) ?  detF[0]  : detF[1];
      *min_detF = (*min_detF < detF[2]) ? *min_detF : detF[2];
      *min_detF = (*min_detF < detF[3]) ? *min_detF : detF[3];
      /*----------------------------------------------------------------*/
      break;
    case DRT::Element::tri3:
      *min_detF = (-xyz[0][0]+xyz[0][1]) * (-xyz[1][0]+xyz[1][2])
                - (-xyz[0][0]+xyz[0][2]) * (-xyz[1][0]+xyz[1][1]);
      if (*min_detF <= 0.0) dserror("Negative JACOBIAN ");
      break;
    default:
      dserror("minimal Jacobian determinant for this distyp not implemented");
      break;
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::static_ke_laplace(DRT::Discretization& dis,
    std::vector<int> &lm, Epetra_SerialDenseMatrix *sys_mat,
    Epetra_SerialDenseVector& residual, std::vector<double>& displacements,
    Teuchos::ParameterList &params)
{
  dserror("We don't know what is really done in the element evaluation"
      "of the Laplace smoothing strategy. Check this CAREFULLY before"
      "using it.");

  const int iel = NumNode();
//  const int nd  = 2 * iel;
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseMatrix xyze(2,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=Nodes()[i]->X()[0];
    xyze(1,i)=Nodes()[i]->X()[1];
  }

  // update spatial configuration
  for(int i=0;i<iel;i++)
  {
    xyze(0,i) += displacements[2*i+0];
    xyze(1,i) += displacements[2*i+1];
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
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix      deriv(2,iel);
  Epetra_SerialDenseMatrix      deriv_xy(2,iel);
  Epetra_SerialDenseMatrix      xjm(2,2);
  Epetra_SerialDenseMatrix      xji(2,2);
//  Epetra_SerialDenseMatrix      bop(3,2*iel);
//  Epetra_SerialDenseMatrix      d(4,4);

  // gaussian points
  const DRT::UTILS::GaussRule2D gaussrule = getOptimalGaussrule(distype);
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule);
  double             min_detF = 0.0;         /* minimal Jacobian determinant   */
  ale2_min_jaco(Shape(),xyze,&min_detF);

  // integration loops
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs4
       &&
       distype != DRT::Element::nurbs9)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct  ,
         deriv  ,
         gp     ,
         myknots,
         weights,
         distype);
    }

    // compute jacobian matrix

    // determine jacobian at point r,s,t
    for (int i=0; i<2; i++)
    {
      for (int j=0; j<2; j++)
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
    const double det = xjm(0,0)*xjm(1,1) - xjm(0,1)*xjm(1,0);
    const double fac = intpoints.qwgt[iquad]*det;

    // calculate operator B

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
    const double k_diff = 1.0/min_detF/min_detF;
    /*------------------------------- sort it into stiffness matrix ---*/
    for (int i=0; i<iel; i++)
    {
       for (int j=0; j<iel; j++)
       {
         (*sys_mat)(i*2,j*2)     += ( deriv_xy(0,i) * deriv_xy(0,j)
                                    + deriv_xy(1,i) * deriv_xy(1,j) )*fac*k_diff;
         (*sys_mat)(i*2+1,j*2+1) += ( deriv_xy(0,i) * deriv_xy(0,j)
                                    + deriv_xy(1,i) * deriv_xy(1,j) )*fac*k_diff;
       }
    }
  }
  residual.Scale(0.0);
  for(int i =0; i< 2*iel; ++i)
  {
    for(int j=0; j<2*iel; ++j)
    {
      residual[i]+=(*sys_mat)(i,j)*displacements[j];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_jacobianmatrix(
    const Epetra_SerialDenseMatrix& xrefe,
    const Epetra_SerialDenseMatrix& deriv, Epetra_SerialDenseMatrix& xjm,
    double* det, const int iel)
{
  /* calculate jacobian matrix (parameter space to reference configuration)
   * xjm = deriv * xrefe^T */
  int err = xjm.Multiply('N', 'T', 1.0, deriv, xrefe, 0.0);
  if (err != 0)
    dserror("Multiply failed");

  // determinant of jacobian
  *det = xjm(0, 0) * xjm(1, 1) - xjm(0, 1) * xjm(1, 0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_jacobianmatrix_cure(
    const Epetra_SerialDenseMatrix& xcure,
    const Epetra_SerialDenseMatrix& deriv, Epetra_SerialDenseMatrix& xjm_cure,
    double* det_cure, const int iel)
{
  /* calculate jacobian matrix (parameter space to reference configuration)
   * xjm_cure = deriv * xcure^T */
  int err = xjm_cure.Multiply('N', 'T', 1.0, deriv, xcure, 0.0);
  if (err != 0)
    dserror("Multiply failed");

  // determinant of jacobian
  *det_cure = xjm_cure(0, 0) * xjm_cure(1, 1) - xjm_cure(0, 1) * xjm_cure(1, 0);
  if (*det_cure <= 0.0)
  {
    dserror("JACOBIAN determinant ZERO or NEGATIVE in ALE element %d: %12.5e\n",
        Id(), *det_cure);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_bop(Epetra_SerialDenseMatrix& bop,
    Epetra_SerialDenseMatrix& deriv, Epetra_SerialDenseMatrix& xjm,
    Epetra_SerialDenseMatrix& xji, const int numnode)
{

  /*
   | Nk,x |
   | Nk,y |   k=0...numnode-1
   */

  // calculate inverse of jacobian
  const double dum = 1.0 / (xjm(0, 0) * xjm(1, 1) - xjm(0, 1) * xjm(1, 0));
  xji(0, 0) = xjm(1, 1) * dum;
  xji(0, 1) = -xjm(0, 1) * dum;
  xji(1, 0) = -xjm(1, 0) * dum;
  xji(1, 1) = xjm(0, 0) * dum;

  // calculate bop = xji * deriv
  int err = bop.Multiply('N', 'N', 1.0, xji, deriv, 0.0);
  if (err != 0)
    dserror("Multiply failed");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Ale2::ale2_greenlag(
    Epetra_SerialDenseMatrix& defgrad,
    Epetra_SerialDenseMatrix& greenlag,
    Epetra_SerialDenseMatrix& deriv,
    Epetra_SerialDenseMatrix& xji,
    const std::vector<double>& my_dispnp)
{
  Epetra_SerialDenseMatrix displacements(NumNode(), 2); // displacements
  Epetra_SerialDenseMatrix rcg(2, 2); // right Cauchy-Green tensor

  // rewrite displacements
  for (int i = 0; i < NumNode(); i++)
    for (int j = 0; j < 2; j++) {
      displacements(i, j) = my_dispnp[2*i+j];
    }

  // determine deformation gradient F dumdum = xji * deriv
  Epetra_SerialDenseMatrix dumdum(2, NumNode());
  int err = dumdum.Multiply('N', 'N', 1.0, xji, deriv, 0.0);
  if (err != 0)
    dserror("Multiply failed");

  // F = I + displacements^T * (xji * deriv)^T = I + displacements^T * dumdum^T
  err = defgrad.Multiply('T', 'T', 1.0, displacements, dumdum, 0.0);
  if (err != 0)
    dserror("Multiply failed");
  defgrad(0, 0) += 1.0;
  defgrad(1, 1) += 1.0;

  // calculate right Cauchy-Green tensor C = F^T * F
  err = rcg.Multiply('T', 'N', 1.0, defgrad, defgrad, 0.0);
  if (err != 0)
    dserror("Multiply failed");

  // calculate Green-Lagrange strain tensor E = 0.5 * (F^T * F - I)
  greenlag(0, 0) = 0.5 * (rcg(0, 0) - 1.0);
  greenlag(0, 1) = 0.5 * rcg(0, 1);
  greenlag(1, 0) = 0.5 * rcg(1, 0);
  greenlag(1, 1) = 0.5 * (rcg(1, 1) - 1.0);

  return;
}

