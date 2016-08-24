/*!----------------------------------------------------------------------
\file so_surface_evaluate.cpp

\brief class for evaluation of equations on the structural surface

\maintainer Michael Gee

\level 1
*----------------------------------------------------------------------*/

#include "so_surface.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_surfstress/drt_surfstress_manager.H"
#include "../drt_potential/drt_potential_manager.H"
#include "Sacado.hpp"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_constraint/springdashpot_new.H"
#include "../drt_mat/structporo.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_immersed_problem/immersed_base.H"
#include "../drt_immersed_problem/immersed_field_exchange_manager.H"
#include "../drt_crack/crackUtils.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"

using UTILS::SurfStressManager;
using POTENTIAL::PotentialManager;
using UTILS::SpringDashpotNew;

/*----------------------------------------------------------------------*
 * Integrate a Surface Neumann boundary condition (public)     gee 04/08|
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::EvaluateNeumann(Teuchos::ParameterList&  params,
                                                      DRT::Discretization&     discretization,
                                                      DRT::Condition&          condition,
                                                      std::vector<int>&        lm,
                                                      Epetra_SerialDenseVector& elevec1,
                                                      Epetra_SerialDenseMatrix* elemat1)
{
  // set the interface ptr in the parent element
  ParentElement()->SetParamsInterfacePtr(params);
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
  double pstime = sdyn.get<double>("PRESTRESSTIME");

  // IMPORTANT: The 'neum_orthopressure' case represents a truly nonlinear follower-load
  // acting on the spatial configuration. Therefore, it needs to be linearized. On the
  // contrary, the simplified 'neum_pseudo_orthopressure' option allows for an approximative
  // modeling of an orthopressure load without the need to do any linearization. However,
  // this can only be achieved by referring the 'neum_pseudo_orthopressure' load to the last
  // converged configuration, which introduces an error as compared with 'neum_orthopressure'.
  bool loadlin = (elemat1!=NULL);

  // type of Neumann conditions
  enum LoadType
  {
    neum_none,
    neum_live,                  // standard Neumann load
    neum_pseudo_orthopressure,  // pseudo-orthopressure load
    neum_orthopressure,         // orthopressure load
    neum_torque                 // torque
  };

  LoadType ltype     = neum_none;

  // type of configurations
  enum Configuration
  {
    config_none,
    config_material,       // material configuration
    config_lastconverged,  // last converged configuration
    config_spatial         // spatial configuration
  };

  Configuration config = config_none;

  // get type of condition
  const std::string* type = condition.Get<std::string>("type");
  if (*type == "neum_live")
  {
    ltype   = neum_live;
    config  = config_material;
  }
  else if (*type == "neum_pseudo_orthopressure")
  {
    ltype  = neum_pseudo_orthopressure;
    config = config_lastconverged;
  }
  else if (*type == "neum_orthopressure")
  {
    ltype  = neum_orthopressure;
    config = config_spatial;
  }
  else if (*type == "neum_torque")
  {
    ltype  = neum_torque;
    config = config_spatial;
  }
  else
  {
    dserror("Unknown type of SurfaceNeumann condition");
  }

  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    spa_func = condition.Get<std::vector<int> >("funct");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  double time = -1.0;
  if (ParentElement()->IsParamsInterface())
    time = ParentElement()->ParamsInterfacePtr()->GetTotalTime();
  else
    time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  const int numdim = 3;

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < numdim)
    dserror("Fewer functions or curves defined than the element has dofs.");

  // find out whether we will use time curves and get the factors
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  std::vector<double> curvefacs(numdim, 1.0);
  for (int i=0; i < numdim; ++i)
  {
    const int curvenum = (curve) ? (*curve)[i] : -1;
    if (curvenum>=0 && usetime)
      curvefacs[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // element geometry update
  const int numnode = NumNode();
  const int numdf   = NumDofPerNode(*Nodes()[0]);
  LINALG::SerialDenseMatrix x(numnode,numdim);
  LINALG::SerialDenseMatrix xc;
  switch (config)
  {
  case config_material:
  {
    // no linearization needed for load in material configuration
    loadlin = false;

    // evaluate material configuration
    MaterialConfiguration(x);

  }
  break;
  case config_lastconverged:
  {
    // initialize last converged configuration
    xc.LightShape(numnode,numdim);

    // in inverse design analysis the spatial configuration (i.e imaged state) is the reference for
    //which we want equilibrium. This true spatial configuration is the material configuration in BACI
    // no exection here for mulf because displacements are reset after each load step, hence last converged state is always
    // the material configuration
    if (pstype==INPAR::STR::prestress_id && time <= pstime)
    {
      // no linearization needed for inverse analysis
      loadlin = false;

      // evaluate material configuration
      MaterialConfiguration(xc);
    }
    // standard forward analysis
    else
    {
      // no linearization needed for load in last converged configuration
      loadlin = false;

      // evaluate last converged configuration
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      SpatialConfiguration(xc,mydisp);
    }
  }
  break;
  case config_spatial:
  {
    // initialize spatial configuration
    xc.LightShape(numnode,numdim);

    // in inverse design analysis the spatial configuration (i.e imaged state) is the reference for
    //which we want equilibrium. This true spatial configuration is the material configuration in BACI
    // same holds for mulf
    if ((pstype==INPAR::STR::prestress_id ||pstype==INPAR::STR::prestress_mulf )&& time <= pstime)
    {
      // no linearization needed for inverse design analysis and mulf
      loadlin = false;

      // evaluate material configuration
      MaterialConfiguration(xc);
    }
    else // standard case
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement new");
      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement new'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      SpatialConfiguration(xc,mydisp);
    }
  }
  break;
  default: dserror("Unknown case of frame");
  break;
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  bool nurbsele=false;

  DRT::NURBS::NurbsDiscretization* nurbsdis =
    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if(nurbsdis!=NULL) nurbsele=true;

  // factor for surface orientation
  double normalfac=1.0;

  // local surface id
  const int surfaceid = LSurfNumber();

  // knot vectors for parent volume and this surface
  std::vector<Epetra_SerialDenseVector> mypknots(3);
  std::vector<Epetra_SerialDenseVector> myknots (2);

  // NURBS control point weights for all nodes, ie. CPs
  Epetra_SerialDenseVector weights(numnode);

  if(nurbsele)
  {
    // --------------------------------------------------
    // get knotvector
    Teuchos::RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();
    bool zero_size=knots->GetBoundaryEleAndParentKnots(mypknots     ,
                                                       myknots      ,
                                                       normalfac    ,
                                                       ParentElement()->Id(),
                                                       surfaceid    );
    // elements that have zero size in knotspan are skipped
    // (only required to enforce interpolation at certain points
    //  using repeated knots)
    if(zero_size) return 0;

    // --------------------------------------------------
    // get node weights for nurbs elements
    for (int inode=0; inode<numnode; inode++)
    {
      DRT::NURBS::ControlPoint* cp =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);
      weights(inode) = cp->W();
    }
  }
  // --------------------------------------------------


  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  funct(numnode);
  LINALG::SerialDenseMatrix  deriv(2,numnode);

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
  for (int gp=0; gp<intpoints.nquad; gp++)
  {
    // set gausspoints from integration rule
    Epetra_SerialDenseVector e(2);
    e(0) = intpoints.qxg[gp][0];
    e(1) = intpoints.qxg[gp][1];

    // get shape functions and derivatives in the plane of the element
    if(!nurbsele)
    {
      DRT::UTILS::shape_function_2D(funct,e(0),e(1),Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,e(0),e(1),Shape());
    }
    else
    {
      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct,deriv,e,myknots,weights,nurbs9);
    }

    //Stuff to get spatial Neumann
    LINALG::SerialDenseMatrix gp_coord(1,numdim);

    switch(ltype)
    {
    case neum_live:
    {

      //check for correct input
      for (int checkdof = numdim; checkdof < int(onoff->size()); ++checkdof)
      {
        if ((*onoff)[checkdof] != 0)
          dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
      }

      LINALG::SerialDenseMatrix dxyzdrs(2,3);
      dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);
      LINALG::SerialDenseMatrix  metrictensor(2,2);
      metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      const double detA = sqrt( metrictensor(0,0)*metrictensor(1,1)
                                -metrictensor(0,1)*metrictensor(1,0));

      double functfac = 1.0;
      int functnum = -1;
      double val_curvefac_functfac;

      for(int dof=0;dof<numdim;dof++)
      {
        if ((*onoff)[dof]) // is this dof activated?
        {

          //factor given by spatial function
          if (spa_func)
            functnum = (*spa_func)[dof];

          if (functnum>0)
          {
            //Calculate reference position of GP
            gp_coord.Multiply('T','N',1.0,funct,x,0.0);
            // write coordinates in another datatype
            double gp_coord2[numdim];
            for(int i=0;i<numdim;i++)
            {
              gp_coord2[i]=gp_coord(0,i);
            }
            const double* coordgpref = &gp_coord2[0]; // needed for function evaluation

            //evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,coordgpref,time,NULL);
          }
          else
            functfac = 1.0;

          val_curvefac_functfac = functfac*curvefacs[dof];
          const double fac = intpoints.qwgt[gp] * detA * (*val)[dof] * val_curvefac_functfac;
          for (int node=0; node < numnode; ++node)
          {
            elevec1[node*numdf+dof]+= funct[node] * fac;
          }
        }

      }
    }
    break;

    case neum_pseudo_orthopressure:
    case neum_orthopressure:
    {
     if ((*onoff)[0] != 1) dserror("orthopressure on 1st dof only!");
      for (int checkdof = 1; checkdof < 3; ++checkdof)
        if ((*onoff)[checkdof] != 0) dserror("orthopressure on 1st dof only!");
      double ortho_value = (*val)[0];
      //if (!ortho_value) dserror("no orthopressure value given!"); // in case of coupling with redairways,
                                                                    //there is a zero orthoval in the beginning!!!!
      std::vector<double> normal(3);
      SurfaceIntegration(normal, xc,deriv);
      //Calculate spatial position of GP
      double functfac = 1.0;
      int functnum = -1;
      double val_curvefac_functfac;

      // factor given by spatial function
      if (spa_func)
        functnum = (*spa_func)[0];

      if (functnum>0)
      {
        gp_coord.Multiply('T','N',1.0,funct,xc,0.0);
        // write coordinates in another datatype
        double gp_coord2[numdim];
        for(int i=0;i<numdim;i++)
        {
          gp_coord2[i]=gp_coord(0,i);
        }
        const double* coordgpref = &gp_coord2[0]; // needed for function evaluation

        // evaluate function at current gauss point
        functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,coordgpref,time,NULL);
      }

      val_curvefac_functfac = curvefacs[0]*functfac;

      const double fac = intpoints.qwgt[gp] * val_curvefac_functfac * ortho_value * normalfac;
      for (int node=0; node < numnode; ++node)
        for(int dim=0 ; dim<3; dim++)
          elevec1[node*numdf+dim] += funct[node] * normal[dim] * fac;

      // load linearization (if necessary)
      if (loadlin)
      {
        Epetra_SerialDenseMatrix Dnormal(numdf,numdf*numnode);
        analytical_DSurfaceIntegration(Dnormal, xc, deriv);     // --> analytical derivative
        //automatic_DSurfaceIntegration(Dnormal, xc, deriv);    // --> automatic derivative (Sacado)

        // build surface element load linearization matrix
        // (CAREFUL: Minus sign due to the fact that external forces enter the global
        // residual vector with a minus sign, too! However, the load linaerization is
        // simply added to the global tangent stiffness matrix, thus we explicitly
        // need to set the minus sign here.)
        for (int node=0; node<numnode; ++node)
          for (int dim=0 ; dim<3; dim++)
            for (int dof=0; dof<elevec1.M(); dof++)
              (*elemat1)(node*numdf+dim,dof) -= funct[node] * Dnormal(dim, dof) * fac;
      }
    }
    break;

    case neum_torque:
    {
      // check whether only first, fourth, fifth and sixth value is set
      if ((*onoff)[0] != 1) dserror("Torque value not provided!");
      if ((*onoff)[3] != 1) dserror("X-coordinate of axis for torque not provided!");
      if ((*onoff)[4] != 1) dserror("Y-coordinate of axis for torque not provided!");
      if ((*onoff)[5] != 1) dserror("Z-coordinate of axis for torque not provided!");
      for (int checkdof = 1; checkdof < 3; ++checkdof)
        if ((*onoff)[checkdof] != 0) dserror("Incorrect value for torque!");

      // get values for torque and coordinates of axis
      double torque_value = (*val)[0];
      LINALG::Matrix<3,1> axis(true);
      axis(0) = (*val)[3];
      axis(1) = (*val)[4];
      axis(2) = (*val)[5];

      // compute normal vector (with area as length)
      std::vector<double> normal(3);
      SurfaceIntegration(normal,xc,deriv);

      // compute cross product of axis and (negative) normal
      LINALG::Matrix<3,1> crossprod(true);
      crossprod(0) = -(axis(1)*normal[2]-axis(2)*normal[1]);
      crossprod(1) = -(axis(2)*normal[0]-axis(0)*normal[2]);
      crossprod(2) = -(axis(0)*normal[1]-axis(1)*normal[0]);

      // factor given by spatial function
      double functfac = 1.0;
      int functnum = -1;
      double val_curvefac_functfac;
      if (spa_func) functnum = (*spa_func)[0];
      {
        if (functnum > 0)
        {
          gp_coord.Multiply('T','N',1.0,funct,xc,0.0);
          // write coordinates in another datatype
          double gp_coord2[numdim];
          for(int i=0;i<numdim;i++)
          {
            gp_coord2[i]=gp_coord(0,i);
          }
          const double* coordgpref = &gp_coord2[0]; // needed for function evaluation

          // evaluate function at current gauss point
          functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(0,coordgpref,time,NULL);
        }
        else functfac = 1.0;
      }
      val_curvefac_functfac = curvefacs[0]*functfac;

      const double fac = intpoints.qwgt[gp] * val_curvefac_functfac * torque_value;
      for (int node=0; node < numnode; ++node)
        for(int dim=0 ; dim<3; dim++)
          elevec1[node*numdf+dim] += funct[node] * crossprod(dim) * fac;
    }
    break;

    default:
      dserror("Unknown type of SurfaceNeumann load");
    break;
    }

  } /* end of loop over integration points gp */

  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate normal at gp (private)                             gee 08/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SurfaceIntegration(std::vector<double>& normal,
                                                          const Epetra_SerialDenseMatrix& x,
                                                          const Epetra_SerialDenseMatrix& deriv)
{
  // note that the length of this normal is the area dA

  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

  normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
  normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
  normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

  return;
}

/*----------------------------------------------------------------------*
 * Evaluate sqrt of determinant of metric at gp (private)      gee 04/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::SurfaceIntegration(double& detA,
                                                          std::vector<double>& normal,
                                                          const Epetra_SerialDenseMatrix& x,
                                                          const Epetra_SerialDenseMatrix& deriv)
{

  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(2,3);
  dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

  /* compute covariant metric tensor G for surface element
  **                        | g11   g12 |
  **                    G = |           |
  **                        | g12   g22 |
  ** where (o denotes the inner product, xyz a vector)
  **
  **       dXYZ   dXYZ          dXYZ   dXYZ          dXYZ   dXYZ
  ** g11 = ---- o ----    g12 = ---- o ----    g22 = ---- o ----
  **        dr     dr            dr     ds            ds     ds
  */
  LINALG::SerialDenseMatrix metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
  detA = sqrt( metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0) );
  normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
  normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
  normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

  return;
}

/*----------------------------------------------------------------------*
 * Calculates dnormal/dx_j with Sacado DFAD                   popp 06/13|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::automatic_DSurfaceIntegration(
                                             Epetra_SerialDenseMatrix& d_normal,
                                             const Epetra_SerialDenseMatrix& x,
                                             const Epetra_SerialDenseMatrix& deriv)
{
  // some parameters
  const int numnode    = NumNode();
  const int numdim     = 3;
  const int numdof     = numdim * numnode;

  // create vectors of Sacado type
  std::vector<Sacado::Fad::DFad<double> > saccado_x(x.N()*x.M());
  std::vector<Sacado::Fad::DFad<double> > saccado_deriv(deriv.N()*deriv.M());
  std::vector<Sacado::Fad::DFad<double> > saccado_g1(numdim);
  std::vector<Sacado::Fad::DFad<double> > saccado_g2(numdim);
  std::vector<Sacado::Fad::DFad<double> > saccado_normal(numdim);

  // copy data of coordinate matrix x
  for(int row=0; row < x.M(); row++)
  {
    for(int column = 0; column < x.N(); column++)
    {
      saccado_x[x.N()*row+column] = x(row,column);
      saccado_x[x.N()*row+column].diff(x.N()*row+column,x.N()*x.M());
    }
  }

  // copy data of shape function derivatives matrix deriv
  for(int row=0; row < deriv.M(); row++)
  {
    for(int column = 0; column < deriv.N(); column++)
    {
      saccado_deriv[deriv.N()*row+column]= deriv(row,column);
    }
  }

  // re-compute local basis vectors g1 and g2
  for(int dim = 0; dim < numdim; dim++)
  {
    for(int column = 0; column < deriv.N(); column++)
    {
      saccado_g1[dim] += saccado_deriv[column]* saccado_x[column*x.N()+dim];
      saccado_g2[dim] += saccado_deriv[column+deriv.N()]* saccado_x[column*x.N()+dim];
    }
  }

  // re-compute normal vector (cross product g1 x g2)
  saccado_normal[0] = saccado_g1[1]*saccado_g2[2]-saccado_g1[2]*saccado_g2[1];
  saccado_normal[1] = saccado_g1[2]*saccado_g2[0]-saccado_g1[0]*saccado_g2[2];
  saccado_normal[2] = saccado_g1[0]*saccado_g2[1]-saccado_g1[1]*saccado_g2[0];

  // direct access to the Sacado derivatives
  for(int dim = 0; dim < numdim; dim++)
  {
    for(int dxyz = 0; dxyz< numdof; dxyz++)
    {
      d_normal(dim, dxyz) = saccado_normal[dim].fastAccessDx(dxyz);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 * Calculates dnormal/dx_j analytically                       popp 06/13|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::analytical_DSurfaceIntegration(
                                           Epetra_SerialDenseMatrix& d_normal,
                                           const Epetra_SerialDenseMatrix& x,
                                           const Epetra_SerialDenseMatrix& deriv)
{
  // some parameters
  const int numnode    = NumNode();
  const int numsurfdim = 2;
  const int numdim     = 3;
  const int numdof     = numdim * numnode;

  // compute dXYZ / drs (defining the two local basis vectors)
  Epetra_SerialDenseMatrix dxyzdrs(numsurfdim,numdim);
  dxyzdrs.Multiply('N','N',1.0,deriv, x, 0.0);

  // basis vectors (just ouf of laziness)
  std::vector<double> g1(numdim);
  std::vector<double> g2(numdim);

  for (int k=0;k<numdim;++k)
  {
    g1[k] = dxyzdrs(0,k);
    g2[k] = dxyzdrs(1,k);
  }

  // linearization of basis vectors
  Epetra_SerialDenseMatrix dg1(numdim,numdof);
  Epetra_SerialDenseMatrix dg2(numdim,numdof);

  for (int node=0;node<numnode;++node)
  {
    for (int k=0;k<numdim;++k)
    {
      dg1(k,node*numdim+k) = deriv(0,node);
      dg2(k,node*numdim+k) = deriv(1,node);
    }
  }

  // linearization of local surface normal vector
  for (int dof=0;dof<numdof;++dof)
  {
    d_normal(0,dof) = dg1(1,dof)*g2[2] + g1[1]*dg2(2,dof) - dg1(2,dof)*g2[1] - g1[2]*dg2(1,dof);
    d_normal(1,dof) = dg1(2,dof)*g2[0] + g1[2]*dg2(0,dof) - dg1(0,dof)*g2[2] - g1[0]*dg2(2,dof);
    d_normal(2,dof) = dg1(0,dof)*g2[1] + g1[0]*dg2(1,dof) - dg1(1,dof)*g2[0] - g1[1]*dg2(0,dof);
  }

  return;
}

/*----------------------------------------------------------------------*
 * Evaluate method for StructuralSurface-Elements               tk 10/07*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::Evaluate(Teuchos::ParameterList&   params,
                                               DRT::Discretization&      discretization,
                                               std::vector<int>&         lm,
                                               Epetra_SerialDenseMatrix& elematrix1,
                                               Epetra_SerialDenseMatrix& elematrix2,
                                               Epetra_SerialDenseVector& elevector1,
                                               Epetra_SerialDenseVector& elevector2,
                                               Epetra_SerialDenseVector& elevector3)
{
  // start with "none"
  DRT::ELEMENTS::StructuralSurface::ActionType act = StructuralSurface::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_constrvol")        act = StructuralSurface::calc_struct_constrvol;
  else if (action=="calc_struct_volconstrstiff")   act = StructuralSurface::calc_struct_volconstrstiff;
  else if (action=="calc_struct_monitarea")        act = StructuralSurface::calc_struct_monitarea;
  else if (action=="calc_struct_constrarea")       act = StructuralSurface::calc_struct_constrarea;
  else if (action=="calc_struct_areaconstrstiff")  act = StructuralSurface::calc_struct_areaconstrstiff;
  else if (action=="calc_init_vol")                act = StructuralSurface::calc_init_vol;
  else if (action=="calc_surfstress_stiff")        act = StructuralSurface::calc_surfstress_stiff;
  else if (action=="calc_potential_stiff")         act = StructuralSurface::calc_potential_stiff;
  else if (action=="calc_brownian_motion")         act = StructuralSurface::calc_brownian_motion;
  else if (action=="calc_brownian_motion_damping") act = StructuralSurface::calc_brownian_motion_damping;
  else if (action=="calc_struct_centerdisp")       act = StructuralSurface::calc_struct_centerdisp;
  else if (action=="calc_struct_rotation")         act = StructuralSurface::calc_struct_rotation;
  else if (action=="calc_undo_struct_rotation")    act = StructuralSurface::calc_undo_struct_rotation;
  else if (action=="calc_struct_area")             act = StructuralSurface::calc_struct_area;
  else if (action=="calc_ref_nodal_normals")       act = StructuralSurface::calc_ref_nodal_normals;
  else if (action=="calc_cur_normal_at_point")     act = StructuralSurface::calc_cur_normal_at_point;
  else if (action=="calc_cur_nodal_normals")       act = StructuralSurface::calc_cur_nodal_normals;
  else if (action=="calc_fluid_traction")          act = StructuralSurface::calc_fluid_traction;
  else if (action=="calc_ecm_traction")            act = StructuralSurface::calc_ecm_traction;
  else if (action=="mark_immersed_elements")       act = StructuralSurface::mark_immersed_elements;
  else if (action=="calc_struct_robinforcestiff")  act = StructuralSurface::calc_struct_robinforcestiff;
  else
  {
    std::cout << action << std::endl;
    dserror("Unknown type of action for StructuralSurface");
  }

  //create communicator
  const Epetra_Comm& Comm = discretization.Comm();
  // what the element has to do
  switch(act)
  {
  //gives the center displacement for SlideALE
  case calc_struct_centerdisp:
  {
    //We are not interested in ghosted elements
    if(Comm.MyPID()==Owner())
    {
      // element geometry update
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacementtotal");
      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacementtotal'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnode = NumNode();
      const int numdf=3;
      LINALG::SerialDenseMatrix xc(numnode,numdf);
      SpatialConfiguration(xc,mydisp);

      //integration of the displacements over the surface
      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector  funct(numnode);
      LINALG::SerialDenseMatrix  deriv(2,numnode);

      /*----------------------------------------------------------------------*
        |               start loop over integration points                     |
        *----------------------------------------------------------------------*/
      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

      Teuchos::RCP<const Epetra_Vector> dispincr = discretization.GetState("displacementincr");
      std::vector<double> edispincr(lm.size());
      DRT::UTILS::ExtractMyValues(*dispincr,edispincr,lm);
      elevector2[0] = 0;

      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        const double e0 = intpoints.qxg[gp][0];
        const double e1 = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
        DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

        std::vector<double> normal(3);
        double detA;
        SurfaceIntegration(detA,normal,xc,deriv);

        elevector2[0] +=  sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );

        for (int dim=0; dim<3; dim++)
        {
          for (int j=0; j<numnode; ++j)
          {
            elevector3[dim] +=  funct[j] * intpoints.qwgt[gp]
                                * edispincr[j*numdf + dim] * detA;
          }
        }

      }

    }
  }
  break;
  case calc_struct_rotation:
  {
    //We are not interested in ghosted elements
    if(Comm.MyPID()==Owner())
    {
      const double maxcoord=params.get<double>("maxcoord");
      INPAR::FSI::SlideALEProj aletype = params.get<INPAR::FSI::SlideALEProj>("aletype");
      const int numnode = NumNode();
      const int numdf = 3;
      double tol = 1.0E-5;

      //  element geometry update for time t_n
      Teuchos::RCP<const Epetra_Vector> dispn = discretization.GetState("displacementnp");
      if (dispn==Teuchos::null)
        dserror("Cannot get state vector 'displacementnp");
      std::vector<double> edispn(lm.size());
      DRT::UTILS::ExtractMyValues(*dispn,edispn,lm);
      LINALG::SerialDenseMatrix xcn(numnode,numdf);
      SpatialConfiguration(xcn,edispn);

      Teuchos::RCP<const Epetra_Vector> dispincr = discretization.GetState("displacementincr");
      if (dispn==Teuchos::null)
        dserror("Cannot get state vector 'displacementincr");
      std::vector<double> edispincr(lm.size());
      DRT::UTILS::ExtractMyValues(*dispincr,edispincr,lm);

      //integration of the displacements over the surface
      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector  funct(numnode);
      LINALG::SerialDenseMatrix  deriv(2,numnode);

      // loop over integration points
      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        const double e0 = intpoints.qxg[gp][0];
        const double e1 = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
        DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

        std::vector<double> normal(3);
        double detA;
        SurfaceIntegration(detA,normal,xcn,deriv);

        LINALG::SerialDenseVector  tangent(3);
        if (aletype==INPAR::FSI::ALEprojection_rot_z || aletype == INPAR::FSI::ALEprojection_rot_zsphere)
        {
          // compute tangential direction in xy-plane from normal
          tangent[0] = -normal[1];
          tangent[1] = normal[0];
          tangent[2] = 0.0;
        }
        else
        {
          dserror("rotation not yet implemented!");
        }

        if (tangent.Norm2() > tol)
        {
          tangent.Scale( 1.0/(tangent.Norm2()));
        }
        elevector2[0] +=  sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );

        for (int node=0; node < numnode; ++node)
        {
          double scalarprod = tangent[0]*edispincr[node*numdf] + tangent[1]*edispincr[node*numdf+1];
          if (aletype==INPAR::FSI::ALEprojection_rot_zsphere)
          {
            double circ(0.0);
            const double val = (1.0-pow(xcn(node,2)/maxcoord,2.0));
            if (val<0.0) // negative doubles can happen due to round-off errors
            {
              if (val> - EPS10) // seems to be a round-off error, we proceed assuming val=0.0
              {
                circ=0.0;
              }
              else // severe error
                dserror("Do not use sqrt() with a negative number");
            }
            else
            {
              circ = sqrt(val);
            }
            if (circ>tol)
              elevector3[0] +=  funct[node] * intpoints.qwgt[gp] * scalarprod * detA / circ;
          }
          else
            elevector3[0] +=  funct[node] * intpoints.qwgt[gp] * scalarprod * detA;
        }
      }
    }
  }
  break;
  case calc_undo_struct_rotation:
  {
    const double maxcoord=params.get<double>("maxcoord");
    INPAR::FSI::SlideALEProj aletype = params.get<INPAR::FSI::SlideALEProj>("aletype");
    const int numnode = NumNode();
    const int numdf = 3;
    double tol = 1.0E-5;

    //  element geometry update for time t_n
    Teuchos::RCP<const Epetra_Vector> dispn = discretization.GetState("displacementnp");
    if (dispn==Teuchos::null)
      dserror("Cannot get state vector 'displacementnp");
    std::vector<double> edispn(lm.size());
    DRT::UTILS::ExtractMyValues(*dispn,edispn,lm);
    LINALG::SerialDenseMatrix xcn(numnode,numdf);
    SpatialConfiguration(xcn,edispn);

    Teuchos::RCP<const Epetra_Vector> dispincr = discretization.GetState("displacementincr");
    if (dispn==Teuchos::null)
      dserror("Cannot get state vector 'displacementincr");
    std::vector<double> edispincr(lm.size());
    DRT::UTILS::ExtractMyValues(*dispincr,edispincr,lm);

    //integration of the displacements over the surface
    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);

    std::vector<double> nodalrot (numnode);
    for (int node=0; node<numnode; node++)
    {
      std::vector<double> normal(3); //normal in element center
      // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_2D(funct,0.0,0.0,Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,0.0,0.0,Shape());

      SurfaceIntegration(normal,xcn,deriv);

      LINALG::SerialDenseVector  tangent(3);
      if (aletype==INPAR::FSI::ALEprojection_rot_z || aletype == INPAR::FSI::ALEprojection_rot_zsphere)
      {
        // compute tangential direction in xy-plane from normal
        tangent[0] = -normal[1];
        tangent[1] = normal[0];
        tangent[2] = 0.0;
      }
      else
      {
        dserror("rotation not yet implemented!");
      }

      if (tangent.Norm2() > tol)
      {
        tangent.Scale(1.0/(tangent.Norm2()*Nodes()[node]->NumElement()));
      }

      if (aletype==INPAR::FSI::ALEprojection_rot_zsphere)
      {
        double circ(0.0);
        const double val = (1.0-pow(xcn(node,2)/maxcoord,2.0));
        if (val<0.0) // negative doubles can happen due to round-off errors
        {
          if (val> - EPS10) // seems to be a round-off error, we proceed assuming val=0.0
          {
            circ=0.0;
          }
          else // severe error
            dserror("Do not use sqrt() with a negative number");
        }
        else
        {
          circ = sqrt(val);
        }
        if (circ>tol)
        {
          for (int dof=0; dof<2;dof++)
            elevector1[node*numdf+dof] = tangent[dof]*circ;
        }
      }
      else
      {
        for (int dof=0; dof<2;dof++)
          elevector1[node*numdf+dof] = tangent[dof];
      }
    }
  }
  break;
  case calc_struct_constrvol:
  {
    //We are not interested in volume of ghosted elements
    if(Comm.MyPID()==Owner())
    {
      // element geometry update
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numdim = 3;
      LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
      SpatialConfiguration(xscurr,mydisp);
      //call submethod for volume evaluation and store rseult in third systemvector
      double volumeele = ComputeConstrVols(xscurr,NumNode());
      elevector3[0]= volumeele;
    }
  }
  break;
  case calc_struct_volconstrstiff:
  {
    // element geometry update
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    const int numdim =3;
    LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
    SpatialConfiguration(xscurr,mydisp);
    double volumeele;
    // first partial derivatives
    Teuchos::RCP<Epetra_SerialDenseVector> Vdiff1 = Teuchos::rcp(new Epetra_SerialDenseVector);
    // second partial derivatives
    Teuchos::RCP<Epetra_SerialDenseMatrix> Vdiff2 = Teuchos::rcp(new Epetra_SerialDenseMatrix);

    //get projection method
    Teuchos::RCP<DRT::Condition> condition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    const std::string* projtype = condition->Get<std::string>("projection");

    if (projtype != NULL)
    {
      //call submethod to compute volume and its derivatives w.r.t. to current displ.
      if (*projtype == "yz")
      {
        ComputeVolDeriv(xscurr, NumNode(),numdim*NumNode(), volumeele, Vdiff1, Vdiff2, 0, 0);
      }
      else if (*projtype == "xz")
      {
        ComputeVolDeriv(xscurr, NumNode(),numdim*NumNode(), volumeele, Vdiff1, Vdiff2, 1, 1);
      }
      else if (*projtype == "xy")
      {
        ComputeVolDeriv(xscurr, NumNode(),numdim*NumNode(), volumeele, Vdiff1, Vdiff2, 2, 2);
      }
      else
      {
        ComputeVolDeriv(xscurr, NumNode(),numdim*NumNode(), volumeele, Vdiff1, Vdiff2);
      }
    }
    else
      ComputeVolDeriv(xscurr, NumNode(),numdim*NumNode(), volumeele, Vdiff1, Vdiff2);

    //update rhs vector and corresponding column in "constraint" matrix
    elevector1 = *Vdiff1;
    elevector2 = *Vdiff1;
    elematrix1 = *Vdiff2;
    //call submethod for volume evaluation and store result in third systemvector
    elevector3[0]=volumeele;
  }
  break;
  case calc_init_vol:
  {
    if(Comm.MyPID()==Owner())
    {
      // the reference volume of the RVE (including inner
      // holes) is calculated by evaluating the following
      // surface integral:
      // V = 1/3*int(div(X))dV = 1/3*int(N*X)dA
      // with X being the reference coordinates and N the
      // normal vector of the surface element (exploiting the
      // fact that div(X)=1.0)

      // this is intended to be used in the serial case (microstructure)

      // NOTE: there must not be any holes penetrating the boundary!

      double V = params.get<double>("V0", 0.0);
      double dV = 0.0;
      const int numnode = NumNode();
      LINALG::SerialDenseMatrix x(numnode,3);
      MaterialConfiguration(x);

      // allocate vector for shape functions and matrix for derivatives
      LINALG::SerialDenseVector  funct(numnode);
      LINALG::SerialDenseMatrix  deriv(2,numnode);

      /*----------------------------------------------------------------------*
           |               start loop over integration points                     |
       *----------------------------------------------------------------------*/
      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        const double e0 = intpoints.qxg[gp][0];
        const double e1 = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
        DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

        std::vector<double> normal(3);
        double detA;
        SurfaceIntegration(detA,normal,x,deriv);
        const double fac = intpoints.qwgt[gp] * detA;

        double temp = 0.0;
        std::vector<double> X(3,0.);

        for (int i=0; i<numnode; i++)
        {
          X[0] += funct[i]*x(i,0);
          X[1] += funct[i]*x(i,1);
          X[2] += funct[i]*x(i,2);
        }

        for (int i=0;i<3;++i)
        {
          temp += normal[i]*normal[i];
        }

        if (temp<0.)
          dserror("calculation of initial volume failed in surface element");
        double absnorm = sqrt(temp);

        for (int i=0;i<3;++i)
        {
          normal[i] /= absnorm;
        }
        for (int i=0;i<3;++i)
        {
          dV += 1/3.0*fac*normal[i]*X[i];
        }
      }
      params.set("V0", V+dV);
    }
  }
  break;
  case calc_surfstress_stiff:
  {
    Teuchos::RCP<SurfStressManager> surfstressman =
        params.get<Teuchos::RCP<SurfStressManager> >("surfstr_man",Teuchos::null);

    if (surfstressman==Teuchos::null)
      dserror("No SurfStressManager in Solid3 Surface available");

    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition",Teuchos::null);
    if (cond==Teuchos::null)
      dserror("Condition not available in Solid3 Surface");

    double time = params.get<double>("total time",-1.0);
    double dt = params.get<double>("delta time",0.0);
    bool newstep = params.get<bool>("newstep", false);

    // element geometry update

    const int numnode = NumNode();
    LINALG::SerialDenseMatrix x(numnode,3);

    Teuchos::RCP<const Epetra_Vector> disn = discretization.GetState("new displacement");
    if (disn==Teuchos::null) dserror("Cannot get state vector 'new displacement'");
    std::vector<double> mydisn(lm.size());
    DRT::UTILS::ExtractMyValues(*disn,mydisn,lm);
    SpatialConfiguration(x,mydisn);

    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

    // set up matrices and parameters needed for the evaluation of current
    // interfacial area and its derivatives w.r.t. the displacements

    int ndof = 3*numnode;                                     // overall number of surface dofs
    double A;                                                 // interfacial area
    // first partial derivatives
    Teuchos::RCP<Epetra_SerialDenseVector> Adiff = Teuchos::rcp(new Epetra_SerialDenseVector);
    // second partial derivatives
    Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2 = Teuchos::rcp(new Epetra_SerialDenseMatrix);

    ComputeAreaDeriv(x, numnode, ndof, A, Adiff, Adiff2);

    if (cond->Type()==DRT::Condition::Surfactant)     // dynamic surfactant model
    {
      int curvenum = cond->GetInt("curve");
      double k1xC = cond->GetDouble("k1xCbulk");
      double k2 = cond->GetDouble("k2");
      double m1 = cond->GetDouble("m1");
      double m2 = cond->GetDouble("m2");
      double gamma_0 = cond->GetDouble("gamma_0");
      double gamma_min = cond->GetDouble("gamma_min");
      double gamma_min_eq = gamma_0-m1;
      double con_quot_max = (gamma_min_eq-gamma_min)/m2+1.;
      double con_quot_eq = (k1xC)/(k1xC+k2);

      surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, elevector1, elematrix1, this->Id(),
          time, dt, 0, 0.0, k1xC, k2, m1, m2, gamma_0, gamma_min, gamma_min_eq, con_quot_max, con_quot_eq, newstep);
    }
    else if (cond->Type()==DRT::Condition::SurfaceTension) // ideal liquid
    {
      int curvenum = cond->GetInt("curve");
      double const_gamma = cond->GetDouble("gamma");
      surfstressman->StiffnessAndInternalForces(curvenum, A, Adiff, Adiff2, elevector1, elematrix1, this->Id(),
          time, dt, 1, const_gamma, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, newstep);
    }
    else
      dserror("Unknown condition type %d",cond->Type());
  }
  break;
  // compute additional stresses due to intermolecular potential forces
  case calc_potential_stiff:
  {
    Teuchos::RCP<PotentialManager> potentialmanager =
        params.get<Teuchos::RCP<PotentialManager> >("pot_man",Teuchos::null);
    if (potentialmanager==Teuchos::null)
      dserror("No PotentialManager in Solid3 Surface available");

    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition",Teuchos::null);
    if (cond==Teuchos::null)
      dserror("Condition not available in Solid3 Surface");

    if (cond->Type()==DRT::Condition::LJ_Potential_Surface) // Lennard-Jones potential
    {
      potentialmanager->StiffnessAndInternalForcesPotential(this, gaussrule_, params, lm, elematrix1, elevector1);
    }
    if (cond->Type()==DRT::Condition::ElectroRepulsion_Potential_Surface) // Electrostatic potential
    {
      potentialmanager->StiffnessAndInternalForcesPotential(this, gaussrule_, params,lm, elematrix1, elevector1);
    }
    if (cond->Type()==DRT::Condition::VanDerWaals_Potential_Surface) // Electrostatic potential
    {
      potentialmanager->StiffnessAndInternalForcesPotential(this, gaussrule_, params,lm, elematrix1, elevector1);
    }
    if( cond->Type()!=DRT::Condition::LJ_Potential_Surface &&
        cond->Type()!=DRT::Condition::ElectroRepulsion_Potential_Surface &&
        cond->Type()!=DRT::Condition::VanDerWaals_Potential_Surface)
      dserror("Unknown condition type %d",cond->Type());
  }
  break;
  // compute stochastical forces due to Brownian Motion
  case calc_brownian_motion:
  {
    dserror("not commited");
  }
  break;
  // compute damping matrix due to Brownian Motion
  case calc_brownian_motion_damping:
  {
    dserror("not yet comitted");
  }
  break;
  //compute the area (e.g. for initialization)
  case calc_struct_monitarea:
  {
    //We are not interested in volume of ghosted elements
    if(Comm.MyPID()==Owner())
    {
      // element geometry update
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numdim = 3;
      LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
      SpatialConfiguration(xscurr,mydisp);

      Teuchos::RCP<DRT::Condition> condition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
      const std::string* projtype = condition->Get<std::string>("projection");

      // To compute monitored area consider required projection method
      // and set according coordinates to zero
      if (*projtype == "yz")
      {
        xscurr(0,0)=0;
        xscurr(1,0)=0;
        xscurr(2,0)=0;
        xscurr(3,0)=0;
      }
      else if (*projtype == "xz")
      {
        xscurr(0,1)=0;
        xscurr(1,1)=0;
        xscurr(2,1)=0;
        xscurr(3,1)=0;
      }
      else if (*projtype == "xy")
      {
        xscurr(0,2)=0;
        xscurr(1,2)=0;
        xscurr(2,2)=0;
        xscurr(3,2)=0;
      }

      double areaele=0.0;
      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
      // allocate matrix for derivatives of shape functions
      LINALG::SerialDenseMatrix  deriv(2,NumNode());

      //Compute area
      for (int gp=0; gp<intpoints.nquad; gp++)
      {
        const double e0 = intpoints.qxg[gp][0];
        const double e1 = intpoints.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

        std::vector<double> normal(3);
        double detA;
        SurfaceIntegration(detA,normal,xscurr,deriv);
        const double fac = intpoints.qwgt[gp] * detA;
        areaele += fac;

      }

      //store result in third systemvector
      elevector3[0]=areaele;
    }

  }
  break;
  case calc_struct_constrarea:
  {
    // element geometry update
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    const int numdim =3;
    LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
    SpatialConfiguration(xscurr,mydisp);
    // initialize variables
    double elearea;
    // first partial derivatives
    Teuchos::RCP<Epetra_SerialDenseVector> Adiff = Teuchos::rcp(new Epetra_SerialDenseVector);
    // second partial derivatives
    Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2 = Teuchos::rcp(new Epetra_SerialDenseMatrix);

    //call submethod
    ComputeAreaDeriv(xscurr, NumNode(),numdim*NumNode(), elearea, Adiff, Adiff2);
    // store result
    elevector3[0] = elearea;

  }
  break;
  case calc_struct_areaconstrstiff:
  {
    // element geometry update
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    const int numdim =3;
    LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
    SpatialConfiguration(xscurr,mydisp);
    // initialize variables
    double elearea;
    // first partial derivatives
    Teuchos::RCP<Epetra_SerialDenseVector> Adiff = Teuchos::rcp(new Epetra_SerialDenseVector);
    // second partial derivatives
    Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2 = Teuchos::rcp(new Epetra_SerialDenseMatrix);

    //call submethod
    ComputeAreaDeriv(xscurr, NumNode(),numdim*NumNode(), elearea, Adiff, Adiff2);
    //update elematrices and elevectors
    elevector1 = *Adiff;
    elevector1.Scale(-1.0);
    elevector2 = elevector1;
    elematrix1 = *Adiff2;
    elematrix1.Scale(-1.0);
    elevector3[0] = elearea;
  }
  break;
  case calc_struct_area:
  {
    const int numnode = NumNode();
    LINALG::SerialDenseMatrix x(numnode,3);
    MaterialConfiguration(x);
    //LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);
    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
    double a = 0.0;
    for (int gp=0; gp<intpoints.nquad; gp++)
    {
      const double e0 = intpoints.qxg[gp][0];
      const double e1 = intpoints.qxg[gp][1];
      DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());
      std::vector<double> normal(3);
      double detA;
      SurfaceIntegration(detA,normal,x,deriv);
      a += (intpoints.qwgt[gp] * detA);
    }
    double atmp = params.get("area",-1.0);
    a += atmp;
    params.set("area",a);
  }
  break;
  case calc_ref_nodal_normals:
  {
    // dummy vector
    std::vector<double> dummy(lm.size());
    BuildNormalsAtNodes(elevector1,dummy,true);
  }
  break;
  case calc_cur_nodal_normals:
  {
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    BuildNormalsAtNodes(elevector1,mydisp,false);
  }
  break;
  case calc_cur_normal_at_point:
  {
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

    const int numnode = NumNode();
    const int numdim = 3;

    LINALG::SerialDenseMatrix x(numnode,3);
    SpatialConfiguration(x,mydisp);

    const double e0 = elevector2(0);
    const double e1 = elevector2(1);

    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    double detA;
    std::vector<double> normal(3);
    SurfaceIntegration(detA,normal,x,deriv);

    for (int j=0; j<numdim; ++j)
    {
      elevector1(j) = normal[j];
    }

  }
  break;
  case calc_fluid_traction:
  {
    DRT::Problem* globalproblem = DRT::Problem::Instance();
    std::string backgrddisname(params.get<std::string>("backgrddisname"));
    std::string immerseddisname(params.get<std::string>("immerseddisname"));
    ////////////////////////////////////////////////////////////////////////////
    // rauch 05/14 originally for immersed method
    //
    // this action integrates the fluid stress over the structural surface.
    // the fluid velocities and pressure are interpolated to the structural
    // integration point externally by the method:
    //
    // InterpolateToImmersedIntPoint(...)
    //
    // which is a generally applicable method for this purpose returning
    // the desired quantitiy from another field interpolated to the current
    // integration point belonging to THIS field.
    //
    // fluiddis  = backgrddis
    // structdis = immerseddis
    //
    // if discretizations named differently are to be dealt with, provide
    // their names via the parameter list.
    //
    // CALC:   /
    //        /
    //        |
    //        | P * N * du  dA = (J * \sigma * F^-T * N, du)_{d B_0}
    //        |
    //       /
    //      / d B_0
    //
    // \sigma can be calculated from interpolated velocity gradient and pressure.
    // J is the usual jacobian determinant.
    // F is the deformation gradient of the structural surface ele.
    // N is the normal in reference configuration of the structural surface ele.
    // du denotes the virtual displacement field over the structural surface ele.
    //
    //
    /////////////////////////////////////////////////////////////////////////////

    // get integration rule
    const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad4>::rule);

    const Teuchos::RCP<DRT::Discretization> backgrddis  = globalproblem->GetDis(backgrddisname);
    const Teuchos::RCP<DRT::Discretization> immerseddis = globalproblem->GetDis(immerseddisname);

#ifdef DEBUG
    if(backgrddis == Teuchos::null)
      dserror("Pointer to background dis empty. Correct disname in parameter list 'params'?");
    if(immerseddis == Teuchos::null)
      dserror("Pointer to immersed dis empty. Correct disname in parameter list 'params'?");
#endif

    const int nen = NumNode();
    const int parent_nen = this->ParentElement()->NumNode();
    const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
    const int globdim = globalproblem->NDim();

    LINALG::Matrix<2,4> deriv;
    LINALG::Matrix<1,4> funct;
    LINALG::Matrix<1,8> parent_funct;
    LINALG::Matrix<3,8> parent_deriv;
    LINALG::Matrix<3,8> parent_deriv_notrafo;

    // get parent location matrix
    Element::LocationArray parent_la(immerseddis->NumDofSets());
    this->ParentElement()->LocationVector(*immerseddis,parent_la,false);

    // get structural state and element displacements (parent element)
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("displacement");
#ifdef DEBUG
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
#endif
    std::vector<double> parenteledisp(lm.size());
    std::vector<double> bdryeledisp(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,bdryeledisp,lm);
    DRT::UTILS::ExtractMyValues(*dispnp,parenteledisp,parent_la[0].lm_);

    // geometry (surface ele)
    LINALG::Matrix<3,4> xrefe; // material coord. of element

    DRT::Node** nodes = this->Nodes();
    for (int i=0; i<nen; ++i)
    {
      const double* x = nodes[i]->X();
      xrefe(0,i) = x[0];
      xrefe(1,i) = x[1];
      xrefe(2,i) = x[2];
    }

    // update element geometry (parent ele)
    LINALG::Matrix<3,8> parent_xrefe; // material coord. of element
    LINALG::Matrix<3,8> parent_xcurr; // current  coord. of element

    nodes = this->ParentElement()->Nodes();
    for (int i=0; i<parent_nen; ++i)
    {
      const double* x = nodes[i]->X();
      parent_xrefe(0,i) = x[0];
      parent_xrefe(1,i) = x[1];
      parent_xrefe(2,i) = x[2];

      parent_xcurr(0,i) = parent_xrefe(0,i) + parenteledisp[i*numdofpernode+0];
      parent_xcurr(1,i) = parent_xrefe(1,i) + parenteledisp[i*numdofpernode+1];
      parent_xcurr(2,i) = parent_xrefe(2,i) + parenteledisp[i*numdofpernode+2];
    }

    // get coordinates of gauss points w.r.t. local parent coordinate system
    LINALG::SerialDenseMatrix parent_xi(intpoints.IP().nquad,globdim);
    Epetra_SerialDenseMatrix derivtrafo(3,3);

    DRT::UTILS::BoundaryGPToParentGP<3>(
        parent_xi,
        derivtrafo,
        intpoints,
        DRT::Element::hex8,
        DRT::Element::quad4,
        this->FaceParentNumber());

    ////////////////////////////////////////////////////////////////////
    /////   gauss point loop
    ///////////////////////////////////////////////////////////////////
    for (int gp=0; gp<intpoints.IP().nquad; gp++)
    {
      std::vector<double> bdryxi(globdim-1);
      bdryxi[0] = intpoints.IP().qxg[gp][0];
      bdryxi[1] = intpoints.IP().qxg[gp][1];

      std::vector<double> interpolationresult (6);
      int action = (int)FLD::interpolate_velgrad_to_given_point;

      IMMERSED::InterpolateToImmersedIntPoint <DRT::Element::hex8,  // source
                                               DRT::Element::quad4> // target
                                                           (backgrddis,
                                                            immerseddis,
                                                            *this,
                                                            bdryxi,
                                                            bdryeledisp,
                                                            action,
                                                            interpolationresult // result
                                                           );

//        //////////////////
//       // Debug output //
//      //////////////////
//      std::cout<<"PROC "<<Comm.MyPID()<<": gradu and press at gp "<<gp<<" on surf ele id "<<this->Id()<<":"<<std::endl;
//      for(int i=0;i<(int)interpolationresult.size();++i)
//        std::cout<<" "<<interpolationresult[i]<<" ";
//      std::cout<<" "<<std::endl;


        // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_2D(funct,bdryxi[0],bdryxi[1],Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,bdryxi[0],bdryxi[1],Shape());

      DRT::UTILS::shape_function_3D(parent_funct,parent_xi(gp,0),parent_xi(gp,1),parent_xi(gp,2),this->ParentElement()->Shape());
      DRT::UTILS::shape_function_3D_deriv1(parent_deriv_notrafo,parent_xi(gp,0),parent_xi(gp,1),parent_xi(gp,2),this->ParentElement()->Shape());
      //parent_deriv.Multiply(derivtrafo,parent_deriv_notrafo);

//        //////////////////
//       // Debug output //
//      //////////////////
//      std::cout<<"PROC "<<Comm.MyPID()<<" : deriv at gp "<<bdryxi[0]<<" "<<bdryxi[1]<<" "<<std::endl;
//      std::cout<<"PROC "<<Comm.MyPID()<<" : "<<deriv<<std::endl;


      ////////////////////////////////////////////////////
      // calc unitnormal N in material configuration
      ///////////////////////////////////////////////////
      LINALG::Matrix<3,1> unitnormal;
      std::vector<double> normal(globdim);

      // note that the length of this normal is the area dA

      // compute dXYZ / drs
      LINALG::Matrix<2,3> dxyzdrs;
      dxyzdrs.MultiplyNT(deriv,xrefe); // to calculate unitnormal in current config. argument must be xcurr

      normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
      normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
      normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

      LINALG::Matrix<2,2> metrictensor;
      metrictensor.MultiplyNT(dxyzdrs,dxyzdrs);
      double detA = sqrt( metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0) );

//         //////////////////
//        // Debug output //
//       //////////////////
//      for(int i=0;i<globdim;++i)
//        std::cout<<" normal["<<i<<"] "<<normal[i]<<std::endl;

      for(int i=0;i<globdim;++i)
        unitnormal(i,0)=normal[i];
      const double norm2 = unitnormal.Norm2();
      unitnormal.Scale(1/norm2);

      //         //////////////////
      //        // Debug output //
      //       //////////////////
//      std::cout<<"_________________________________________________________________________________"<<std::endl;
//      std::cout<<unitnormal<<std::endl;

      /////////////////////////////////////////////////////////////////////
      // extract cauchy stress tensor from result vector of interpolation
      /////////////////////////////////////////////////////////////////////
      LINALG::Matrix<3,3> cauchystress;

      cauchystress(0,0)=interpolationresult[0];
      cauchystress(1,1)=interpolationresult[1];
      cauchystress(2,2)=interpolationresult[2];
      cauchystress(0,1)=interpolationresult[3];
      cauchystress(1,2)=interpolationresult[4];
      cauchystress(0,2)=interpolationresult[5];
      cauchystress(1,0)=cauchystress(0,1);
      cauchystress(2,1)=cauchystress(1,2);
      cauchystress(2,0)=cauchystress(0,2);

//        //////////////////
//       // Debug output //
//      //////////////////
//      std::cout<<"[ "<<cauchystress(0,0)<<" "<<cauchystress(0,1)<<" "<<cauchystress(0,2)<<" ]"<<std::endl;
//      std::cout<<"[ "<<cauchystress(1,0)<<" "<<cauchystress(1,1)<<" "<<cauchystress(1,2)<<" ]"<<std::endl;
//      std::cout<<"[ "<<cauchystress(2,0)<<" "<<cauchystress(2,1)<<" "<<cauchystress(2,2)<<" ]"<<std::endl;
//      std::cout<<"\n"<<std::endl;



      /* get the inverse of the Jacobian matrix which looks like:
      **            [ x_,r  y_,r  z_,r ]^-1
      **     J^-1 = [ x_,s  y_,s  z_,s ]
      **            [ x_,t  y_,t  z_,t ]
      */
      // compute derivatives N_XYZ at gp w.r.t. material coordinates
      // by N_XYZ = J^-1 * N_rst
      LINALG::Matrix<3,3> invJ(true);
      LINALG::Matrix<3,8> N_XYZ(true);
      invJ.MultiplyNT(parent_deriv_notrafo,parent_xrefe);
      invJ.Invert();
      N_XYZ.Multiply(invJ,parent_deriv_notrafo);

      // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
      LINALG::Matrix<3,3> defgrd_inv;
      // deformation gradient
      defgrd_inv.MultiplyNT(parent_xcurr,N_XYZ);
      // Jacobian determinant
      double J = defgrd_inv.Determinant();
      // invert deformation gradient
      defgrd_inv.Invert();

      LINALG::Matrix<3,1> tempvec;
      LINALG::Matrix<3,3> tempmat;
      tempmat.MultiplyNT(cauchystress,defgrd_inv);

//        //////////////////
//       // Debug output //
//      //////////////////
//      std::cout<<"[ "<<tempmat(0,0)<<" "<<tempmat(0,1)<<" "<<tempmat(0,2)<<" ]"<<std::endl;
//      std::cout<<"[ "<<tempmat(1,0)<<" "<<tempmat(1,1)<<" "<<tempmat(1,2)<<" ]"<<std::endl;
//      std::cout<<"[ "<<tempmat(2,0)<<" "<<tempmat(2,1)<<" "<<tempmat(2,2)<<" ]"<<std::endl;
//
//      std::cout<<"_________________________________________________________________________________"<<std::endl;

      tempvec.MultiplyNN(tempmat,unitnormal);

      double gpweight = intpoints.IP().qwgt[gp];
      // fill element vector
      for(int node = 0; node < nen; node++)
      {
        for (int dof=0;dof<globdim;dof++)
        {
          elevector1(node*numdofpernode+dof) += (gpweight*detA)*J*tempvec(dof,0)*funct(node);
        }
      }

    }// gauss point loop
  }
  break;
  case calc_ecm_traction:
  {
    DRT::Problem* globalproblem = DRT::Problem::Instance();
    std::string backgrddisname(params.get<std::string>("backgrddisname"));
    std::string immerseddisname(params.get<std::string>("immerseddisname"));
    int during_init = params.get<int>("during_init",0);
    ////////////////////////////////////////////////////////////////////////////
    // rauch 05/15 for immersed cell migration
    //
    // this action integrates the ecm interaction force for cell confinement
    // module over the cell surface.
    //
    // CALC:   /
    //        /
    //        |
    //        | f(phi) da
    //        |
    //       /
    //      / d B_t
    //
    // phi    : porosity at integration point
    // f(phi) : modeled interaction force depending on porosity phi
    //
//#define VELOCITYDEPENENCY // todo get disp_old and build scalarproduct_sip with incrementvector d_n+1 - dn !!!
//#define ECM_INTERACTION_SIMPLE
    /////////////////////////////////////////////////////////////////////////////

    // get exchange manager pointer
    DRT::ImmersedFieldExchangeManager* immersedmanager = DRT::ImmersedFieldExchangeManager::Instance();
    //Teuchos::RCP<Epetra_Vector> cell_penalty_traction = immersedmanager->GetPointerToECMPenaltyTraction();
    //Teuchos::RCP<Epetra_Vector> cell_current_nodal_normals = immersedmanager->GetPointerToCurrentNodalNormals();
    std::map<int,std::vector<LINALG::Matrix<3,1> > >* penalty_traction_at_gp = immersedmanager->GetPointerToPenaltyTractionAtGPMap();

    int isdummycall = params.get<int>("dummy_call");

    std::map<int,std::vector<double> >* porosity_at_gp_ = immersedmanager->GetPointerToPorosityAtGP();
    std::map<int,std::vector<double> >* porosity_at_gp_old = immersedmanager->GetPointerToPorosityAtGPOldTimestep();

    // DEBUG output
    //std::cout<<(porosity_at_gp_->at(this->Id()))[3]<<std::endl;

    // get integration rule
    const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad4>::rule);

    double tol = 1.0e-06;

    Teuchos::RCP<DRT::Discretization> backgrddis  = globalproblem->GetDis(backgrddisname);
    Teuchos::RCP<DRT::Discretization> immerseddis = globalproblem->GetDis(immerseddisname);

    // experimental pre cmbe
    Teuchos::RCP<DRT::Discretization> scatradis  = globalproblem->GetDis("scatra");

#ifdef DEBUG
    if(backgrddis == Teuchos::null)
      dserror("Pointer to background dis empty. Correct disname in parameter list 'params'?");
    if(immerseddis == Teuchos::null)
      dserror("Pointer to immersed dis empty. Correct disname in parameter list 'params'?");
#endif

    const int nen = NumNode();
    const int parent_nen = this->ParentElement()->NumNode();
    const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
    const int globdim = globalproblem->NDim();

    bool is_leadingedge_or_rear=false;
    //bool contains_leading_edge_node = false;

    if(not during_init)
    {
      int num_conditioned_nodes = 0;
      std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator condition_iterator;
      static std::multimap<std::string,Teuchos::RCP<DRT::Condition> > conditions = immerseddis->GetAllConditions();
      DRT::Node** bdrynodes = Nodes();

      for (condition_iterator=conditions.begin(); condition_iterator!=conditions.end(); ++condition_iterator)
      {
        if(condition_iterator->first == (std::string)"SurfaceNeumann")
        {
          for(int numnode = 0; numnode < nen; numnode++)
          {
            bool containsnode = condition_iterator->second->ContainsNode(bdrynodes[0]->Id());
            if(containsnode)
              num_conditioned_nodes++;
            else
              break;
          }
        } // if surfaceneumann condition
      } // loop over all conditions

      if(num_conditioned_nodes==nen)
      {
        is_leadingedge_or_rear=true;
      }
      else if(num_conditioned_nodes>0 and num_conditioned_nodes<nen)
      {
        //contains_leading_edge_node = true;
      }
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    } // not during init

    LINALG::Matrix<2,4> deriv;
    LINALG::Matrix<1,4> funct;
    LINALG::Matrix<1,8> parent_funct;
    LINALG::Matrix<3,8> parent_deriv;
    LINALG::Matrix<3,8> parent_deriv_notrafo;

    // get parent location matrix
    Element::LocationArray parent_la(immerseddis->NumDofSets());
    this->ParentElement()->LocationVector(*immerseddis,parent_la,false);

    // get structural state and element displacements of cell (parent element)
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("displacement");
#ifdef DEBUG
    if (dispnp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
#endif
    std::vector<double> parenteledisp(lm.size());
    std::vector<double> bdryeledisp(lm.size());
    //std::vector<double> bdryeletraction(lm.size());
    std::vector<double> bdryelenodalnormals(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,bdryeledisp,lm);
    DRT::UTILS::ExtractMyValues(*dispnp,parenteledisp,parent_la[0].lm_);
    //DRT::UTILS::ExtractMyValues(*cell_penalty_traction,bdryeletraction,lm);
    //DRT::UTILS::ExtractMyValues(*cell_current_nodal_normals,bdryelenodalnormals,lm);

    // geometry (surface ele)
    LINALG::Matrix<3,4> xrefe; // material coord. of element
    LINALG::Matrix<3,4> xcurr; // spatial  coord. of element

    DRT::Node** nodes = this->Nodes();
    for (int i=0; i<nen; ++i)
    {
      const double* x = nodes[i]->X();
      xrefe(0,i) = x[0];
      xrefe(1,i) = x[1];
      xrefe(2,i) = x[2];
    }

    for (int i=0; i<nen; ++i)
    {
      xcurr(0,i) = xrefe(0,i) + bdryeledisp[i*numdofpernode+0];
      xcurr(1,i) = xrefe(1,i) + bdryeledisp[i*numdofpernode+1];
      xcurr(2,i) = xrefe(2,i) + bdryeledisp[i*numdofpernode+2];
    }

    // update element geometry (parent ele)
    LINALG::Matrix<3,8> parent_xrefe; // material coord. of element
    LINALG::Matrix<3,8> parent_xcurr; // current  coord. of element

    nodes = this->ParentElement()->Nodes();
    for (int i=0; i<parent_nen; ++i)
    {
      const double* x = nodes[i]->X();
      parent_xrefe(0,i) = x[0];
      parent_xrefe(1,i) = x[1];
      parent_xrefe(2,i) = x[2];

      parent_xcurr(0,i) = parent_xrefe(0,i) + parenteledisp[i*numdofpernode+0];
      parent_xcurr(1,i) = parent_xrefe(1,i) + parenteledisp[i*numdofpernode+1];
      parent_xcurr(2,i) = parent_xrefe(2,i) + parenteledisp[i*numdofpernode+2];

    }

//    // DEBUG output
//    for(int lm_it=0;lm_it<lm.size();lm_it++)
//    {
//      if(lm[lm_it]>80231 and lm[lm_it]<80244)
//      {
//        std::cout<<"PROC "<<backgrddis->Comm().MyPID()<<": Bdry Ele with GID="<<this->Id()<<" owner="<<this->Owner()<<"peleid="<<this->ParentElementId()<<" peleowner="<<this->ParentElement()->Owner()<<" has dof="<<lm[lm_it]<<" at pos of node 0=["<<xcurr(0,0)<<" "<<xcurr(1,0)<<" "<<xcurr(2,0)<<"]"<<std::endl;
//      }
//    }

    // get coordinates of gauss points w.r.t. local parent coordinate system
    LINALG::SerialDenseMatrix parent_xi(intpoints.IP().nquad,globdim);
    Epetra_SerialDenseMatrix derivtrafo(3,3);

    DRT::UTILS::BoundaryGPToParentGP<3>(
        parent_xi,
        derivtrafo,
        intpoints,
        DRT::Element::hex8,
        DRT::Element::quad4,
        this->FaceParentNumber());


#ifdef VELOCITYDEPENENCY
    LINALG::Matrix<3,4> velcurr; // current velocity at t_{n+1}
    LINALG::Matrix<3,4> dispcurr; // current displacement at t_{n+1}
    if(not during_init)
    {
      // get structural state and element velocity of cell (parent element)
      Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
#ifdef DEBUG
      if (velnp==Teuchos::null) dserror("Cannot get state vector 'velnp'");
#endif

      std::vector<double> bdryelevelnp(lm.size());
      DRT::UTILS::ExtractMyValues(*velnp,bdryelevelnp,lm);

      for (int i=0; i<nen; ++i)
      {
        velcurr(0,i) = bdryelevelnp[i*numdofpernode+0];
        velcurr(1,i) = bdryelevelnp[i*numdofpernode+1];
        velcurr(2,i) = bdryelevelnp[i*numdofpernode+2];

        dispcurr(0,i) = bdryeledisp[i*numdofpernode+0];
        dispcurr(1,i) = bdryeledisp[i*numdofpernode+1];
        dispcurr(2,i) = bdryeledisp[i*numdofpernode+2];
      }
    }
#endif


    ////////////////////////////////////////////////////////////////////
    /////   gauss point loop
    LINALG::Matrix<3,1> spatial_coord_gp;
    ///////////////////////////////////////////////////////////////////
    for (int gp=0; gp<intpoints.IP().nquad; gp++)
    {
      std::vector<double> bdryxi(globdim-1);
      bdryxi[0] = intpoints.IP().qxg[gp][0];
      bdryxi[1] = intpoints.IP().qxg[gp][1];

      std::vector<double> porosity_gp(1);
      porosity_gp[0]=-1234.0;


      // experimental pre cmbe
      if (during_init)
      {
      int action = 81311;
      IMMERSED::InterpolateToImmersedIntPoint <DRT::Element::hex8,  // source
                                               DRT::Element::quad4> // target
                                                           (backgrddis,
                                                            immerseddis,
                                                            *this,
                                                            bdryxi,
                                                            bdryeledisp,
                                                            action,
                                                            porosity_gp // result
                                                           );
      }
      else
      {
        int action = 813112;
        IMMERSED::InterpolateToImmersedIntPoint <DRT::Element::hex8,  // source
                                                 DRT::Element::quad4> // target
                                                             (scatradis, // experimental pre cmbe // backgrddis
                                                              immerseddis,
                                                              *this,
                                                              bdryxi,
                                                              bdryeledisp,
                                                              action,
                                                              porosity_gp // result
                                                             );
      }

      if(porosity_gp[0]==-1234.0)
        dserror("interpolation of porosity failed");
      // DEBUG output
      //std::cout<<porosity_gp[0]<<std::endl;

      // porosity_new - porosity_old
      double delta_porosity_at_gp = porosity_gp[0]-(porosity_at_gp_old->at(this->Id()))[gp];
      if(abs(delta_porosity_at_gp)<tol)
        delta_porosity_at_gp = 0.0;

      if(abs(immersedmanager->GetDeltaPorosityMax())<abs(delta_porosity_at_gp))
      immersedmanager->SetDeltaPorosityMax(delta_porosity_at_gp);

      // fill map of porosities (will be used to update porosity of old timestep after convergence in global algorithm)
      (porosity_at_gp_->at(this->Id()))[gp] = porosity_gp[0];


        // get shape functions and derivatives in the plane of the element
      DRT::UTILS::shape_function_2D(funct,bdryxi[0],bdryxi[1],Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,bdryxi[0],bdryxi[1],Shape());

      DRT::UTILS::shape_function_3D(parent_funct,parent_xi(gp,0),parent_xi(gp,1),parent_xi(gp,2),this->ParentElement()->Shape());
      DRT::UTILS::shape_function_3D_deriv1(parent_deriv_notrafo,parent_xi(gp,0),parent_xi(gp,1),parent_xi(gp,2),this->ParentElement()->Shape());

      ////////////////////////////////////////////////////
      // calc unitnormal N in material configuration
      ///////////////////////////////////////////////////
      LINALG::Matrix<3,1> unitnormal;
      std::vector<double> normal(globdim);

      // note that the length of this normal is the area dA

      // compute dXYZ / drs
      LINALG::Matrix<2,3> dxyzdrs;
      dxyzdrs.MultiplyNT(deriv,xrefe); // to calculate unitnormal in current config. argument must be xcurr (xrefe instead)

      normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
      normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
      normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

      LINALG::Matrix<2,2> metrictensor;
      metrictensor.MultiplyNT(dxyzdrs,dxyzdrs);
      double detA = sqrt( metrictensor(0,0)*metrictensor(1,1)-metrictensor(0,1)*metrictensor(1,0) );

      for(int i=0;i<globdim;++i)
        unitnormal(i,0)=normal[i];
      const double norm2 = unitnormal.Norm2();
      unitnormal.Scale(1.0/norm2);

      // get space point from where line is originatin
      spatial_coord_gp.MultiplyNT(xcurr,funct);

      // initialize pore (void) diameter
      double voiddiameter = -1234.0;
      //double voiddiameter_old = -1234.0;
      ///////////////////////////////////////////////////////////////
      //
      // Calculate pore diameter from local porosity
      //
      // Model descriptions:
      //
      // 1) Chaterjee2012:
      //
      //    Assumptions: ECM is isotropic and homogenuous.
      //    R := uniform fiber radius.
      //    phi := porosity (void fraction)
      //    Eta := solid fraction (1 - phi)
      //    alpha = ln(1/(1-Eta)) = ln(1/phi)
      //    <r> := mean pore radius
      //    <r^2> := mean squared pore radius
      //
      //    <r>  = (R/2.0)*(sqrt(pi*alpha)/Eta)*(1.0-Erf(sqrt(alpha))
      //    <r^2>= R^2*alpha*((1.0-eta)/eta)^2*(1.0-exp(alpha)*sqrt(pi*alpha)*(1.0-Erf(sqrt(alpha))))
      //
      ///////////////////////////////////////////////////////////////
      static double R = DRT::Problem::Instance()->CellMigrationParams().sublist("CONFINEMENT MODULE").get<double>("ECM_FIBER_RADIUS");
      double penalty = params.get<double>("penalty",-1234.0);
      if(penalty == -1234.0)
        dserror("could not get penalty parameter from parameter list");
      double Eta = 1.0 - porosity_gp[0];
      double alpha = log(1.0/(1.0-Eta));

      //double Eta_old = 1.0 - (porosity_at_gp_old->at(this->Id()))[gp];
      //double alpha_old = log(1.0/(1.0-Eta_old));

      // 2 x equation (5) in Chaterjee2012
      voiddiameter = (R)*(sqrt(PI*alpha)/Eta)*(1.0-erf(sqrt(alpha)));
      //voiddiameter_old = (R)*(sqrt(PI*alpha_old)/Eta_old)*(1.0-erf(sqrt(alpha_old)));

      double gpweight = intpoints.IP().qwgt[gp];

            /* inverse of the Jacobian matrix:
            **            [ x_,r  y_,r  z_,r ]^-1
            **     J^-1 = [ x_,s  y_,s  z_,s ]
            **            [ x_,t  y_,t  z_,t ]
            */
            // compute derivatives N_XYZ at gp w.r.t. material coordinates
            // by N_XYZ = J^-1 * N_rst
            LINALG::Matrix<3,3> invJ(true);
            LINALG::Matrix<3,8> N_XYZ(true);
            invJ.MultiplyNT(parent_deriv_notrafo,parent_xrefe);
            invJ.Invert();
            N_XYZ.Multiply(invJ,parent_deriv_notrafo);

      // (material) deformation gradient F = d xcurr / d xrefe = xcurr * N_XYZ^T
      LINALG::Matrix<3,3> defgrd_inv;
      // deformation gradient
      defgrd_inv.MultiplyNT(parent_xcurr,N_XYZ);
      // Jacobian determinant
      double J = defgrd_inv.Determinant();
      // invert deformation gradient
      defgrd_inv.Invert();

      LINALG::Matrix<3,1> tempvec;
      tempvec.MultiplyTN(defgrd_inv,unitnormal);

      // searchdirection is showing inwards the discretization (inward normal in !!! current !!! configuration)
      LINALG::Matrix<3,1> searchdirection(tempvec);
      const double norm2_currnormal = searchdirection.Norm2();
      searchdirection.Scale(-1.0/norm2_currnormal);

      double gap = 0.0;
      double celldiameter = -1234.0;
      double eps = 2.0e-06;

#ifdef ECM_INTERACTION_SIMPLE
      double simple_ecm_traction_norm=-1234.0;
      double slope = -1234.0;
      if(during_init)
      {
#endif
      //std::cout<<"PROC "<<backgrddis->Comm().MyPID()<<": Evaluate cell diameter from point ["<<spatial_coord_gp(0)<<" "<<spatial_coord_gp(1)<<" "<<spatial_coord_gp(2)<<"] on ele with GID="<<Id()<<std::endl;
      // get intersection point of cell interface and line: x(xi_gp) - unitnormal * lambda (i.e. 'space point' + 'direction' * 'free parameter')
      celldiameter = IMMERSED::IntersectLineWithBoundary<DRT::Element::quad4, // source
                                          DRT::Element::quad4> // target
                                         (immerseddis,
                                          dispnp,
                                          params,
                                          searchdirection,
                                          spatial_coord_gp,
                                          "IMMERSEDCoupling",
                                          -1);

      // DEBUG safety check
      if(celldiameter==-1234.0)
      {
        std::cout<<"PROC "<<backgrddis->Comm().MyPID()<<": Warning! Evaluation of cell diameter failed!\n"
                                                  "        Pertubate space point by very small value and repeat evaluation."<<std::endl;
        //dserror("evaluation of cell diameter failed");
        spatial_coord_gp(0)+=eps;
        spatial_coord_gp(1)+=eps;
        spatial_coord_gp(2)+=eps;
        // get intersection point of cell interface and line: x(xi_gp) - unitnormal * lambda (i.e. 'space point' + 'direction' * 'free parameter')
        // returns cell diameter perpendicular to spatial_coord_gp
        celldiameter = IMMERSED::IntersectLineWithBoundary<DRT::Element::quad4, // source
                                            DRT::Element::quad4> // target
                                           (immerseddis,
                                            dispnp,
                                            params,
                                            searchdirection,
                                            spatial_coord_gp,
                                            "IMMERSEDCoupling",
                                            -1);
        if(celldiameter==-1234.0)
          std::cout<<"Warning! Found no intersection for evaluation of cell diameter!"<<std::endl;

        celldiameter = voiddiameter;
      }

      // calculate gap
      gap = voiddiameter - celldiameter;
      if(abs(gap)<tol)
        gap=0.0;

#ifdef ECM_INTERACTION_SIMPLE
      }
      else if (not during_init)
      {
        gap = 0.0;
        slope = immersedmanager->GetSimpleECMInteractionConstant();
        simple_ecm_traction_norm = slope*(1.0-porosity_gp[0]);
        //std::cout<<"simple_ecm_traction_norm="<<simple_ecm_traction_norm<<std::endl;
      }
#endif

      double multiplicator = -1234.0;
#ifdef VELOCITYDEPENENCY
      double scalarproduct =  1234.0;
      double scalarproduct_disp =  1234.0;
      if(not during_init)
      {
        LINALG::Matrix<3,1> velnp_gp;
        LINALG::Matrix<3,1> dispnp_gp;
        LINALG::Matrix<3,1> unit_inward_normal(tempvec);
        double tempvec_norm = tempvec.Norm2();
        unit_inward_normal.Scale(-1.0/tempvec_norm);

        // get velocity at gp
        // get space point from where line is originating
        velnp_gp.MultiplyNT(velcurr,funct);
        dispnp_gp.MultiplyNT(dispcurr,funct);

        // test if velocity has component in direction of outward normal
        scalarproduct = tempvec.Dot(velnp_gp);
        scalarproduct_disp = unit_inward_normal.Dot(dispnp_gp); // component of displacement in direction inward normal (positive)

        // interface is moving in negative outward normal direction
        if(scalarproduct<-1.0e-06 and abs(delta_porosity_at_gp)>0.0)
        {
          multiplicator = 0.0;
          // effective gap
          gap = ((voiddiameter-voiddiameter_old)/2.0) + scalarproduct_disp; //deltaD_void/2.0 + scalarproduct_disp (zero if change in void space is equal to displacement)
          std::cout<<"Id()="<<Id()<<"  voiddiameter="<<voiddiameter<<"  voiddiameter_old="<<voiddiameter_old<<"  scalarprod="<<scalarproduct_disp<<std::endl;
        }
        else
          multiplicator=1.0;
      }
      else
        multiplicator=1.0;
#else
      multiplicator=1.0;
#endif

          if(multiplicator==-1234.0)
            dserror("determination of variable multiplicator failed");



// defines for penalty ecm force
#define SMOOTHED_PENALTY
#define LOG_LAW
//#define POLYNOMIAL_LAW
//#define LINEAR_LAW

#ifndef SMOOTHED_PENALTY
      // fill element vector
      for(int node = 0; node < nen; node++)
      {
        for (int dof=0;dof<globdim;dof++)
        {
          if (gap<0.0)
            elevector1(node*numdofpernode+dof) += (gpweight*detA)*J*penalty*gap*tempvec(dof)*funct(node);
          else
            elevector1(node*numdofpernode+dof) += 0.0;
        }
      }
#else

      //////////////////////////////////////////////////////////////////////////////////////////////
      // project traction vector at integration point (gp) to current normal at integration point
      //////////////////////////////////////////////////////////////////////////////////////////////
      LINALG::Matrix<3,1> previous_traction_gp(true);
      LINALG::Matrix<3,1> current_traction_gp(true);
      //LINALG::Matrix<3,1> unit_normal(tempvec);
      double norm_tempvec = tempvec.Norm2();
      //unit_normal.Scale(1.0/norm_tempvec);

      previous_traction_gp.Update((penalty_traction_at_gp->at(Id())).at(gp));
      double norm_previous_traction_gp = previous_traction_gp.Norm2();

      current_traction_gp.Update(-norm_previous_traction_gp,tempvec);
      //std::cout<<"norm_previous_traction_gp="<<norm_previous_traction_gp<<std::endl;
      //std::cout<<"norm_tempvec="<<norm_tempvec<<std::endl;

      // now curr_tract_gp contains traction in direct. of new normal
      ///////////////////////////////////////////////////////////////////////////////////////////////

      double gap_max=0.9*celldiameter;

      double delta_penalty_force = -1234.0;

#ifdef LOG_LAW
      // logarithmic law
      //std::cout<<"PROC "<<backgrddis->Comm().MyPID()<<" gap_max = "<<std::setprecision(6)<<gap_max<<"for ele with GID="<<Id()<<std::endl;
      delta_penalty_force = log((1.0/gap_max)*gap+1.0) * penalty; // < 0.0 if celldiameter > voiddiameter
#endif

#ifdef POLYNOMIAL_LAW
      // polynomial (gap^3 + gap) law
      delta_penalty_force = ((1.0/(2.0*gap_max))*(pow(gap,3)) + ((1.0/(gap_max))*gap)) * penalty; // < 0.0 if celldiameter > voiddiameter
#endif

#ifdef LINEAR_LAW
      // linear gap law
      delta_penalty_force = gap * penalty;
#endif

      if(not is_leadingedge_or_rear and (not during_init))
      {
        if     (gap<0.0 and delta_porosity_at_gp<0.0)  // cell > void ; new_void < old_void
          delta_penalty_force+=0.0; // compress
        else if(gap>0.0 and delta_porosity_at_gp<0.0)  // cell < void ; new_void < old_void
          delta_penalty_force=0.0; // relax
        else if(gap<0.0 and delta_porosity_at_gp==0.0) // cell > void ; new_void = old_void
          delta_penalty_force=0.0; // no change in porosity -> already equilibrium force on surface
        else if(gap>0.0 and delta_porosity_at_gp==0.0) // cell < void ; new_void = old_void
          delta_penalty_force=0.0; // no change in porosity -> already equilibrium force on surface
        else if(gap<0.0 and delta_porosity_at_gp>0.0)  // cell > void ; new_void > old_void
          delta_penalty_force+=0.0;
        else if(gap>0.0 and delta_porosity_at_gp>0.0)  // cell < void ; new_void > old_void
          delta_penalty_force=0.0; // relax
      }
      else if(not is_leadingedge_or_rear and during_init) // cell initialization
        delta_penalty_force+=0.0;
      else if(is_leadingedge_or_rear and (not during_init)) // do not load leading edge
      {
#ifdef ECM_INTERACTION_SIMPLE
        simple_ecm_traction_norm = slope*(1.0-0.8);
#endif
        delta_penalty_force=0.0;
      }
      else if (is_leadingedge_or_rear and during_init)
        dserror("ERROR : is_leadingedge_or_rear == true during cell initialization");


      // total current traction at gp may not be pulling traction
#ifdef ECM_INTERACTION_SIMPLE
      if(during_init)
      {
#endif

      if((-norm_previous_traction_gp+delta_penalty_force)<1.0e-12 and not isdummycall) // traction equal or less than zero (good case)
        current_traction_gp.Update(delta_penalty_force,tempvec,1.0);
      else if((-norm_previous_traction_gp+delta_penalty_force)>1.0e-12) // new traction would be pulling traction (bad case)
        current_traction_gp.PutScalar(0.0);    // -> set traction zero

#ifdef ECM_INTERACTION_SIMPLE
      }
      else if(not during_init)
      {
        current_traction_gp.PutScalar(0.0); // safety
        current_traction_gp.Update(-simple_ecm_traction_norm,tempvec);
      }
#endif

      // fill element vector
      for(int node = 0; node < nen; node++)
      {
        for (int dof=0;dof<globdim;dof++)
        {
            elevector1(node*numdofpernode+dof) += (gpweight*detA)*J*current_traction_gp(dof)*funct(node); // compress
        } // fill vector dof-loop

        // fill gap vector to visualize gap over interface in the global algorithm
        if(not is_leadingedge_or_rear)
        {
          elevector2(node*numdofpernode+0) += delta_porosity_at_gp;
          elevector2(node*numdofpernode+1) += (gpweight*detA)*J*gap*funct(node)*multiplicator;
          elevector2(node*numdofpernode+2) += multiplicator;
        }

      } // fill vector node-loop
#endif // SMOOTHED_PENALTY

      /////////////////////////
      //// UPDATE
      /////////////////////////
#ifdef ECM_INTERACTION_SIMPLE
      if(during_init)
      {
#endif

        // update current traction at current gp
        ((penalty_traction_at_gp->at(Id())).at(gp)).Update((1.0/norm_tempvec),current_traction_gp); // curr_tract_gp contains old force in new normal direction + delta force -> total force

#ifdef ECM_INTERACTION_SIMPLE
      }
      else if (not during_init)
        ((penalty_traction_at_gp->at(Id())).at(gp)).PutScalar(0.0); // needed no more after init
#endif


      if(abs(delta_penalty_force)>1.0e-08 and (not is_leadingedge_or_rear))
      {
        // only gaps that contribute to force are reported

        if(immersedmanager->GetGapMax()<gap)
        {
          immersedmanager->SetGapMax(gap);
          immersedmanager->SetMaxGapSpacePoint(spatial_coord_gp);
          LINALG::Matrix<3,1> aux_vec(true);
          aux_vec.Update(1.0,spatial_coord_gp,gap,searchdirection);
          immersedmanager->SetMaxGapSearchDirection(aux_vec);
        }

        if(immersedmanager->GetGapMin()>gap)
        {
          immersedmanager->SetGapMin(gap);
          immersedmanager->SetMinGapSpacePoint(spatial_coord_gp);
          LINALG::Matrix<3,1> aux_vec(true);
          aux_vec.Update(1.0,spatial_coord_gp,gap,searchdirection);
          immersedmanager->SetMinGapSearchDirection(aux_vec);
        }


        if(immersedmanager->GetVoidMax()<voiddiameter)
          immersedmanager->SetVoidMax(voiddiameter);

        if(immersedmanager->GetVoidMin()>voiddiameter)
          immersedmanager->SetVoidMin(voiddiameter);

      }

    } // gp loop
  } // calc_ecm_traction
  break;
  case mark_immersed_elements:
  {

    // obsolete

  }
  break;

  // a surface element routine that calculates the force and stiff contributions of a structural
  // Robin boundary condition (a spring and/or dashpot per unit reference area) - mhv 08/2016
  case calc_struct_robinforcestiff:
  {

    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("displacement");
    Teuchos::RCP<const Epetra_Vector> velonp = discretization.GetState("velocity");

    // time-integration factor for stiffness contribution of dashpot, d(v_{n+1})/d(d_{n+1})
    const double time_fac = params.get("time_fac",0.0);

    const std::vector<int>*    onoff         = params.get<const std::vector<int>* >   ("onoff");
    const std::vector<double>* springstiff   = params.get<const std::vector<double>* >("springstiff");
    const std::vector<double>* dashpotvisc   = params.get<const std::vector<double>* >("dashpotvisc");
    const std::vector<double>* disploffset   = params.get<const std::vector<double>* >("disploffset");

    // type of Robin conditions
    enum RobinType
    {
      none,
      xyz,
      refsurfnormal,
      cursurfnormal
    };

    RobinType rtype = none;

    // get type of Robin condition
    const std::string* direction = params.get<const std::string*>("direction");
    if (*direction == "xyz")
      rtype = xyz;
    else if (*direction == "refsurfnormal")
      rtype = refsurfnormal;
    else if (*direction == "cursurfnormal")
      rtype = cursurfnormal;
    else
      dserror("Unknown type of Robin condition");


    // element geometry update
    const int numdim = 3;
    const int numnode = NumNode();

    const int numdf   = NumDofPerNode(*Nodes()[0]);
    LINALG::SerialDenseMatrix x(numnode,numdim);

    std::vector<double> mydisp(lm.size());
    std::vector<double> myvelo(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydisp,lm);
    DRT::UTILS::ExtractMyValues(*velonp,myvelo,lm);

    // set material configuration
    MaterialConfiguration(x);

    // --------------------------------------------------
    // Now do the nurbs specific stuff
    bool nurbsele=false;

    DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    if(nurbsdis!=NULL) nurbsele=true;

    // knot vectors for parent volume and this surface
    std::vector<Epetra_SerialDenseVector> mypknots(3);
    std::vector<Epetra_SerialDenseVector> myknots (2);

    // NURBS control point weights for all nodes, ie. CPs
    Epetra_SerialDenseVector weights(numnode);

    if(nurbsele)
    {
      // --------------------------------------------------
      // get node weights for nurbs elements
      for (int inode=0; inode<numnode; inode++)
      {
       DRT::NURBS::ControlPoint* cp =
         dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);
       weights(inode) = cp->W();
      }
    }
    // --------------------------------------------------


    std::vector<double> mydisp_refnormal(lm.size());
    std::vector<double> myvelo_refnormal(lm.size());
    Epetra_SerialDenseMatrix N_otimes_N;
    N_otimes_N.Shape(lm.size(),lm.size());

    if (rtype == refsurfnormal)
    {
      std::vector<double> dummy(lm.size()); // dummy vector - we only want the reference normals!
      BuildNormalsAtNodes(elevector2,dummy,true);

      // norm of nodal subvectors of element normal vector
      Epetra_SerialDenseVector norm_refnormal_sq;
      norm_refnormal_sq.Size(numnode);
      for (int node=0; node<numnode; ++node)
        for (int dim=0; dim<numdim; dim++)
          norm_refnormal_sq[node] += elevector2[node*numdf+dim] * elevector2[node*numdf+dim];

      // normalize nodal subvectors of element normal vector
      for (int node=0; node<numnode; ++node)
        for (int dim=0; dim<numdim; dim++)
          elevector2[node*numdf+dim] /= sqrt(norm_refnormal_sq[node]);

      // build nodal N \otimes N matrix
      for (int node=0; node<numnode; ++node)
        for (int dim1=0; dim1<numdf; dim1++)
          for (int dim2=0; dim2<numdf; dim2++)
            N_otimes_N(node*numdf+dim1,node*numdf+dim2) = elevector2[node*numdf+dim1]*elevector2[node*numdf+dim2];

      // (N \otimes N) disp, (N \otimes N) velo
      for (int node=0; node<numnode; ++node)
        for (int dim1=0; dim1<numdim; dim1++)
          for (int dim2=0; dim2<numdim; dim2++)
          {
            mydisp_refnormal[node*numdf+dim1] += N_otimes_N(node*numdf+dim1,node*numdf+dim2) * mydisp[node*numdf+dim2];
            myvelo_refnormal[node*numdf+dim1] += N_otimes_N(node*numdf+dim1,node*numdf+dim2) * myvelo[node*numdf+dim2];
          }

    }

    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);


    /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
    *----------------------------------------------------------------------*/
    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
    for (int gp=0; gp<intpoints.nquad; gp++)
    {
      // set gausspoints from integration rule
      Epetra_SerialDenseVector e(2);
      e(0) = intpoints.qxg[gp][0];
      e(1) = intpoints.qxg[gp][1];

      // get shape functions and derivatives in the plane of the element
      if(!nurbsele)
      {
       DRT::UTILS::shape_function_2D(funct,e(0),e(1),Shape());
       DRT::UTILS::shape_function_2D_deriv1(deriv,e(0),e(1),Shape());
      }
      else
      {
       DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
         (funct,deriv,e,myknots,weights,nurbs9);
      }

      //check for correct input
      for (int checkdof = numdim; checkdof < int(onoff->size()); ++checkdof)
      {
       if ((*onoff)[checkdof] != 0)
         dserror("Number of dimensions in Robin evaluation is 3. Further DoFs are not considered.");
      }

      LINALG::SerialDenseMatrix dxyzdrs(2,3);
      dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);
      LINALG::SerialDenseMatrix  metrictensor(2,2);
      metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      const double detA = sqrt( metrictensor(0,0)*metrictensor(1,1)
                               -metrictensor(0,1)*metrictensor(1,0));

      switch(rtype)
      {
      case xyz:
      {
        for (int dim=0; dim<numdim; dim++)
        {
          if ((*onoff)[dim]) // is this dof activated?
          {
            const double fac_d = intpoints.qwgt[gp] * detA * (*springstiff)[dim];
            const double fac_v = intpoints.qwgt[gp] * detA * (*dashpotvisc)[dim];

            // displacement and velocity at Gauss point
            double dispnp_gp = 0.;
            double velonp_gp = 0.;
            for (int node=0; node<numnode; ++node)
            {
              dispnp_gp += funct[node] * mydisp[node*numdf+dim];
              velonp_gp += funct[node] * myvelo[node*numdf+dim];
            }

            for (int node=0; node<numnode; ++node)
              elevector1[node*numdf+dim] +=
                  funct[node] * (fac_d * (dispnp_gp-(*disploffset)[dim]) + fac_v * velonp_gp);
          }
        }

        for (int dim=0 ; dim<numdim; dim++)
          if ((*onoff)[dim]) // is this dof activated?
          {
            const double fac_d = intpoints.qwgt[gp] * detA * (*springstiff)[dim];
            const double fac_v = intpoints.qwgt[gp] * detA * (*dashpotvisc)[dim];
            for (int node1=0; node1<numnode; ++node1)
              for (int node2=0; node2<numnode; ++node2)
                (elematrix1)(node1*numdf+dim,node2*numdf+dim) +=
                    funct[node1]*funct[node2] * (fac_d + fac_v * time_fac);
          }

      }
      break;

      case refsurfnormal:
      {
        if ((*onoff)[0] != 1) dserror("refsurfnormal Robin condition on 1st dof only!");
          for (int checkdof = 1; checkdof < 3; ++checkdof)
            if ((*onoff)[checkdof] != 0) dserror("refsurfnormal Robin condition on 1st dof only!");

        const double ref_stiff = (*springstiff)[0];
        const double ref_visco = (*dashpotvisc)[0];
        const double ref_disploff = (*disploffset)[0];

        const double fac_d = intpoints.qwgt[gp] * detA * ref_stiff;
        const double fac_v = intpoints.qwgt[gp] * detA * ref_visco;

        for (int dim=0; dim<numdim; dim++)
        {
          // displacement and velocity in normal direction at Gauss point
          double dispnp_refnormal_gp = 0.;
          double velonp_refnormal_gp = 0.;
          for (int node=0; node<numnode; ++node)
          {
            dispnp_refnormal_gp += funct[node] * mydisp_refnormal[node*numdf+dim];
            velonp_refnormal_gp += funct[node] * myvelo_refnormal[node*numdf+dim];
          }

          for (int node=0; node<numnode; ++node)
            elevector1[node*numdf+dim] +=
                funct[node] * (fac_d * (dispnp_refnormal_gp-ref_disploff) + fac_v * velonp_refnormal_gp);
        }

        for (int dim1=0 ; dim1<numdim; dim1++)
          for (int dim2=0 ; dim2<numdim; dim2++)
            for (int node1=0; node1<numnode; ++node1)
              for (int node2=0; node2<numnode; ++node2)
                (elematrix1)(node1*numdf+dim1,node2*numdf+dim2) +=
                    funct[node1]*funct[node2] * (fac_d + fac_v * time_fac) * N_otimes_N(node2*numdf+dim1,node2*numdf+dim2);

      }
      break;

      case cursurfnormal:
      {
        dserror("cursurfnormal option not (yet) implemented in calc_struct_robinforcestiff routine!");
      }
      break;

      default:
        dserror("Unknown type of Robin direction");
        break;
      }
    }
  }
  break;

  default:
    dserror("Unimplemented type of action for StructuralSurface");
    break;
  }
  return 0;
}

/*----------------------------------------------------------------------*
 * Evaluate method for StructuralSurface-Elements               tk 10/07*
 * ---------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralSurface::Evaluate(Teuchos::ParameterList&   params,
                                               DRT::Discretization&      discretization,
                                               LocationArray&            la,
                                               Epetra_SerialDenseMatrix& elematrix1,
                                               Epetra_SerialDenseMatrix& elematrix2,
                                               Epetra_SerialDenseVector& elevector1,
                                               Epetra_SerialDenseVector& elevector2,
                                               Epetra_SerialDenseVector& elevector3)
{
  if (la.Size()==1)
  {
    return Evaluate(
      params,
      discretization,
      la[0].lm_, // location vector is build by the first column of la
      elematrix1,
      elematrix2,
      elevector1,
      elevector2,
      elevector3
      );
  }

  // start with "none"
  DRT::ELEMENTS::StructuralSurface::ActionType act = StructuralSurface::none;

  // get the required action
  std::string action = params.get<std::string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_area_poro")        act = StructuralSurface::calc_struct_area_poro;
  else if (action=="calc_cur_nodal_normals")       act = StructuralSurface::calc_cur_nodal_normals;
  else if (action=="calc_ref_nodal_normals")       act = StructuralSurface::calc_ref_nodal_normals;
  else if (action=="calc_cell_growth")             act = StructuralSurface::calc_cell_growth;
  else
    dserror("Unknown type of action for StructuralSurface");

  // what the element has to do
  switch(act)
  {
  case calc_struct_area_poro:
  {
    CalculateSurfacePorosity(params,discretization,la);
  }
  break;
  case calc_ref_nodal_normals:
  case calc_cur_nodal_normals:
  {
    Evaluate(
          params,
          discretization,
          la[0].lm_, // location vector is build by the first column of la
          elematrix1,
          elematrix2,
          elevector1,
          elevector2,
          elevector3
          );
  }
  break;
  case calc_cell_growth:
  {
    // frequently used constants
    double delta = 2.7e-3;   //< actin monomer size [micrometer]
    double beta = 0.61;      //< 35° in rad --> half branching angle, for calc of growth value to gp normal
    //double k_bT = 0.004114;//< unused

    //////////////////////////////////////////////////////////////
    // get parameters, pointers, and perform safety checks
    /////////////////////////////////////////////////////////////
    DRT::Problem* globalproblem = DRT::Problem::Instance();    //< the problem instance
    int numNode = this->NumNode();                             //< number of nodes of FaceElement
    int numdofpernode = this->NumDofPerNode(*this->Nodes()[0]);//< numdof of FaceElemet nodes

    // get displacement state from structure discretization
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState("displacement");
    if (dispnp == Teuchos::null)
      dserror("Cannot get displacement vector of structure");

    // get mapping from struct bdry ele id to scatra bdry ele id
    const Teuchos::RCP<std::map<int,int> > structscatraelemap = params.get<Teuchos::RCP<std::map<int,int> > >("structscatraelemap");
    if(structscatraelemap==Teuchos::null)
      dserror("could not get structscatraelemap from parameterlist");

    // get matched node maps of structure and scatra dis.
    Teuchos::RCP<Epetra_Map> structnodemap = params.get<Teuchos::RCP<Epetra_Map> >("masternodemap");
    if(structnodemap == Teuchos::null)
      dserror("could not get structnodemap");
    Teuchos::RCP<Epetra_Map> scatranodemap = params.get<Teuchos::RCP<Epetra_Map> >("slavenodemap");
    if(scatranodemap == Teuchos::null)
      dserror("could not get scatranodemap");

    // get scatra dis
    std::string scatradisname(params.get<std::string>("scatradisname"));
    Teuchos::RCP<DRT::Discretization> scatradis  = globalproblem->GetDis(scatradisname);

    // get condition string
    std::string condstring(params.get<std::string>("condstring"));

    // get time step
    const double dt = params.get<double>("dt");

    // get theta params from input file
    double theta = -1.0;
    theta = globalproblem->StructuralDynamicParams().sublist("ONESTEPTHETA").get<double>("THETA");
    if (theta != 1.0)
      dserror("wrong theta!");

    // get pointer to immersedmanager
    DRT::ImmersedFieldExchangeManager* immersedmanager = DRT::ImmersedFieldExchangeManager::Instance();
    if(immersedmanager==NULL)
      dserror("failed to get immersed data manager");

    // dof number of pointed end (of actin filaments) concentration
    static int numofcaps = globalproblem->CellMigrationParams().sublist("PROTRUSION MODULE").get<int>("NUMDOF_PE");
    // dof number of actin monomer concentration
    static int numofcams = globalproblem->CellMigrationParams().sublist("PROTRUSION MODULE").get<int>("NUMDOF_ACTIN");
    if(numofcaps == -1 or numofcams==-1)
      dserror("provide the dof numbers for pointed ends and actin monomers in section ---CELL DYNAMIC/PROTRUSION MODULE !");


    //////////////////////////////////////////////////////////////
    // prerequisites for evaluation at integration point
    /////////////////////////////////////////////////////////////

    // get parent element (Hex8)
    DRT::Element* parentele = this->ParentElement();
    static const int numnodeparentele = parentele->NumNode();
    if (numnodeparentele != 8)
      dserror("only tested for hex8 elements.\n"
              "check this implementation before using other element types !");

    DRT::Element::LocationArray  parentele_la(discretization.NumDofSets());
    parentele->LocationVector(discretization, parentele_la, false);

    LINALG::Matrix<3,8> myeledispnp;
    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,8> >(*dispnp, myeledispnp, parentele_la[0].lm_);

    // get displacement field of face element
    DRT::Node** elenodes = this->Nodes();
    LINALG::Matrix<3,4> eledispnp;
    for (int i=0; i<numNode; i++)
    {
      int elenodeID = elenodes[i]->Id();          // node ID face element node
      for (int k=0; k<numnodeparentele; k++)
      {
        int nodeID = parentele->Nodes()[k]->Id(); // node ID parent element node
        if (elenodeID == nodeID)
          {
              eledispnp(0,i) = myeledispnp(0,k);
              eledispnp(1,i) = myeledispnp(1,k);
              eledispnp(2,i) = myeledispnp(2,k);
          }
      }
    }

    // integration points and weights for boundary (!) gp --> quad4
    const DRT::UTILS::IntPointsAndWeights<2> intpoints (DRT::ELEMENTS::DisTypeToOptGaussRule<DRT::Element::quad4>::rule);

    // get coordinates of gauss point with respect to local parent coordinate system
    LINALG::SerialDenseMatrix pqxg(intpoints.IP().nquad, 3);
    LINALG::Matrix<3,1> pxsi(true);
    LINALG::Matrix<3,3> derivtrafo(true);

    DRT::UTILS::BoundaryGPToParentGP<3>(pqxg,
                                        derivtrafo,
                                        intpoints,
                                        parentele->Shape(),
                                        DRT::Element::quad4,
                                        this->FaceMasterNumber());

    // define reference configuration of coordinates of face element
    LINALG::Matrix<3,4> xrefe;
    LINALG::Matrix<3,4> xcurre;

    // element geometry of face element
    for (int i=0; i<numNode; i++)
    {
      const double* elex = elenodes[i]->X();
        xrefe(0,i) = elex[0];
        xrefe(1,i) = elex[1];
        xrefe(2,i) = elex[2];

        xcurre(0,i) = xrefe(0,i) + eledispnp(0,i);
        xcurre(1,i) = xrefe(1,i) + eledispnp(1,i);
        xcurre(2,i) = xrefe(2,i) + eledispnp(2,i);
    }
    // define current and reference configuration coordinates of parent element
    LINALG::Matrix<3,8> xref;  //reference configuration coordinates
    LINALG::Matrix<3,8> xcurr;  //current configuration coordinates
    DRT::Node** nodes = parentele->Nodes();

    // element geometry of parent element
    for (int i=0; i<numnodeparentele; i++)
    {
      const double* x = nodes[i]->X();
      for ( int k=0; k<3; k++)
      {
        xref(k,i) = x[k];
        xcurr(k,i) = xref(k,i) + myeledispnp(k,i);
      }
    }

    Teuchos::RCP<MAT::So3Material> mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(parentele->Material(0));
    if (mat == Teuchos::null)
      dserror("cast to MAT::So3Material failed");

    // get scatra element matching to ParentElement from matched element map
    int scatraeleID = structscatraelemap->find(this->Id())->second;
    // get scatra element corresponding to structure element
    DRT::Element* scatraele = scatradis->gElement(scatraeleID);
    DRT::Element::LocationArray scatra_la(1);
    scatraele->LocationVector(*scatradis, scatra_la, false);
    const std::vector<int>& scatra_lm = scatra_la[0].lm_;
    // number of scatra element nodes
    const int numofnodes_scatra = scatraele->NumNode();
    // number of scatra element node dofs
    const int numdofpernode_scatra = scatraele->NumDofPerNode(*scatraele->Nodes()[0]);
    int scatranodeID = -1234;

    // get phinp pointer from scatra
    Teuchos::RCP<Epetra_MultiVector> phinp = immersedmanager->GetPointerToPhinps();
    if(phinp == Teuchos::null)
      dserror("phinp = Teuchos::null");
    // get phin pointer from scatra
    Teuchos::RCP<Epetra_MultiVector> phin = immersedmanager->GetPointerToPhins();
    if (phin == Teuchos::null)
      dserror("phin = Teuchos::null");
    // get rates from scatra
    Teuchos::RCP<Epetra_MultiVector> rates = immersedmanager->GetPointerToRates();
    if (rates == Teuchos::null)
      dserror("rates = Teuchos::null");

    /////////////////////////////////////////////////////////////////////////////////
    // loop over all integration points
    ////////////////////////////////////////////////////////////////////////////////
    for (int iquad = 0; iquad<intpoints.IP().nquad; iquad++)
    {
      // coordinates of the current integration point in parent coordinate system --> HEX8
      for(int idim=0; idim<3; idim++)
      {
        pxsi(idim) = pqxg (iquad, idim); //< pxsi := coordinates of current gauss point for parent element (HEX8) [3x1]
      }
      // coordinates of current integration point in face element coordinate system --> QUAD4
      LINALG::Matrix<2,1> xsi(true);
      xsi(0) = intpoints.IP().qxg[iquad][0];
      xsi(1) = intpoints.IP().qxg[iquad][1];

      const int numdim = 3;
      const int numnod = 8;

      // derivatives of parent element shape functions in parent element coordinates system
      LINALG::Matrix<numdim,numnod> pderiv_loc(true);
      DRT::UTILS::shape_function_deriv1<DRT::Element::hex8>(pxsi, pderiv_loc);

      // shapefunct and derivates of face element in face element coordinate system
      LINALG::Matrix<4, 1> shapefunct;
      DRT::UTILS::shape_function<DRT::Element::quad4>(xsi,shapefunct);
      LINALG::Matrix<2,4> deriv(true);
      DRT::UTILS::shape_function_deriv1<DRT::Element::quad4>(xsi, deriv);

      ////////////////////////////////////////////////////////////////////////////////////////
      // calc surface traction on growth surface, influencing polymerisation rate
      ////////////////////////////////////////////////////////////////////////////////////////

      /*get inverse of Jacobian-Determinant
       *             [x_r  y_r    z_r]^-1
       *     J^-1 =  [x_s  y_s    z_s]
       *             [x_t  y_t    z_t]
       */
      // invJ = N_rst*xref^T
      LINALG::Matrix<numdim,numdim> invJ;
      invJ.MultiplyNT(pderiv_loc, xref);

      // compute derivatives N_XYZ at gp with respect to material coordinates by N_XYZ = invJ*N_rst
      LINALG::Matrix<numdim, numnod> N_XYZ;
      N_XYZ.Multiply(invJ, pderiv_loc);

      //material deformation gradient F = dxcurr/dxref = xcurr^T*N_XYZ^T
      LINALG::Matrix<numdim, numdim> defgrd(false);
      defgrd.MultiplyNT(xcurr, N_XYZ);

      // right Cauchy-Green tensor C = F^{T} \cdot F
      LINALG::Matrix<numdim, numdim> cauchygreen;
      cauchygreen.MultiplyTN(defgrd, defgrd);

      // Green-Lagrange strain glstrain = [E11, E22, E33, 2*E12, 2*E23, 2*E31]
//      Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
//      LINALG::Matrix<MAT::NUM_STRESS_3D,1> glstrain(glstrain_epetra.A(), true);
//      // E=1/2*(C-I)
//      glstrain(0) = 0.5 * (cauchygreen(0,0) -1.0);
//      glstrain(1) = 0.5 * (cauchygreen(1,1) -1.0);
//      glstrain(2) = 0.5 * (cauchygreen(2,2) -1.0);
//      glstrain(3) = cauchygreen(0,1);
//      glstrain(4) = cauchygreen(1,2);
//      glstrain(5) = cauchygreen(2,0);

      // evaluate material to get surface stress
//      LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> cmat(true);
//      LINALG::Matrix<MAT::NUM_STRESS_3D, 1> stress(true);
//      mat->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, parentele->Id());

      // stress = PK2
      // convert PK2 to cauchy stress
//      LINALG::Matrix<3,3> cauchystress(true);
//      LINALG::Matrix<3,3> PK2stress;
//      PK2stress(0,0) = stress(0);
//      PK2stress(0,1) = stress(3);
//      PK2stress(0,2) = stress(5);
//      PK2stress(1,0) = PK2stress(0,1);
//      PK2stress(1,1) = stress(1);
//      PK2stress(1,2) = stress(4);
//      PK2stress(2,0) = PK2stress(0,2);
//      PK2stress(2,1) = PK2stress(1,2);
//      PK2stress(2,2) = stress(2);

      // \sigma = \frac{1}{J} * F \cdot PK2 \cdot F^{T}
//      double detF = defgrd.Determinant();
//      LINALG::Matrix<3,3> temp(true);
//      temp.MultiplyNN(defgrd, PK2stress);
//      cauchystress.MultiplyNT(temp, defgrd);
//      cauchystress.Scale(1.0/detF);

      // calc unitnormal N in current configuration
      LINALG::Matrix<3,1> unitnormal;
      LINALG::Matrix<3,1> normal(3);
      // note that the length of this normal is the area dA
      // compute dXYZ / drs
      LINALG::Matrix<2,3> dxyzdrs;
      dxyzdrs.MultiplyNT(deriv,xcurre); // to calculate unitnormal in current config. argument must be xcurr

      normal(0) = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
      normal(1) = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
      normal(2) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

      for(int i=0;i<3;++i)
        unitnormal(i,0)=normal(i);
      const double norm2 = unitnormal.Norm2();
      unitnormal.Scale(1.0/norm2);

      // calculation of traction at gp t = \sigma \cdot n
//      LINALG::Matrix<3,1> traction(true);
//      traction.Multiply(cauchystress,unitnormal);

      ////////////////////////////////////////////////////////////////////////////////////////
      // calc biochemo-mechano dependent growth at gp
      ////////////////////////////////////////////////////////////////////////////////////////

      // initialize growth value
      double growth = -1234.0;

      // PHI CONCENTRATION FOR CALC
      // phinp concentrations at nodes
      LINALG::Matrix<4,1> capen_nd(true);   //< nodal actin pointed end concentration
      LINALG::Matrix<4,1> camnp_nd(true);   //< nodal actin monomer concentration
      LINALG::Matrix<4,1> rate_am_nd(true); //< nodal actin monomer rate
      LINALG::Matrix<4,1> rate_be_nd(true); //< nodal barbed end rate (equal to nodal rate of change of branching points)

      // concentrations at n+1
      std::vector<double> myphinp(scatra_lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,scatra_lm);

      // concentrations at n
      std::vector<double> myphin(scatra_lm.size());
      DRT::UTILS::ExtractMyValues(*phin,myphin,scatra_lm);

      // rates at n+1
      std::vector<double> myrates(scatra_lm.size());
      DRT::UTILS::ExtractMyValues(*rates,myrates,scatra_lm);

      // loop over all scatra element nodes to get actin monomer and actin pointed end concentration at nodes (HEX8)
      for (int node=0; node<numofnodes_scatra; node++)
      {
        // find corresponding struct - scatra nodes
        scatranodeID = scatraele->NodeIds()[node];
        int lid = scatranodemap->LID(scatranodeID);
        int gid = structnodemap->GID(lid);
        // loop over all structure element nodes (QUAD4)
        for (int nd=0; nd<numNode; nd++)
        {
          // check gid
          int fgid = this->NodeIds()[nd];
          if (fgid == gid)
          {
            const int dof_cams = node*numdofpernode_scatra + numofcams;
            const int dof_caps = node*numdofpernode_scatra + numofcaps;

            capen_nd(nd) = myphin[dof_caps];
            camnp_nd(nd) = myphinp[dof_cams];

            rate_am_nd(nd) = myrates[dof_cams];
            rate_be_nd(nd) = myrates[dof_caps];
          }
        }// end structure element node loop
      }// end scatra element node loop

      // concentration at gp
      double capsn_gp   = 0.0;
      double camsnp_gp  = 0.0;
      double rate_am_gp = 0.0;
      double rate_be_gp = 0.0;

      for (int i=0; i<numNode; i++)
      {
        capsn_gp += shapefunct(i)  * capen_nd(i);
        camsnp_gp += shapefunct(i) * camnp_nd(i);

        rate_am_gp += shapefunct(i) * rate_am_nd(i);
        rate_be_gp += shapefunct(i) * rate_be_nd(i);
      }

      if ((capsn_gp > 0.0) and (camsnp_gp > 0.0))
      {
//        std::cout<<"rate actin at gp "<<iquad<<"      ="<<rate_am_gp<<std::endl;
//        std::cout<<"2 * rate barbed ends at gp "<<iquad<<" ="<<rate_be_gp<<std::endl;
        growth = -1.0*( ( rate_am_gp - (2.0 * rate_be_gp) )*dt*delta*cos(beta) )/(2.0*capsn_gp);
//        std::cout<<"growth at gp "<<iquad<<" ="<<growth<<std::endl;
      }
      else
      {
        std::cout<<"WARNING!!! "<<"capsn_gp="<<capsn_gp<<"   camsnp_gp="<<camsnp_gp<<std::endl;
        growth = 0.0;
      }


      ////////////////////////////////////////////////////////////////////////////////////////
      //  LEAST SQUARES METHOD
      //
      //  Calc the least squares error minimized nodal growth from the
      //  gauss point growth, evaluated in this action.
      //
      //  Here we assemble the least squares system matrix and the right-hand side.
      //  The system is solved globally in ssi_partitioned_2wc_protrusionformation.cpp.
      //
      ///////////////////////////////////////////////////////////////////////////////////////

      // sum over all gp: growth of gp * shapefunction
      // loop over element nodes
      for (int node=0; node<numNode; node++)
      {
        for(int dof=0;dof<numdofpernode;dof++)
        {
          elevector1[node*numdofpernode+dof] += growth*unitnormal(dof)*shapefunct(node);
        } // loop dof

        // fill elematrix
        // sum over all gp: combination of shapefunctions
        for (int i=0; i<numNode; i++)
        {
          int dofrow1 = i*numdofpernode+0;
          int dofrow2 = i*numdofpernode+1;
          int dofrow3 = i*numdofpernode+2;

          int dofcol1 = node*numdofpernode +0;
          int dofcol2 = node*numdofpernode +1;
          int dofcol3 = node*numdofpernode +2;

          elematrix1[dofrow1][dofcol1] += shapefunct(node)*shapefunct(i);
          elematrix1[dofrow1][dofcol2] += 0.0;
          elematrix1[dofrow1][dofcol3] += 0.0;

          elematrix1[dofrow2][dofcol1] += 0.0;
          elematrix1[dofrow2][dofcol2] += shapefunct(node)*shapefunct(i);
          elematrix1[dofrow2][dofcol3] += 0.0;

          elematrix1[dofrow3][dofcol1] += 0.0;
          elematrix1[dofrow3][dofcol2] += 0.0;
          elematrix1[dofrow3][dofcol3] += shapefunct(node)*shapefunct(i);
        }// loop i
      }// loop node

    }//end gp loop

}
  break;
  default:
    dserror("Unimplemented type of action for StructuralSurface");
    break;
  }
  return 0;
}

/*----------------------------------------------------------------------*
 * Compute Volume enclosed by surface.                          tk 10/07*
 * ---------------------------------------------------------------------*/
double DRT::ELEMENTS::StructuralSurface::ComputeConstrVols
(
    const LINALG::SerialDenseMatrix& xc,
    const int numnode
)
{
  double V = 0.0;

  //Volume is calculated by evaluating the integral
  // 1/3*int_A(x dydz + y dxdz + z dxdy)

  // we compute the three volumes separately
  for (int indc = 0; indc < 3; indc++)
  {
    //split current configuration between "ab" and "c"
    // where a!=b!=c and a,b,c are in {x,y,z}
    LINALG::SerialDenseMatrix ab= xc;
    LINALG::SerialDenseVector c (numnode);
    for (int i = 0; i < numnode; i++)
    {
      ab(i,indc) = 0.0; // project by z_i = 0.0
      c(i) = xc(i,indc); // extract z coordinate
    }
    // index of variables a and b
    int inda = (indc+1)%3;
    int indb = (indc+2)%3;

    // get gaussrule
    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
    int ngp = intpoints.nquad;

    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);

    /*----------------------------------------------------------------------*
     |               start loop over integration points                     |
     *----------------------------------------------------------------------*/
    for (int gpid = 0; gpid < ngp; ++gpid)
    {
      const double e0 = intpoints.qxg[gpid][0];
      const double e1 = intpoints.qxg[gpid][1];

      // get shape functions and derivatives of shape functions in the plane of the element
      DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

      double detA;
      // compute "metric tensor" deriv*ab, which is a 2x3 matrix with zero indc'th column
      LINALG::SerialDenseMatrix metrictensor(2,3);
      metrictensor.Multiply('N','N',1.0,deriv,ab,0.0);
      //LINALG::SerialDenseMatrix metrictensor(2,2);
      //metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      detA =  metrictensor(0,inda)*metrictensor(1,indb)-metrictensor(0,indb)*metrictensor(1,inda);
      const double dotprodc = funct.Dot(c);
      // add weighted volume at gausspoint
      V -= dotprodc*detA*intpoints.qwgt[gpid];

    }
  }
  return V/3.0;
}

/*----------------------------------------------------------------------*
 * Compute volume and its first and second derivatives          tk 02/09*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::ComputeVolDeriv
(
  const LINALG::SerialDenseMatrix& xc,
  const int numnode,
  const int ndof,
  double& V,
  Teuchos::RCP<Epetra_SerialDenseVector> Vdiff1,
  Teuchos::RCP<Epetra_SerialDenseMatrix> Vdiff2,
  const int minindex,
  const int maxindex
)
{
  // necessary constants
  const int numdim = 3;
  const double invnumind = 1.0/(maxindex-minindex+1.0);

  // initialize
  V = 0.0;
  Vdiff1->Size(ndof);
  if (Vdiff2!=Teuchos::null) Vdiff2->Shape(ndof, ndof);

  //Volume is calculated by evaluating the integral
  // 1/3*int_A(x dydz + y dxdz + z dxdy)

  // we compute the three volumes separately
  for (int indc = minindex; indc < maxindex+1; indc++)
  {
    //split current configuration between "ab" and "c"
    // where a!=b!=c and a,b,c are in {x,y,z}
    LINALG::SerialDenseMatrix ab= xc;
    LINALG::SerialDenseVector c (numnode);
    for (int i = 0; i < numnode; i++)
    {
      ab(i,indc) = 0.0; // project by z_i = 0.0
      c(i) = xc(i,indc); // extract z coordinate
    }
    // index of variables a and b
    int inda = (indc+1)%3;
    int indb = (indc+2)%3;

    // get gaussrule
    const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
    int ngp = intpoints.nquad;

    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);

    /*----------------------------------------------------------------------*
     |               start loop over integration points                     |
     *----------------------------------------------------------------------*/
    for (int gpid = 0; gpid < ngp; ++gpid)
    {
      const double e0 = intpoints.qxg[gpid][0];
      const double e1 = intpoints.qxg[gpid][1];

      // get shape functions and derivatives of shape functions in the plane of the element
      DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

      // evaluate Jacobi determinant, for projected dA*
      std::vector<double> normal(numdim);
      double detA;
      // compute "metric tensor" deriv*xy, which is a 2x3 matrix with zero 3rd column
      LINALG::SerialDenseMatrix metrictensor(2,numdim);
      metrictensor.Multiply('N','N',1.0,deriv,ab,0.0);
      //metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      detA =  metrictensor(0,inda)*metrictensor(1,indb)-metrictensor(0,indb)*metrictensor(1,inda);
      const double dotprodc = funct.Dot(c);
      // add weighted volume at gausspoint
      V -= dotprodc*detA*intpoints.qwgt[gpid];

      //-------- compute first derivative
      for (int i = 0; i < numnode ; i++)
      {
        (*Vdiff1)[3*i+inda] += invnumind*intpoints.qwgt[gpid]*dotprodc*(deriv(0,i)*metrictensor(1,indb)-metrictensor(0,indb)*deriv(1,i));
        (*Vdiff1)[3*i+indb] += invnumind*intpoints.qwgt[gpid]*dotprodc*(deriv(1,i)*metrictensor(0,inda)-metrictensor(1,inda)*deriv(0,i));
        (*Vdiff1)[3*i+indc] += invnumind*intpoints.qwgt[gpid]*funct[i]*detA;
      }

      //-------- compute second derivative
      if (Vdiff2!=Teuchos::null)
      {
        for (int i = 0; i < numnode ; i++)
        {
          for (int j = 0; j < numnode ; j++)
          {
            //"diagonal" (dV)^2/(dx_i dx_j) = 0, therefore only six entries have to be specified
            (*Vdiff2)(3*i+inda,3*j+indb) += invnumind*intpoints.qwgt[gpid]*dotprodc*(deriv(0,i)*deriv(1,j)-deriv(1,i)*deriv(0,j));
            (*Vdiff2)(3*i+indb,3*j+inda) += invnumind*intpoints.qwgt[gpid]*dotprodc*(deriv(0,j)*deriv(1,i)-deriv(1,j)*deriv(0,i));
            (*Vdiff2)(3*i+inda,3*j+indc) += invnumind*intpoints.qwgt[gpid]*funct[j]*(deriv(0,i)*metrictensor(1,indb)-metrictensor(0,indb)*deriv(1,i));
            (*Vdiff2)(3*i+indc,3*j+inda) += invnumind*intpoints.qwgt[gpid]*funct[i]*(deriv(0,j)*metrictensor(1,indb)-metrictensor(0,indb)*deriv(1,j));
            (*Vdiff2)(3*i+indb,3*j+indc) += invnumind*intpoints.qwgt[gpid]*funct[j]*(deriv(1,i)*metrictensor(0,inda)-metrictensor(1,inda)*deriv(0,i));
            (*Vdiff2)(3*i+indc,3*j+indb) += invnumind*intpoints.qwgt[gpid]*funct[i]*(deriv(1,j)*metrictensor(0,inda)-metrictensor(1,inda)*deriv(0,j));
          }
        }
      }

    }
  }
  V*=invnumind;
  return;
}


/*----------------------------------------------------------------------*
 * Compute surface area and its first and second derivatives    lw 05/08*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::ComputeAreaDeriv(const LINALG::SerialDenseMatrix& x,
                                                        const int numnode,
                                                        const int ndof,
                                                        double& A,
                                                        Teuchos::RCP<Epetra_SerialDenseVector> Adiff,
                                                        Teuchos::RCP<Epetra_SerialDenseMatrix> Adiff2)
{
  // initialization
  A = 0.;
  Adiff->Size(ndof);

  if (Adiff2!=Teuchos::null) Adiff2->Shape(ndof, ndof);

  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

  int ngp = intpoints.nquad;

  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseMatrix  deriv(2,numnode);
  LINALG::SerialDenseMatrix  dxyzdrs(2,3);

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid = 0; gpid < ngp; ++gpid)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get derivatives of shape functions in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    std::vector<double> normal(3);
    double detA;
    SurfaceIntegration(detA,normal,x,deriv);
    A += detA*intpoints.qwgt[gpid];

    LINALG::SerialDenseMatrix ddet(3,ndof,true);
    LINALG::SerialDenseMatrix ddet2(3*ndof,ndof,true);
    LINALG::SerialDenseVector jacobi_deriv(ndof,true);

    dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);

    /*--------------- derivation of minor determiants of the Jacobian
     *----------------------------- with respect to the displacements */
    for (int i=0;i<numnode;++i)
    {
      ddet(0,3*i)   = 0.;
      ddet(0,3*i+1) = deriv(0,i)*dxyzdrs(1,2)-deriv(1,i)*dxyzdrs(0,2);
      ddet(0,3*i+2) = deriv(1,i)*dxyzdrs(0,1)-deriv(0,i)*dxyzdrs(1,1);

      ddet(1,3*i)   = deriv(1,i)*dxyzdrs(0,2)-deriv(0,i)*dxyzdrs(1,2);
      ddet(1,3*i+1) = 0.;
      ddet(1,3*i+2) = deriv(0,i)*dxyzdrs(1,0)-deriv(1,i)*dxyzdrs(0,0);

      ddet(2,3*i)   = deriv(0,i)*dxyzdrs(1,1)-deriv(1,i)*dxyzdrs(0,1);
      ddet(2,3*i+1) = deriv(1,i)*dxyzdrs(0,0)-deriv(0,i)*dxyzdrs(1,0);
      ddet(2,3*i+2) = 0.;

      jacobi_deriv(i*3)   = 1/detA*(normal[2]*ddet(2,3*i  )+normal[1]*ddet(1,3*i  ));
      jacobi_deriv(i*3+1) = 1/detA*(normal[2]*ddet(2,3*i+1)+normal[0]*ddet(0,3*i+1));
      jacobi_deriv(i*3+2) = 1/detA*(normal[0]*ddet(0,3*i+2)+normal[1]*ddet(1,3*i+2));
    }

    /*--- calculation of first derivatives of current interfacial area
     *----------------------------- with respect to the displacements */
    for (int i=0;i<ndof;++i)
    {
      (*Adiff)[i] += jacobi_deriv(i)*intpoints.qwgt[gpid];
    }

    if (Adiff2!=Teuchos::null)
    {
      /*--------- second derivates of minor determiants of the Jacobian
       *----------------------------- with respect to the displacements */
      for (int n=0;n<numnode;++n)
      {
        for (int o=0;o<numnode;++o)
        {
          ddet2(n*3+1,o*3+2) = deriv(0,n)*deriv(1,o)-deriv(1,n)*deriv(0,o);
          ddet2(n*3+2,o*3+1) = - ddet2(n*3+1,o*3+2);

          ddet2(ndof+n*3  ,o*3+2) = deriv(1,n)*deriv(0,o)-deriv(0,n)*deriv(1,o);
          ddet2(ndof+n*3+2,o*3  ) = - ddet2(ndof+n*3,o*3+2);

          ddet2(2*ndof+n*3  ,o*3+1) = ddet2(n*3+1,o*3+2);
          ddet2(2*ndof+n*3+1,o*3  ) = - ddet2(2*ndof+n*3,o*3+1);
        }
      }

      /*- calculation of second derivatives of current interfacial areas
       *----------------------------- with respect to the displacements */
      for (int i=0;i<ndof;++i)
      {
        int var1, var2;

        if (i%3==0)           // displacement in x-direction
        {
          var1 = 1;
          var2 = 2;
        }
        else if ((i-1)%3==0)  // displacement in y-direction
        {
          var1 = 0;
          var2 = 2;
        }
        else if ((i-2)%3==0)  // displacement in z-direction
        {
          var1 = 0;
          var2 = 1;
        }
        else
        {
          dserror("calculation of second derivatives of interfacial area failed");
          exit(1);
        }

        for (int j=0;j<ndof;++j)
        {
          (*Adiff2)(i,j) += (-1/detA*jacobi_deriv(j)*jacobi_deriv(i)+1/detA*
                             (ddet(var1,i)*ddet(var1,j)+normal[var1]*ddet2(var1*ndof+i,j)+
                              ddet(var2,i)*ddet(var2,j)+normal[var2]*ddet2(var2*ndof+i,j)))*intpoints.qwgt[gpid];
        }
      }
    }
  }

  return;
}


void DRT::ELEMENTS::StructuralSurface::BuildNormalsAtNodes(Epetra_SerialDenseVector& nodenormals,
                                                           std::vector<double> mydisp,
                                                           bool refconfig)
{
  const int numnode = NumNode();
  const int numdim = 3;

  LINALG::SerialDenseMatrix x(numnode,3);
  if (refconfig)
    MaterialConfiguration(x);
  else
  {
    SpatialConfiguration(x,mydisp);
  }

  for (int i=0; i<numnode; ++i)
  {
    LINALG::Matrix<3,1> loc_coor;
    loc_coor = DRT::UTILS::getNodeCoordinates(i,Shape());

    const double e0 = loc_coor(0);
    const double e1 = loc_coor(1);

    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    double detA;
    std::vector<double> normal(3);
    SurfaceIntegration(detA,normal,x,deriv);

    for (int j=0; j<numdim; ++j)
    {
      nodenormals(numdim*i+j) = normal[j];
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralSurface::CalculateSurfacePorosity(
    Teuchos::ParameterList& params,
    DRT::Discretization&      discretization,
    LocationArray&            la)
{
  // get the parent element
  DRT::Element* parentele = ParentElement();
  const int nenparent = parentele->NumNode();
  // get element location vector and ownerships
  std::vector<int> lmpar;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  parentele->LocationVector(discretization,lmpar,lmowner,lmstride);

  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
  const int ngp = intpoints.nquad;
  Teuchos::RCP<Epetra_SerialDenseVector> poro = Teuchos::rcp(new Epetra_SerialDenseVector(ngp));
  const int numdim = 3;
  const int numnode = NumNode();
  const int noddof = NumDofPerNode(*(Nodes()[0]));

  // element geometry update
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
  std::vector<double> mydisp(lmpar.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lmpar);

  // update element geometry
  Epetra_SerialDenseMatrix xrefe(numdim,nenparent); // material coord. of element
  Epetra_SerialDenseMatrix xcurr(numdim,nenparent); // current  coord. of element

  DRT::Node** nodes = parentele->Nodes();
  for (int i=0; i<nenparent; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(0,i) = x[0];
    xrefe(1,i) = x[1];
    xrefe(2,i) = x[2];

    xcurr(0,i) = xrefe(0,i) + mydisp[i*noddof+0];
    xcurr(1,i) = xrefe(1,i) + mydisp[i*noddof+1];
    xcurr(2,i) = xrefe(2,i) + mydisp[i*noddof+2];
  }

  const int numdofpernode = 4;

  Teuchos::RCP<const Epetra_Vector> velnp = discretization.GetState(1,"fluidvel");
  if (velnp==Teuchos::null) dserror("Cannot get state vector 'fluidvel'");
  // extract local values of the global vectors
  std::vector<double> myvelpres(la[1].lm_.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelpres,la[1].lm_);

  Epetra_SerialDenseVector mypres(numnode);
  for (int inode=0; inode<numnode; ++inode) // number of nodes
  {
    (mypres)(inode,0) = myvelpres[numdim+(inode*numdofpernode)];
  }

  // get coordinates of gauss points w.r.t. local parent coordinate system
  Epetra_SerialDenseMatrix pqxg(intpoints.nquad,3);
  Epetra_SerialDenseMatrix derivtrafo(3,3);

  DRT::UTILS::SurfaceGPToParentGP(pqxg     ,
                                  derivtrafo,
                                  intpoints,
                                  parentele->Shape() ,
                                  Shape()  ,
                                  LSurfNumber());

  Teuchos::RCP< MAT::StructPoro > structmat =
      Teuchos::rcp_dynamic_cast<MAT::StructPoro>(parentele->Material(1));

  for (int gp=0; gp<ngp; ++gp)
  {
    // get shape functions and derivatives in the plane of the element
   // LINALG::SerialDenseVector  funct(nenparent);
    LINALG::SerialDenseMatrix  deriv(3,nenparent);
   // DRT::UTILS::shape_function_3D(funct,pqxg(gp,0),pqxg(gp,1),pqxg(gp,2),parentele->Shape());
    DRT::UTILS::shape_function_3D_deriv1(deriv,pqxg(gp,0),pqxg(gp,1),pqxg(gp,2),parentele->Shape());

    LINALG::SerialDenseVector  funct2D(numnode);
    DRT::UTILS::shape_function_2D(funct2D,intpoints.qxg[gp][0],intpoints.qxg[gp][1],Shape());

    // pressure at integration point
    double press = funct2D.Dot(mypres);

    // get Jacobian matrix and determinant w.r.t. spatial configuration
    //! transposed jacobian "dx/ds"
    LINALG::SerialDenseMatrix xjm(numdim,numdim);
    xjm.Multiply('N','T',1.0,deriv,xcurr,0.0);
    LINALG::SerialDenseMatrix Jmat(numdim,numdim);
    Jmat.Multiply('N','T',1.0,deriv,xrefe,0.0);

    double det = 0.0;
    double detJ = 0.0;

    if (numdim == 3)
    {
      det = xjm(0,0)*(xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2)) + xjm(0,1)*(- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2)) + xjm(0,2)*(xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1));
      detJ = Jmat(0,0)*(Jmat(1,1)*Jmat(2,2) - Jmat(2,1)*Jmat(1,2)) + Jmat(0,1)*(- Jmat(1,0)*Jmat(2,2) + Jmat(2,0)*Jmat(1,2)) + Jmat(0,2)*(Jmat(1,0)*Jmat(2,1) - Jmat(2,0)*Jmat(1,1));
    }
    else dserror("not implemented");

    const double J = det/detJ;

    double porosity=0.0;

    structmat->ComputeSurfPorosity( params,
                                       press,
                                       J,
                                       LSurfNumber(),
                                       gp,
                                       porosity,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       NULL,
                                       true);
  }
}


/*------------------------------------------------------------------------------------------------------------------------*
 * Evaluate the surface integral contribution associated with the J-integral                                      sudhakar 01/15
 * This is computed only when the crack surfaces are subjected to external traction
 *------------------------------------------------------------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::StructuralSurface::Evaluate_J_integral_crack( Teuchos::RCP<DRT::Discretization>& discret,
                                                                                 Teuchos::RCP<Epetra_Vector>& dispvec,
                                                                                 std::map<int,double> supp_func,
                                                                                 LINALG::Matrix<3,1>& normalcrack,
                                                                                 LINALG::Matrix<3,1>& tangentcrack )
{
  std::vector<double> Jvec;
  Jvec.push_back(0.0);
  Jvec.push_back(0.0);

  // get the displacement at each node of parent element
  // This is necessary to compute displacement gradient
  DRT::Element * parele = ParentElement();
  const int pnnode = parele->NumNode();
  Epetra_SerialDenseMatrix pxref(3,pnnode);   // current  coord. of parent
  LINALG::SerialDenseMatrix pardisp(3,pnnode); // displacement of parent
  LINALG::SerialDenseVector  suppele(pnnode);
  {
    const int * pnids = parele->NodeIds();
    DRT::Node** pnodes = parele->Nodes();
    for( int pn = 0; pn < pnnode; pn++ )
    {
      const double* x = pnodes[pn]->X();
      std::vector<double> dn = DRT::CRACK::UTILS::getDisplacementNode( discret, pnids[pn], dispvec );

      pardisp(0,pn) = dn[0];
      pardisp(1,pn) = dn[1];
      pardisp(2,pn) = dn[2];

      pxref(0,pn) = x[0];
      pxref(1,pn) = x[1];
      pxref(2,pn) = x[2];

      suppele(pn) = supp_func[pnids[pn]];
    }
  }

  // get number of nodes and nodeids of surface element
  const int numnode = NumNode();
  DRT::Node ** snodes = Nodes();

  // to get shape functions and its derivatives of surfacelement at Gauss points
  LINALG::SerialDenseVector  xsi(numnode);
  LINALG::SerialDenseMatrix  deriv(2,numnode);

  // to get shape functions and its derivatives of parent element at Gauss points
  LINALG::SerialDenseVector  pxsi(pnnode);
  LINALG::SerialDenseMatrix  pderiv(3,pnnode);

  // get Gauss integration rule
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);
  const int ngp = intpoints.nquad;

  // get nodes of boundary cell in reference configuration
  LINALG::SerialDenseMatrix xref(numnode,3);
  for( int nodno = 0; nodno < numnode; nodno++ )
  {
    const double* x = snodes[nodno]->X();
    xref(nodno,0) = x[0];
    xref(nodno,1) = x[1];
    xref(nodno,2) = x[2];
  }

  Epetra_SerialDenseMatrix gps(ngp,2);  // coordinates of bcell integration points
  for (int gp=0; gp<ngp; gp++)
  {
    gps(gp,0) = intpoints.qxg[gp][0];
    gps(gp,1) = intpoints.qxg[gp][1];
  }
  Epetra_SerialDenseMatrix pqxg(ngp,3); // coordinates of parent element integration points
  // Get the corresponding coordinates in parent element's local coordinate system
  DRT::UTILS::BoundaryGPToParentGP3(pqxg,gps,parele->Shape(),Shape(),LSurfNumber());

  // Compute surface area of this element
  double surfarea = 0.0;
  for (int gp=0; gp<intpoints.nquad; gp++)
  {
    const double e0 = intpoints.qxg[gp][0];
    const double e1 = intpoints.qxg[gp][1];

    // get shape functions and derivatives in the plane of the element
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    std::vector<double> normal(3);
    double detA;
    SurfaceIntegration(detA,normal,xref,deriv);
    const double fac = intpoints.qwgt[gp] * detA;
    surfarea += fac;

  }

  for (int gp=0; gp<ngp; gp++)
  {
    const double e0 = gps(gp,0);
    const double e1 = gps(gp,1);

    // get shape functions and derivatives of surface
    DRT::UTILS::shape_function_2D(xsi,e0,e1,Shape());
    DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

    std::vector<double> norm(3);
    double detA;
    SurfaceIntegration(detA,norm,xref,deriv);
    double norscale = sqrt( norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2] );
    LINALG::SerialDenseVector normal_bc(3);
    normal_bc(0) = norm[0]/norscale;
    normal_bc(1) = norm[1]/norscale;
    normal_bc(2) = norm[2]/norscale;

    double weight = detA * intpoints.qwgt[gp];

    // Find Jacobian and its inverse w.r. to parent ele
    DRT::UTILS::shape_function_3D(pxsi,pqxg(gp,0),pqxg(gp,1),pqxg(gp,2),parele->Shape());
    DRT::UTILS::shape_function_3D_deriv1(pderiv,pqxg(gp,0),pqxg(gp,1),pqxg(gp,2),parele->Shape());
    LINALG::SerialDenseMatrix xjm(3,3);
    LINALG::SerialDenseMatrix xji(3,3);
    xjm.Multiply('N','T',1.0,pderiv,pxref,0.0);

    double det = xjm(0,0)*(xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2)) + xjm(0,1)*(- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2)) + xjm(0,2)*(xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1));
    xji(0,0) = (xjm(1,1)*xjm(2,2)-xjm(2,1)*xjm(1,2))/det;
    xji(1,1) = (xjm(0,0)*xjm(2,2)-xjm(2,0)*xjm(0,2))/det;
    xji(2,2) = (xjm(0,0)*xjm(1,1)-xjm(1,0)*xjm(0,1))/det;
    xji(0,1) = (xjm(2,1)*xjm(0,2)-xjm(0,1)*xjm(2,2))/det;
    xji(0,2) = (xjm(0,1)*xjm(1,2)-xjm(1,1)*xjm(0,2))/det;
    xji(1,0) = (xjm(2,0)*xjm(1,2)-xjm(1,0)*xjm(2,2))/det;
    xji(1,2) = (xjm(1,0)*xjm(0,2)-xjm(0,0)*xjm(1,2))/det;
    xji(2,0) = (xjm(1,0)*xjm(2,1)-xjm(2,0)*xjm(1,1))/det;
    xji(2,1) = (xjm(2,0)*xjm(0,1)-xjm(0,0)*xjm(2,1))/det;

    LINALG::SerialDenseMatrix pderxy(3,pnnode);
    LINALG::SerialDenseMatrix defgrad(3,3);       // deformation gradient
    LINALG::SerialDenseMatrix dispgrad(3,3);      // usual displacement gradient
    pderxy.Multiply('N','N',1.0,xji,pderiv,0.0);
    dispgrad.Multiply('N','T',1.0,pardisp,pderxy,0.0);
    for( int ii=0; ii<3; ii++ )
    {
      for( int jj=0; jj<3; jj++ )
        defgrad(ii,jj) = dispgrad(ii,jj);
    }
    for( int dim = 0; dim < 3; dim++ )
      defgrad(dim,dim) = dispgrad(dim,dim) + 1.0;

    //---
    // Compute support function
    //---
    double supp_func = pxsi.Dot(suppele);
    /*double glo_gaus[3] = {0.0, 0.0, 0.0};
    for( int nn = 0; nn < pnnode; nn++)
    {
      for( int dim = 0; dim < 3; dim++ )
        glo_gaus[dim] += pxref(dim,nn) *  pxsi(nn);
    }
    double tipdist = pow( (glo_gaus[0]-tipcoo_2d(0,0)), 2 ) + pow( (glo_gaus[1]-tipcoo_2d(1,0)), 2 );
    tipdist = sqrt( tipdist );

    double supp_func = 0.0;
    if( tipdist < min_outer_dist )
    {
      supp_func = 1.0 - tipdist / min_outer_dist;
    }*/

    //---
    // Compute Green-Lagrange strain
    //---
    LINALG::SerialDenseMatrix cauchygreen(3,3);
    cauchygreen.Multiply('T','N',1.0,defgrad,defgrad,0.0);
    LINALG::Matrix<6,1> glstrain;
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    //---
    // Compute 2nd Piola-Kirchhoff stress
    //---
    Teuchos::ParameterList params;
    params.set<int>("gp",gp);

    LINALG::Matrix<6,1> stress2pk;
    LINALG::Matrix<6,6> cmat;
    params.set<int>("iostress", INPAR::STR::stress_2pk); // needed for activefiber material; if output is requested only active stresses are written
    LINALG::Matrix<3,3> defgrad_temp(defgrad.A(),true);
    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(ParentElement()->Material());
    so3mat->Evaluate(&defgrad_temp,&glstrain,params,&stress2pk,&cmat,ParentElement()->Id());

    //---
    // Find first Piola Kirchhoff stress tensor
    // P1 = F(P2)
    //---
    LINALG::SerialDenseMatrix first_piola(3,3);
    for( int dim=0; dim<3; dim++ )
    {
      first_piola(dim,0) = defgrad(dim,0)*stress2pk(0,0) + defgrad(dim,1)*stress2pk(3,0) + defgrad(dim,2)*stress2pk(5,0);
      first_piola(dim,1) = defgrad(dim,0)*stress2pk(3,0) + defgrad(dim,1)*stress2pk(1,0) + defgrad(dim,2)*stress2pk(4,0);
      first_piola(dim,2) = defgrad(dim,0)*stress2pk(5,0) + defgrad(dim,1)*stress2pk(4,0) + defgrad(dim,2)*stress2pk(2,0);
    }

    //---
    // Compute traction vector
    //---
    LINALG::SerialDenseVector traction(3);
    traction.Multiply('N','N',1.0,first_piola,normal_bc,0.0);

    //---
    // Compute J-integrand at Gauss point
    //---
    LINALG::SerialDenseVector j_xyz(3);
    j_xyz.Multiply('T','N',1.0,dispgrad,traction,0.0);


    for( int dim=0; dim<3; dim++ )
      j_xyz(dim) = j_xyz(dim) * weight * supp_func;

    for( int dim=0; dim<3; dim++ )
    {
      Jvec[0] += j_xyz(dim) * normalcrack(dim,0);
      Jvec[1] += j_xyz(dim) * tangentcrack(dim,0);
    }

  }
  return Jvec;
}

