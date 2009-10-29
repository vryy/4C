/*!----------------------------------------------------------------------
\file fluid3_surface_evaluate.cpp
\brief

Integrate a Surface Neumann boundary condition on a given boundary
element (tri or quad)

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "fluid3_weak_dbc.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

#include <blitz/array.h>

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3Surface::Evaluate(     ParameterList&            params,
                                                DRT::Discretization&      discretization,
                                                vector<int>&              lm,
                                                Epetra_SerialDenseMatrix& elemat1,
                                                Epetra_SerialDenseMatrix& elemat2,
                                                Epetra_SerialDenseVector& elevec1,
                                                Epetra_SerialDenseVector& elevec2,
                                                Epetra_SerialDenseVector& elevec3)
{
    DRT::ELEMENTS::Fluid3Surface::ActionType act = Fluid3Surface::none;
    string action = params.get<string>("action","none");
    if (action == "none") dserror("No action supplied");
    else if (action == "integrate_Shapefunction")
        act = Fluid3Surface::integrate_Shapefunction;
    else if (action == "area calculation")
        act = Fluid3Surface::areacalc;
    else if (action == "flowrate calculation")
        act = Fluid3Surface::flowratecalc;
    else if (action == "flowrate_deriv")
        act = Fluid3Surface::flowratederiv;
    else if (action == "Outlet impedance")
        act = Fluid3Surface::Outletimpedance;
    else if (action == "calc_node_normal")
        act = Fluid3Surface::calc_node_normal;
    else if (action == "calc_node_curvature")
        act = Fluid3Surface::calc_node_curvature;
    else if (action == "calc_surface_tension")
        act = Fluid3Surface::calc_surface_tension;
    else if (action == "enforce_weak_dbc")
        act = Fluid3Surface::enforce_weak_dbc;
    else if (action == "conservative_outflow_bc")
        act = Fluid3Surface::conservative_outflow_bc;
    else if (action == "calc_Neumann_inflow")
        act = Fluid3Surface::calc_Neumann_inflow;
    else dserror("Unknown type of action for Fluid3_Surface");

    switch(act)
    {
    case integrate_Shapefunction:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      if (parent_->IsAle())
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp!=null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
      }

      IntegrateShapeFunction(params,discretization,lm,elevec1,mydispnp);
      break;
    }
    case areacalc:
    {
        AreaCaculation(params);
        break;
    }
    case flowratecalc:
    {
        FlowRateParameterCalculation(params,discretization,lm);
        break;
    }
    case flowratederiv:
    {
      FlowRateDeriv(params,discretization,lm,elemat1,elemat2,elevec1,elevec2,elevec3);
      break;
    }
    case Outletimpedance:
    {
        ImpedanceIntegration(params,discretization,lm,elevec1);
        break;
    }
    case calc_node_normal:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      if (parent_->IsAle())
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp!=null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
      }
      ElementNodeNormal(params,discretization,lm,elevec1,mydispnp);
      break;
    }
    case calc_node_curvature:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      if (parent_->IsAle())
      {
        dispnp = discretization.GetState("dispnp");
        if (dispnp!=null)
        {
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
        }
      }

      RefCountPtr<const Epetra_Vector> normals;
      vector<double> mynormals;

      normals = discretization.GetState("normals");
      if (normals!=null)
      {
        mynormals.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*normals,mynormals,lm);
      }

      ElementMeanCurvature(params,discretization,lm,elevec1,mydispnp,mynormals);
      break;
    }
    case enforce_weak_dbc:
    {
      return DRT::ELEMENTS::Fluid3SurfaceWeakDBCInterface::Impl(this)->EvaluateWeakDBC(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;
    }
    case conservative_outflow_bc:
    {
      SurfaceConservativeOutflowConsistency(
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;
    }
    case calc_Neumann_inflow:
    {
      NeumannInflow(
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;
    }
    case calc_surface_tension:
    {
      // employs the divergence theorem acc. to Saksono eq. (24) and does not
      // require second derivatives.

      RCP<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }

      vector<double> mynormals;
      vector<double> mycurvature;

      ElementSurfaceTension(params,discretization,lm,elevec1,mydispnp,mynormals,mycurvature);
      break;
    }
    default:
        dserror("Unknown type of action for Fluid3_Surface");
    } // end of switch(act)

    return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3Surface::EvaluateNeumann(
                                           ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // get time-curve factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

  // get time parameter
  const double thsl = params.get("thsl",0.0);

  // get constant density (only relevant for incompressible flow)
  //const double inc_dens = params.get("inc_density",0.0);

  // get flag for low-Mach-number solver
  const bool loma  = params.get<bool>("low-Mach-number solver",false);

  // get discretization type
  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel = this->NumNode();

  // local surface id
  const int surfaceid =lsurface_;

  // Gaussian points
  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4: case nurbs4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9: case nurbs9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints(gaussrule);

  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix  deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix  xyze(3,iel);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  double                    drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  double normalfac=1.0;

  std::vector<Epetra_SerialDenseVector> mypknots(3);
  std::vector<Epetra_SerialDenseVector> myknots (2);

  Epetra_SerialDenseVector weights(iel);

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(Shape()==Fluid3::nurbs4 || Shape()==Fluid3::nurbs9)
  {
    // --------------------------------------------------
    // get knotvector

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

    bool zero_size=knots->GetBoundaryEleAndParentKnots(mypknots     ,
                                                       myknots      ,
                                                       normalfac    ,
                                                       parent_->Id(),
                                                       surfaceid    );

    if(zero_size)
    {
      return 0;
    }

    // --------------------------------------------------
    // get node weights for nurbs elements
    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  // get scalar vector
  RCP<const Epetra_Vector> scanp = discretization.GetState("scanp");
  if (scanp==null) dserror("Cannot get state vector 'scanp'");

  // extract local values from global vector
  vector<double> myscanp(lm.size());
  DRT::UTILS::ExtractMyValues(*scanp,myscanp,lm);

  // create vector for scalar array
  Epetra_SerialDenseVector escanp(iel);

  // insert scalar into element array
  for (int i=0;i<iel;++i)
  {
    escanp(i) = myscanp[3+(i*4)];
  }

  // This is a hack for low-Mach-number flow with temperature
  // equation until material data will be available here
  // get thermodynamic pressure and its time derivative or history
  double thermpress = params.get<double>("thermpress at n+1",0.0);
  double gasconstant = 287.0;

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    if(!(distype == DRT::Element::nurbs9))
    {
      // ------------------------------------------------
      // shape function derivs of boundary element at gausspoint
      DRT::UTILS::shape_function_2D       (funct,e0,e1,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,distype);
    }
    else
    {
      // this is just a temporary work-around
      Epetra_SerialDenseVector gp(2);
      gp(0)=e0;
      gp(1)=e1;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct  ,
         deriv  ,
         gp     ,
         myknots,
         weights,
         distype);
    }

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    // compute temperature and density for low-Mach-number flow
    double dens = 1.0;
    if (loma)
    {
      double temp = 0.0;
      for (int i=0;i<3;++i)
      {
        temp += funct(i)*escanp(i);
      }
      dens = thermpress/(gasconstant*temp);
    }

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor, the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.) and density
    const double fac = intpoints.qwgt[gpid] * drs * curvefac * thsl * dens;

    // factor given by spatial function
    double functfac = 1.0;
    // determine coordinates of current Gauss point
    double coordgp[3];
    coordgp[0]=0.0;
    coordgp[1]=0.0;
    coordgp[2]=0.0;
    for (int i = 0; i< iel; i++)
    {
      coordgp[0] += xyze(0,i) * funct[i];
      coordgp[1] += xyze(1,i) * funct[i];
      coordgp[2] += xyze(2,i) * funct[i];
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        if (func) functnum = (*func)[dim];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgpref,time,NULL);
          }
          else
            functfac = 1.0;
        }
        elevec1[node*numdf+dim]+= funct[node]*(*onoff)[dim]*(*val)[dim]*fac*functfac;
      }
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | apply outflow boundary condition which is necessary for the          |
 | conservative element formulation (since the convective term was      |
 | partially integrated)                                                |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::SurfaceConservativeOutflowConsistency(
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1,
    Epetra_SerialDenseVector&  elevec1)
{
  // ------------------------------------
  //     GET TIME INTEGRATION DATA
  // ------------------------------------
  // we use two timefacs for matrix and right hand side to be able to
  // use the method for both time integrations
  const double timefac_mat = params.get<double>("timefac_mat");
  const double timefac_rhs = params.get<double>("timefac_rhs");

  // ------------------------------------
  //     GET GENERAL ELEMENT DATA
  // ------------------------------------

  // local surface id
  int surfaceid =lsurface_;

  // get distype
  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel = this->NumNode();

  // Gaussian points
  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9: case nurbs9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints(gaussrule);

  // vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix  deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix  xyze(3,iel);

  // the element's normal vector
  Epetra_SerialDenseVector  norm(3);

  // dyadic product of element's normal vector and velocity
  Epetra_SerialDenseMatrix  n_x_u(3,3);

  // velocity at gausspoint
  Epetra_SerialDenseVector  velint(3);

  // 3 temp vector
  Epetra_SerialDenseVector  tempvec(3);

  // 3x3 temp array
  Epetra_SerialDenseMatrix  temp(3,3);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  Epetra_SerialDenseMatrix  dxyzdrs(2,3);
  double                    drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  // ------------------------------------
  // get statevectors from discretisation

  // displacements
  RCP<const Epetra_Vector>      dispnp;
  vector<double>                mydispnp;

  if (parent_->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    for (int i=0;i<iel;++i)
    {
      const int fi=4*i;

      xyze(0,i)+=mydispnp[  fi];
      xyze(1,i)+=mydispnp[1+fi];
      xyze(2,i)+=mydispnp[2+fi];
    }
  }

  // velocities
  RCP<const Epetra_Vector> vel = discretization.GetState("u and p (trial)");
  if (vel==null) dserror("Cannot get state vector 'u and p (trial)'");

  // extract local values from the global vectors
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

  // create object for element array and insert velocity into element array
  Epetra_SerialDenseMatrix  evel(3,iel);

  for (int i=0;i<iel;++i)
  {
    const int fi=4*i;
    evel(0,i) = myvel[  fi];
    evel(1,i) = myvel[1+fi];
    evel(2,i) = myvel[2+fi];
  }

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  double normalfac=0.0;

  std::vector<Epetra_SerialDenseVector> mypknots(3);
  std::vector<Epetra_SerialDenseVector> myknots (2);

  const int piel= parent_->NumNode();

  Epetra_SerialDenseVector weights(iel);
  Epetra_SerialDenseVector pweights(piel);


  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(Shape()==Fluid3::nurbs4 || Shape()==Fluid3::nurbs9)
  {
    // --------------------------------------------------
    // get knotvector

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

    bool zero_size=knots->GetBoundaryEleAndParentKnots(mypknots     ,
                                                       myknots      ,
                                                       normalfac    ,
                                                       parent_->Id(),
                                                       surfaceid    );
    if(zero_size)
    {
      return;
    }

    // --------------------------------------------------
    // get node weights for nurbs elements
    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }

    // extract node coords
    for(int i=0;i<piel;++i)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (parent_->Nodes()[i]);

      pweights(i) = cp->W();
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    if(!(distype == DRT::Element::nurbs9))
    {
      // ------------------------------------------------
      // shape function derivs of boundary element at gausspoint
      DRT::UTILS::shape_function_2D       (funct,e0,e1,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,distype);
    }
    else
    {
      // this is just a temporary work-around
      Epetra_SerialDenseVector gp(2);
      gp(0)=e0;
      gp(1)=e1;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (funct  ,
         deriv  ,
         gp     ,
         myknots,
         weights,
         distype);
    }

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    /*
      |                                              0 1 2
      |                                             +-+-+-+
      |       0 1 2              0...iel-1          | | | | 0
      |      +-+-+-+             +-+-+-+-+          +-+-+-+
      |      | | | | 1           | | | | | 0        | | | | .
      |      +-+-+-+       =     +-+-+-+-+       *  +-+-+-+ .
      |      | | | | 2           | | | | | 1        | | | | .
      |      +-+-+-+             +-+-+-+-+          +-+-+-+
      |                                             | | | | iel-1
      |                                             +-+-+-+
      |
      |       dxyzdrs             deriv              xyze^T
      |
      |
      |                                 +-            -+
      |                                 | dx   dy   dz |
      |                                 | --   --   -- |
      |                                 | dr   dr   dr |
      |     yields           dxyzdrs =  |              |
      |                                 | dx   dy   dz |
      |                                 | --   --   -- |
      |                                 | ds   ds   ds |
      |                                 +-            -+
      |
    */
    dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);

    // values are multiplied by the product from inf. area element,
    // and the gauss weight
    const double fac = intpoints.qwgt[gpid] * drs;


    // compute this element's normal vector scaled by infinitesimal area
    // element and gaussweight


    // ------------------------------------------------
    // compute normal
    if(distype!=DRT::Element::nurbs9)
    {
      double length = 0.0;
      norm(0) = (xyze(1,1)-xyze(1,0))*(xyze(2,2)-xyze(2,0))
        -
        (xyze(2,1)-xyze(2,0))*(xyze(1,2)-xyze(1,0));
      norm(1) = (xyze(2,1)-xyze(2,0))*(xyze(0,2)-xyze(0,0))
        -
        (xyze(0,1)-xyze(0,0))*(xyze(2,2)-xyze(2,0));
      norm(2) = (xyze(0,1)-xyze(0,0))*(xyze(1,2)-xyze(1,0))
        -
        (xyze(1,1)-xyze(1,0))*(xyze(0,2)-xyze(0,0));

      length = norm.Norm2();

      norm(0) = (1.0/length)*norm(0);
      norm(1) = (1.0/length)*norm(1);
      norm(2) = (1.0/length)*norm(2);
    }
    else
    {
      /*
      |
      |                      +-  -+     +-  -+
      |                      | dx |     | dx |
      |                      | -- |     | -- |
      |                      | dr |     | ds |
      |                      |    |     |    |
      |             1.0      | dy |     | dy |
      |    n  =  --------- * | -- |  X  | -- |
      |                      | dr |     | ds |
      |          ||.....||   |    |     |    |
      |                      | dz |     | dz |
      |                      | -- |     | -- |
      |                      | dr |     | ds |
      |                      +-  -+     +-  -+
      |
    */
      norm(0) = dxyzdrs(0,1)*dxyzdrs(1,2)-dxyzdrs(1,1)*dxyzdrs(0,2);
      norm(1) = dxyzdrs(0,2)*dxyzdrs(1,0)-dxyzdrs(1,2)*dxyzdrs(0,0);
      norm(2) = dxyzdrs(0,0)*dxyzdrs(1,1)-dxyzdrs(1,0)*dxyzdrs(0,1);

      const double length = norm.Norm2()*normalfac;

      norm(0) = (1.0/length)*norm(0);
      norm(1) = (1.0/length)*norm(1);
      norm(2) = (1.0/length)*norm(2);
    }

    for(int rr=0;rr<3;++rr)
    {
      norm(rr) *= fac;
    }


    /* interpolate velocities to integration point
    //
    //                 +-----
    //                  \
    //        vel(x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    */
    for(int rr=0;rr<3;++rr)
    {
      velint(rr)=funct(0)*evel(rr,0);
      for(int nn=1;nn<iel;++nn)
      {
        velint(rr)+=funct(nn)*evel(rr,nn);
      }
    }

    // compute normal flux
    const double u_o_n = velint(0)*norm(0)+velint(1)*norm(1)+velint(2)*norm(2);

    // rescaled flux (according to time integration)
    const double timefac_mat_u_o_n = timefac_mat*u_o_n;

    // dyadic product of u and n
    n_x_u(0,0) = timefac_mat*velint(0)*norm(0);
    n_x_u(0,1) = timefac_mat*velint(0)*norm(1);
    n_x_u(0,2) = timefac_mat*velint(0)*norm(2);

    n_x_u(1,0) = timefac_mat*velint(1)*norm(0);
    n_x_u(1,1) = timefac_mat*velint(1)*norm(1);
    n_x_u(1,2) = timefac_mat*velint(1)*norm(2);

    n_x_u(2,0) = timefac_mat*velint(2)*norm(0);
    n_x_u(2,1) = timefac_mat*velint(2)*norm(1);
    n_x_u(2,2) = timefac_mat*velint(2)*norm(2);


    for (int ui=0; ui<iel; ++ui) // loop columns
    {
      const int fui   =4*ui;
      const int fuip  =fui+1;
      const int fuipp =fui+2;

      temp(0,0) = n_x_u(0,0)*funct(ui);
      temp(0,1) = n_x_u(0,1)*funct(ui);
      temp(0,2) = n_x_u(0,2)*funct(ui);

      temp(1,0) = n_x_u(1,0)*funct(ui);
      temp(1,1) = n_x_u(1,1)*funct(ui);
      temp(1,2) = n_x_u(1,2)*funct(ui);

      temp(2,0) = n_x_u(2,0)*funct(ui);
      temp(2,1) = n_x_u(2,1)*funct(ui);
      temp(2,2) = n_x_u(2,2)*funct(ui);

      const double timefac_mat_u_o_n_funct_ui = timefac_mat_u_o_n*funct(ui);

      for (int vi=0; vi<iel; ++vi)  // loop rows
      {
        const int fvi   =4*vi;
        const int fvip  =fvi+1;
        const int fvipp =fvi+2;

        /*


                  /                \
                 |                  |
               + |  Du o n , u o v  |
                 |                  |
                  \                /
        */

        elemat1(fvi  ,fui  ) += temp(0,0)*funct(vi);
        elemat1(fvi  ,fuip ) += temp(0,1)*funct(vi);
        elemat1(fvi  ,fuipp) += temp(0,2)*funct(vi);

        elemat1(fvip ,fui  ) += temp(1,0)*funct(vi);
        elemat1(fvip ,fuip ) += temp(1,1)*funct(vi);
        elemat1(fvip ,fuipp) += temp(1,2)*funct(vi);

        elemat1(fvipp,fui  ) += temp(2,0)*funct(vi);
        elemat1(fvipp,fuip ) += temp(2,1)*funct(vi);
        elemat1(fvipp,fuipp) += temp(2,2)*funct(vi);

        /*


                  /                \
                 |                  |
               + |  u o n , Du o v  |
                 |                  |
                  \                /
        */

        const double timefac_mat_u_o_n_funct_ui_funct_vi
          =
          timefac_mat_u_o_n_funct_ui*funct(vi);

        elemat1(fvi  ,fui  ) += timefac_mat_u_o_n_funct_ui_funct_vi;
        elemat1(fvip ,fuip ) += timefac_mat_u_o_n_funct_ui_funct_vi;
        elemat1(fvipp,fuipp) += timefac_mat_u_o_n_funct_ui_funct_vi;

      } // vi
    } // ui

    tempvec(0)=timefac_rhs*u_o_n*velint(0);
    tempvec(1)=timefac_rhs*u_o_n*velint(1);
    tempvec(2)=timefac_rhs*u_o_n*velint(2);

    for (int ui=0; ui<iel; ++ui) // loop rows  (test functions)
    {
      int fui=4*ui;

      /*


                  /               \
                 |                 |
               + |  u o n , u o v  |
                 |                 |
                  \               /
      */

      elevec1(fui++) -= tempvec(0)*funct(ui);
      elevec1(fui++) -= tempvec(1)*funct(ui);
      elevec1(fui  ) -= tempvec(2)*funct(ui);

    } // ui
  } // end gaussloop

  return;
}// DRT::ELEMENTS::Fluid3Surface::SurfaceConservativeOutflowConsistency


/*----------------------------------------------------------------------*
 | compute potential Neumann inflow                            vg 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::NeumannInflow(
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1,
    Epetra_SerialDenseVector&  elevec1)
{
  //----------------------------------------------------------------------
  // get control parameters for time integration
  //----------------------------------------------------------------------
  // check whether we have a generalized-alpha time-integration scheme
  const bool is_genalpha = params.get<bool>("using generalized-alpha time integration");

  // get timefactor for left hand side
  // One-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = (alpha_F/alpha_M) * gamma * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // get flag for low-Mach-number solver
  const bool loma = params.get<bool>("low-Mach-number solver");

  // get discretization type
  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel = this->NumNode();

  // Gaussian points
  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints(gaussrule);

  // (density-weighted) shape functions and first derivatives
  Epetra_SerialDenseVector funct(iel);
  Epetra_SerialDenseVector densfunct(iel);
  Epetra_SerialDenseMatrix deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix xyze(3,iel);

  // the element's normal vector
  Epetra_SerialDenseVector normal(3);

  // velocity and momentum at gausspoint
  Epetra_SerialDenseVector velint(3);
  Epetra_SerialDenseVector momint(3);

  // metric tensor and area of infinitesimal surface element
  Epetra_SerialDenseMatrix  metrictensor(2,2);
  double                    drs;

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  // determine outward-pointing normal to this element
  normal(0) = (xyze(1,1)-xyze(1,0))*(xyze(2,2)-xyze(2,0))-(xyze(2,1)-xyze(2,0))*(xyze(1,2)-xyze(1,0));
  normal(1) = (xyze(2,1)-xyze(2,0))*(xyze(0,2)-xyze(0,0))-(xyze(0,1)-xyze(0,0))*(xyze(2,2)-xyze(2,0));
  normal(2) = (xyze(0,1)-xyze(0,0))*(xyze(1,2)-xyze(1,0))-(xyze(1,1)-xyze(1,0))*(xyze(0,2)-xyze(0,0));

  // length of normal
  double length = 0.0;
  length = sqrt(normal(0)*normal(0)+normal(1)*normal(1)+normal(2)*normal(2));

  // outward-pointing normal of unit length
  for(int inode=0;inode<3;inode++)
  {
    normal(inode) = normal(inode)/length;
  }

  // get velocity and scalar vector at time n+alpha_F/n+1
  RCP<const Epetra_Vector> velaf = discretization.GetState("velaf");
  RCP<const Epetra_Vector> scaaf = discretization.GetState("scaaf");
  if (velaf==null or scaaf==null)
    dserror("Cannot get state vector 'velaf' and/or 'scaaf'");

  // extract local values from global vector
  vector<double> myvelaf(lm.size());
  vector<double> myscaaf(lm.size());
  DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
  DRT::UTILS::ExtractMyValues(*scaaf,myscaaf,lm);

  // create Epetra objects for scalar array and velocities
  Epetra_SerialDenseMatrix evelaf(3,iel);
  Epetra_SerialDenseVector escaaf(iel);

  // insert velocity and scalar into element array
  for (int i=0;i<iel;++i)
  {
    evelaf(0,i) = myvelaf[0+(i*4)];
    evelaf(1,i) = myvelaf[1+(i*4)];
    evelaf(2,i) = myvelaf[2+(i*4)];

    escaaf(i) = myscaaf[3+(i*4)];
  }

  // This is a hack for low-Mach-number flow with temperature
  // equation until material data will be available here
  // get thermodynamic pressure and its time derivative or history
  double thermpress = params.get<double>("thermpress at n+alpha_F/n+1",0.0);
  double gasconstant = 287.0;

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct,e0,e1,distype);

    // compute velocity and normal velocity
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double normvel = 0.0;
    for (int j=0;j<3;++j)
    {
      velint(j) = 0.0;
      for (int i=0;i<iel;++i)
      {
        velint(j) += funct(i)*evelaf(j,i);
      }
      normvel += velint(j)*normal(j);
    }

    // computation only required for negative normal velocity
    if (normvel<-0.0001)
    {
      // compute measure tensor for surface element, infinitesimal area
      // element drs for the integration and integration factor
      shape_function_2D_deriv1(deriv,e0,e1,distype);
      DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);
      const double fac = intpoints.qwgt[gpid] * drs;

      // compute temperature and density for low-Mach-number flow
      double dens = 1.0;
      if (loma)
      {
        double temp = 0.0;
        for (int i=0;i<3;++i)
        {
          temp += funct(i)*escaaf(i);
        }
        dens = thermpress/(gasconstant*temp);
      }

      // integration factor for left- and right-hand side
      const double lhsfac = dens*normvel*timefac*fac;
      double rhsfac = dens*normvel*fac;
      if (not is_genalpha) rhsfac *= timefac;

      // matrix
      for (int vi=0; vi<iel; ++vi)
      {
        const double vlhs = lhsfac*funct(vi);

        const int fvi   = 4*vi;
        const int fvip  = fvi+1;
        const int fvipp = fvi+2;

        for (int ui=0; ui<iel; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;

          elemat1(fvi  ,fui  ) -= vlhs*funct(ui);
          elemat1(fvip ,fuip ) -= vlhs*funct(ui);
          elemat1(fvipp,fuipp) -= vlhs*funct(ui);
        }
      }

      // rhs
      const double vrhs0 = rhsfac*velint(0);
      const double vrhs1 = rhsfac*velint(1);
      const double vrhs2 = rhsfac*velint(2);
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi   =4*vi;
        const int fvip  =fvi+1;
        const int fvipp =fvi+2;

        elevec1(fvi  ) += funct(vi)*vrhs0;
        elevec1(fvip ) += funct(vi)*vrhs1;
        elevec1(fvipp) += funct(vi)*vrhs2;
      }
    }
  }

  return;
}// DRT::ELEMENTS::Fluid3Surface::NeumannInflow


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)            gjb 07/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::IntegrateShapeFunction(ParameterList& params,
                  DRT::Discretization&       discretization,
                  vector<int>&               lm,
                  Epetra_SerialDenseVector&  elevec1,
                  const std::vector<double>& edispnp)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel   = this->NumNode();

  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }

    // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector      funct       (iel);
  Epetra_SerialDenseMatrix      deriv       (2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix      metrictensor(2,2);
  double                        drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int i=0;i<iel;i++)
    {
      xyze(0,i) += edispnp[4*i];
      xyze(1,i) += edispnp[4*i+1];
      xyze(2,i) += edispnp[4*i+2];
    }
  }

  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  const IntegrationPoints2D  intpoints(gaussrule);

  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration

    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    // values are multiplied by the product from inf. area element,
    // the gauss weight

    const double fac = intpoints.qwgt[gpid] * drs;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        elevec1[node*numdf+dim]+= funct[node] * fac;
      }
    }

  } /* end of loop over integration points gpid */


return;
} // DRT::ELEMENTS::Fluid3Surface::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::ElementNodeNormal(ParameterList& params,
                                                     DRT::Discretization&       discretization,
                                                     vector<int>&               lm,
                                                     Epetra_SerialDenseVector&  elevec1,
                                                     const std::vector<double>& edispnp)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel   = this->NumNode();

  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }

  // allocate vector for shape functions and matrix for derivatives
  Epetra_SerialDenseVector      funct       (iel);
  Epetra_SerialDenseMatrix      deriv       (2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix      xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  Epetra_SerialDenseMatrix 	metrictensor(2,2);
  double                        drs;

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int i=0;i<iel;i++)
    {
      xyze(0,i) += edispnp[numdf*i];
      xyze(1,i) += edispnp[numdf*i+1];
      xyze(2,i) += edispnp[numdf*i+2];
    }
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  const IntegrationPoints2D  intpoints(gaussrule);

  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration

    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    // values are multiplied by the product from inf. area element and
    // the gauss weight

    const double fac = drs * intpoints.qwgt[gpid];

    Epetra_SerialDenseMatrix dxyzdrs(2,3);
    dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);

    // normal vector in gausspoint
    Epetra_SerialDenseVector   norm(3);
    norm[0] = dxyzdrs(0,1)*dxyzdrs(1,2) - dxyzdrs(0,2)*dxyzdrs(1,1);
    norm[1] = dxyzdrs(0,2)*dxyzdrs(1,0) - dxyzdrs(0,0)*dxyzdrs(1,2);
    norm[2] = dxyzdrs(0,0)*dxyzdrs(1,1) - dxyzdrs(0,1)*dxyzdrs(1,0);
    const double length = norm.Norm2();
    norm[0] = norm[0] / length;
    norm[1] = norm[1] / length;
    norm[2] = norm[2] / length;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        elevec1[node*numdf+dim]+= norm[dim] * funct[node] * fac;
      }
      elevec1[node*numdf+3] = 0.0;
    }
  } /* end of loop over integration points gpid */

  return;

} // DRT::ELEMENTS::Fluid3Surface::ElementNodeNormal


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::ElementMeanCurvature(ParameterList& params,
                                                        DRT::Discretization& discretization,
                                                        vector<int>& lm,
                                                        Epetra_SerialDenseVector& elevec1,
                                                        const std::vector<double>& edispnp,
                                                        std::vector<double>& enormals)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  // set number of nodes
  const int iel   = this->NumNode();

  // gaussian points
  const DiscretizationType distype = this->Shape();
  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints(gaussrule);

  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix   deriv(2,iel);

  // the metric tensor and its determinant
  Epetra_SerialDenseMatrix      metrictensor(2,2);
  double sqrtdetg;

  // node coordinates
  Epetra_SerialDenseMatrix 	xyze(3,iel);

  // get node coordinates
  for(int node=0;node<iel;node++)
  {
    xyze(0,node)=this->Nodes()[node]->X()[0];
    xyze(1,node)=this->Nodes()[node]->X()[1];
    xyze(2,node)=this->Nodes()[node]->X()[2];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int node=0;node<iel;node++)
    {
      xyze(0,node) += edispnp[numdf*node     ];
      xyze(1,node) += edispnp[numdf*node + 1 ];
      xyze(2,node) += edispnp[numdf*node + 2 ];
    }
  }

  // node normals
  Epetra_SerialDenseMatrix 	norm_elem(3,iel);

  //set normal vectors to length = 1.0
  for (int node=0;node<iel;++node)
  {
    double length = 0.0;
    for (int dim=0;dim<3;dim++)
    {
      norm_elem(dim,node) = enormals[numdf*node+dim];
      length += norm_elem(dim,node)*norm_elem(dim,node);
    }
    length = sqrt(length);
    for (int dim=0;dim<3;dim++)
    {
      norm_elem(dim,node) = (1.0/length) * norm_elem(dim,node);
    }
  }

  // ============================== loop over nodes ==========================
  for (int node=0;node<iel;node++)
  {
    // get local coordinates of this code
    double this_e0 = 0.0;
    double this_e1 = 0.0;
    switch(distype)
    {
    case quad4: case quad8: case quad9:
      switch(node)
      {
      case 0:
        this_e0 = -1.0;
        this_e1 = -1.0;
        break;
      case 1:
        this_e0 = +1.0;
        this_e1 = -1.0;
        break;
      case 2:
        this_e0 = +1.0;
        this_e1 = +1.0;
        break;
      case 3:
        this_e0 = -1.0;
        this_e1 = +1.0;
        break;
      case 4:
        this_e0 =  0.0;
        this_e1 = -1.0;
        break;
      case 5:
        this_e0 = +1.0;
        this_e1 =  0.0;
        break;
      case 6:
        this_e0 =  0.0;
        this_e1 = +1.0;
        break;
      case 7:
        this_e0 = -1.0;
        this_e1 =  0.0;
        break;
      case 8:
        this_e0 =  0.0;
        this_e1 =  0.0;
        break;
      default:
        dserror("Strange element!\n");
      }
      break;
    case tri3 : case tri6:
      switch(node)
      {
      case 0:
        this_e0 =  0.0;
        this_e1 =  0.0;
        break;
      case 1:
        this_e0 =  1.0;
        this_e1 =  0.0;
        break;
      case 2:
        this_e0 =  0.0;
        this_e1 =  1.0;
        break;
      case 3:
        this_e0 =  0.5;
        this_e1 =  0.0;
        break;
      case 4:
        this_e0 =  0.5;
        this_e1 =  0.5;
        break;
      case 5:
        this_e0 =  0.0;
        this_e1 =  0.5;
        break;
      default:
        dserror("Strange element!\n");
      }
      break;
    default:
      dserror("shape type unknown!\n");
    }

    const double e0 = this_e0;
    const double e1 = this_e1;

    // get shape derivatives at this node
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute surface element's measure tensor and the squareroot of its determinant
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&sqrtdetg);

    // compute mapping
    Epetra_SerialDenseMatrix dxyzdrs(2,3);
    dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);

    // calculate mean curvature H at node.
    double H = 0.0;
    Epetra_SerialDenseMatrix dn123drs(2,3);

    for (int i=0;i<2;i++)
    {
      for (int dim=0;dim<3;dim++)
      {
        for (int node=0;node<iel;node++)
        {
          dn123drs(i,dim) += deriv(i,node) * norm_elem(dim,node);
        }
      }
    }

    //Acc. to Bronstein ..."mittlere Kruemmung":
    double L = 0.0, twoM = 0.0, N = 0.0;
    for (int i=0;i<3;i++)
    {
      L += (-1.0) * dxyzdrs(0,i) * dn123drs(0,i);
      twoM += (-1.0) * dxyzdrs(0,i) * dn123drs(1,i) - dxyzdrs(1,i) * dn123drs(0,i);
      N += (-1.0) * dxyzdrs(1,i) * dn123drs(1,i);
    }
    H = 0.5 *
        (metrictensor(0,0)*N - twoM*metrictensor(0,1) + metrictensor(1,1)*L)
        / (sqrtdetg*sqrtdetg);


    // get the number of elements adjacent to this node. Find out how many
    // will contribute to the interpolated mean curvature value.
    int contr_elements = 0;
    DRT::Node* thisNode = (this->Nodes())[node];
#ifdef DEBUG
    if (thisNode == NULL) dserror("No node!\n");
#endif
    int NumElement = thisNode->NumElement();
    DRT::Element** ElementsPtr = thisNode->Elements();

    // loop over adjacent Fluid3 elements
    for (int ele=0;ele<NumElement;ele++)
    {
      DRT::Element* Element = ElementsPtr[ele];

      // get surfaces
      vector< RCP< DRT::Element > > surfaces = Element->Surfaces();

      // loop over surfaces: how many free surfaces with this node on it?
      for (unsigned int surf=0; surf<surfaces.size(); surf++)
      {
        RCP< DRT::Element > surface = surfaces[surf];
        DRT::Node** NodesPtr = surface->Nodes();
        int numfsnodes = 0;
        bool hasthisnode = false;

        for (int surfnode = 0; surfnode < surface->NumNode(); surfnode++)
        {
          DRT::Node* checkNode = NodesPtr[surfnode];
          // check whether a free surface condition is active on this node
          if (checkNode->GetCondition("FREESURFCoupling") != NULL)
          {
            numfsnodes++;
          }
          if (checkNode->Id() == thisNode->Id())
          {
            hasthisnode = true;
          }
        }

        if (numfsnodes == surface->NumNode() and hasthisnode)
        {
          // this is a free surface adjacent to this node.
          contr_elements++;
        }

      }

    }
#ifdef DEBUG
    if (!contr_elements) dserror("No contributing elements found!\n");
#endif

    for(int dim=0;dim<3;dim++)
    {
      elevec1[node*numdf+dim] = H / contr_elements;
    }
    elevec1[node*numdf+3] = 0.0;
  } // END: loop over nodes

} // DRT::ELEMENTS::Fluid3Surface::ElementMeanCurvature



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::ElementSurfaceTension(ParameterList& params,
                                                         DRT::Discretization& discretization,
                                                         vector<int>& lm,
                                                         Epetra_SerialDenseVector& elevec1,
                                                         const std::vector<double>& edispnp,
                                                         std::vector<double>& enormals,
                                                         std::vector<double>& ecurvature)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  // set number of nodes
  const int iel   = this->NumNode();

  // get time integration parameters
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // isotropic and isothermal surface tension coefficient
  double SFgamma = 0.0;
  // get material data
  RCP<MAT::Material> mat = parent_->Material();
  if (mat==null)
    dserror("no mat from parent!");
  else if (mat->MaterialType()==INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    SFgamma = actmat->Gamma();
  }
  else
    dserror("Newtonian fluid material expected but got type %d", mat->MaterialType());

  // gaussian points
  const DiscretizationType distype = this->Shape();
  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints(gaussrule);

  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix   deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix 	xyze(3,iel);

  // get node coordinates
  for(int node=0;node<iel;node++)
  {
    xyze(0,node)=this->Nodes()[node]->X()[0];
    xyze(1,node)=this->Nodes()[node]->X()[1];
    xyze(2,node)=this->Nodes()[node]->X()[2];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int node=0;node<iel;node++)
    {
      xyze(0,node) += edispnp[numdf*node     ];
      xyze(1,node) += edispnp[numdf*node + 1 ];
      xyze(2,node) += edispnp[numdf*node + 2 ];
    }
  }

  // the metric tensor and its determinant
  Epetra_SerialDenseMatrix      metrictensor(2,2);
  double sqrtdetg;

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/

  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&sqrtdetg);

    // values are multiplied by the product of the determinant of the metric
    // tensor and the gauss weight
    const double fac = intpoints.qwgt[gpid] * sqrtdetg * timefac;

    Epetra_SerialDenseMatrix dxyzdrs(2,3);
    dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);

    double abs_dxyzdr = 0.0;
    double abs_dxyzds = 0.0;
    double pointproduct_rs = 0.0;

    for (int dim=0;dim<3;dim++)
    {
      abs_dxyzdr += dxyzdrs(0,dim) * dxyzdrs(0,dim);
      abs_dxyzds += dxyzdrs(1,dim) * dxyzdrs(1,dim);
      pointproduct_rs += dxyzdrs(0,dim) * dxyzdrs(1,dim);
    }
    abs_dxyzdr = sqrt(abs_dxyzdr);
    abs_dxyzds = sqrt(abs_dxyzds);

    for (int node=0;node<iel;++node)
    {
      for (int dim=0;dim<3;dim++)
      {
        // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
        // function) on Gamma_FS
        // See Saksono eq. (26)
        // discretized as surface gradient * ( Shapefunction-Matrix
        // transformed )

        // This uses a surface_gradient extracted from gauss general
        // formula for 2H...
        // this gives convincing results with TET elements, but HEX
        // elements seem more difficult -> due to edge problems?
        // too many nonlinear iterations
        elevec1[node*numdf+dim] += SFgamma *
                                   (-1.0) / (
                                     sqrtdetg * sqrtdetg  //= abs_dxyzdr * abs_dxyzdr * abs_dxyzds * abs_dxyzds - pointproduct_rs * pointproduct_rs
                                     )
                                   *
                                   (
                                     abs_dxyzds * abs_dxyzds * deriv(0,node) * dxyzdrs(0,dim)
                                     - pointproduct_rs * deriv(0,node) * dxyzdrs(1,dim)
                                     - pointproduct_rs * deriv(1,node) * dxyzdrs(0,dim)
                                     + abs_dxyzdr * abs_dxyzdr * deriv(1,node) * dxyzdrs(1,dim)
                                     )
                                   * fac;

      }
      elevec1[node*numdf+3] = 0.0;
    }
  } /* end of loop over integration points gpid */

} // DRT::ELEMENTS::Fluid3Surface::ElementSurfaceTension

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::AreaCaculation(ParameterList& params)
{
  const int iel   = this->NumNode();
  const DiscretizationType distype = this->Shape();
  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  	funct       (iel);
  LINALG::SerialDenseMatrix  	deriv       (2,iel);

  // node coordinates
  LINALG::SerialDenseMatrix  	xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  LINALG::SerialDenseMatrix 	metrictensor  (2,2);
  double                        drs;

  // get the required material information
  RCP<MAT::Material> mat = parent_->Material();
  double density=0.0, viscosity=0.0;

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    density = actmat->Density();
    viscosity =  actmat->Viscosity();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    density = actmat->Density();
    dserror("How to extract viscosity from Carreau Yasuda material law for artery tree??");
  }
  else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
    density = actmat->Density();
    dserror("How to extract viscosity from modified power law material for artery tree??");
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());

  // set required material data for communication
  params.set<double>("density", density);
  params.set<double>("viscosity", viscosity);

  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
    gaussrule = intrule_quad_4point;
    break;
  case quad8: case quad9:
    gaussrule = intrule_quad_9point;
    break;
  case tri3 :
    gaussrule = intrule_tri_3point;
    break;
  case tri6:
    gaussrule = intrule_tri_6point;
    break;
  default:
    dserror("shape type unknown!\n");
  }

  double area        = params.get<double>("Area calculation");

  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    //Calculate infinitesimal area of element (drs)
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    const double fac = intpoints.qwgt[gpid] * drs;

    area += fac;
  }

  params.set<double>("Area calculation", area);
}//DRT::ELEMENTS::Fluid3Surface::AreaCaculation



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::FlowRateParameterCalculation(ParameterList& params,
                                                                DRT::Discretization&       discretization,
                                                                vector<int>&               lm)
{
  const int iel   = this->NumNode();
  const DiscretizationType distype = this->Shape();
  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  	funct       (iel);
  LINALG::SerialDenseMatrix  	deriv       (2,iel);

  // node coordinates
  LINALG::SerialDenseMatrix  	xyze        (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  LINALG::SerialDenseMatrix 	metrictensor  (2,2);
  double                        drs;

  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
    gaussrule = intrule_quad_4point;
    break;
  case quad8: case quad9:
    gaussrule = intrule_quad_9point;
    break;
  case tri3 :
    gaussrule = intrule_tri_3point;
    break;
  case tri6:
    gaussrule = intrule_tri_6point;
    break;
  default:
    dserror("shape type unknown!\n");
  }

  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");


  // extract local values from the global vectors
  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  double flowrate    = params.get<double>("Outlet flowrate");

  // create blitz objects for element arrays
  const int numnode = NumNode();

  Epetra_SerialDenseVector eprenp(numnode);
  Epetra_SerialDenseMatrix evelnp(3,numnode);
  Epetra_SerialDenseMatrix evhist(3,numnode);

  // split velocity and pressure, insert into element arrays
  for (int i=0;i<numnode;++i)
  {
    evelnp(0,i) = myvelnp[0+(i*4)];
    evelnp(1,i) = myvelnp[1+(i*4)];
    evelnp(2,i) = myvelnp[2+(i*4)];
    eprenp(i)   = myvelnp[3+(i*4)];
  }

  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  // Determine normal to this element
  std::vector<double> dist1(3), dist2(3), normal(3);
  double length;

  for (int i=0; i<3; i++)
  {
    dist1[i] = xyze(i,1)-xyze(i,0);
    dist2[i] = xyze(i,2)-xyze(i,0);
  }

  normal[0] = dist1[1]*dist2[2] - dist1[2]*dist2[1];
  normal[1] = dist1[2]*dist2[0] - dist1[0]*dist2[2];
  normal[2] = dist1[0]*dist2[1] - dist1[1]*dist2[0];

  length = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );

  // here we need an outward normal!!!
  for (int i=0; i<3; i++)
    normal[i] = normal[i] / length;

  const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    //Calculate infinitesimal area of element (drs)
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    const double fac = intpoints.qwgt[gpid] * drs;

    //Compute elment flowrate
    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
	flowrate += funct[node] * evelnp(dim,node)*normal[dim] *fac;
      }
    }
  }
  params.set<double>("Outlet flowrate", flowrate);
}//DRT::ELEMENTS::Fluid3Surface::FlowRateParameterCalculation



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::FlowRateDeriv(ParameterList& params,
                                                 DRT::Discretization&       discretization,
                                                 vector<int>&               lm,
                                                 Epetra_SerialDenseMatrix& elemat1,
                                                 Epetra_SerialDenseMatrix& elemat2,
                                                 Epetra_SerialDenseVector& elevec1,
                                                 Epetra_SerialDenseVector& elevec2,
                                                 Epetra_SerialDenseVector& elevec3)
{
  RefCountPtr<const Epetra_Vector> dispnp;
  vector<double> edispnp;

  if (parent_->IsAle())
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    edispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,edispnp,lm);
  }

  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  // set number of nodes
  const int iel   = this->NumNode();

  // gaussian points
  const DiscretizationType distype = this->Shape();

  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
  case quad4:
      gaussrule = intrule_quad_4point;
      break;
  case quad8: case quad9:
      gaussrule = intrule_quad_9point;
      break;
  case tri3 :
      gaussrule = intrule_tri_3point;
      break;
  case tri6:
      gaussrule = intrule_tri_6point;
      break;
  default:
      dserror("shape type unknown!\n");
  }
  const IntegrationPoints2D  intpoints(gaussrule);

  // order of accuracy of grid velocity determination
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  int gridvel = fdyn.get<int>("GRIDVEL");

  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix   deriv(2,iel);

  // node coordinates
  Epetra_SerialDenseMatrix   xyze(3,iel);

  // get node coordinates
  for(int i=0;i<iel;i++)
  {
    xyze(0,i)=this->Nodes()[i]->X()[0];
    xyze(1,i)=this->Nodes()[i]->X()[1];
    xyze(2,i)=this->Nodes()[i]->X()[2];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int i=0;i<iel;i++)
    {
      xyze(0,i) += edispnp[numdf*i];
      xyze(1,i) += edispnp[numdf*i+1];
      xyze(2,i) += edispnp[numdf*i+2];
    }
  }

  // get nodal velocities and pressures
  RefCountPtr<const Epetra_Vector> convelnp = discretization.GetState("convectivevel");

  if (convelnp==null)
    dserror("Cannot get state vector 'convectivevel'");

  // extract local values from the global vectors
  vector<double> myconvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*convelnp,myconvelnp,lm);

  // extract velocities
  blitz::Array<double, 2> evelnp(3,iel,blitz::ColumnMajorArray<2>());

  for (int i=0;i<iel;++i)
  {
    evelnp(0,i) = myconvelnp[0+(i*4)];
    evelnp(1,i) = myconvelnp[1+(i*4)];
    evelnp(2,i) = myconvelnp[2+(i*4)];
  }


  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // outward gp normal
    Epetra_SerialDenseMatrix dxyzdrs (2,3);
    dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);
    vector<double> normal(3);
    normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);

    //-------------------------------------------------------------------
    //  Q
    const double fac = intpoints.qwgt[gpid];
    blitz::Array<double,1> u(3);
    u = 0.;
    for (int dim=0;dim<3;++dim)
        for (int node=0;node<iel;++node)
          u(dim) += funct[node] * evelnp(dim,node);

    for(int dim=0;dim<3;++dim)
      elevec3[0] += u(dim) * normal[dim] * fac;

    if (params.get<bool>("flowrateonly", false)==false)
    {
      //-------------------------------------------------------------------
      // dQ/du
      for (int node=0;node<iel;++node)
      {
        for (int dim=0;dim<3;++dim)
          elevec1[node*numdf+dim] += funct[node] * normal[dim] * fac;
        elevec1[node*numdf+3] = 0.0;
      }

      //-------------------------------------------------------------------
      // dQ/dd

      // determine derivatives of surface normals wrt mesh displacements
      blitz::Array<double,2> normalderiv(3,iel*3);

      for (int node=0;node<iel;++node)
      {
        normalderiv(0,3*node)   = 0.;
        normalderiv(0,3*node+1) = deriv(0,node)*dxyzdrs(1,2)-deriv(1,node)*dxyzdrs(0,2);
        normalderiv(0,3*node+2) = deriv(1,node)*dxyzdrs(0,1)-deriv(0,node)*dxyzdrs(1,1);

        normalderiv(1,3*node)   = deriv(1,node)*dxyzdrs(0,2)-deriv(0,node)*dxyzdrs(1,2);
        normalderiv(1,3*node+1) = 0.;
        normalderiv(1,3*node+2) = deriv(0,node)*dxyzdrs(1,0)-deriv(1,node)*dxyzdrs(0,0);

        normalderiv(2,3*node)   = deriv(0,node)*dxyzdrs(1,1)-deriv(1,node)*dxyzdrs(0,1);
        normalderiv(2,3*node+1) = deriv(1,node)*dxyzdrs(0,0)-deriv(0,node)*dxyzdrs(1,0);
        normalderiv(2,3*node+2) = 0.;
      }

      for (int node=0;node<iel;++node)
      {
        for (int dim=0;dim<3;++dim)
          for (int iterdim=0;iterdim<3;++iterdim)
            elevec2[node*numdf+dim] += u(iterdim) * normalderiv(iterdim,3*node+dim) * fac;
        elevec2[node*numdf+3] = 0.0;
      }

      // consideration of grid velocity
      if (parent_->IsAle())
      {
        // get time step size
        const double dt = params.get<double>("dt", -1.0);
        if (dt < 0.) dserror("invalid time step size");

        if (gridvel == 1)  // BE time discretization
        {
          for (int node=0;node<iel;++node)
          {
            for (int dim=0;dim<3;++dim)
              elevec2[node*numdf+dim] -= 1.0/dt * funct[node] * normal[dim] * fac;
          }
        }
        else
          dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(du dd)

      for (int unode=0;unode<iel;++unode)
      {
        for (int udim=0;udim<numdf;++udim)
        {
          for (int nnode=0;nnode<iel;++nnode)
          {
            for (int ndim=0;ndim<numdf;++ndim)
            {
              if (udim == 3 or ndim == 3)
                elemat1(unode*numdf+udim,nnode*numdf+ndim) = 0.0;
              else
                elemat1(unode*numdf+udim,nnode*numdf+ndim) = funct[unode] * normalderiv(udim,3*nnode+ndim) * fac;
            }
          }
        }
      }

      //-------------------------------------------------------------------
      // (d^2 Q)/(dd)^2

      // determine second derivatives of surface normals wrt mesh displacements
      blitz::Array<double,3> normalderiv2(3,iel*3,iel*3);
      normalderiv2 = 0.;

      for (int node1=0;node1<iel;++node1)
      {
        for (int node2=0;node2<iel;++node2)
        {
          double temp = deriv(0,node1)*deriv(1,node2)-deriv(1,node1)*deriv(0,node2);

          normalderiv2(0,node1*3+1,node2*3+2) = temp;
          normalderiv2(0,node1*3+2,node2*3+1) = - temp;

          normalderiv2(1,node1*3  ,node2*3+2) = - temp;
          normalderiv2(1,node1*3+2,node2*3  ) = temp;

          normalderiv2(2,node1*3  ,node2*3+1) = temp;
          normalderiv2(2,node1*3+1,node2*3  ) = - temp;
        }
      }

      for (int node1=0;node1<iel;++node1)
      {
        for (int dim1=0;dim1<numdf;++dim1)
        {
          for (int node2=0;node2<iel;++node2)
          {
            for (int dim2=0;dim2<numdf;++dim2)
            {
              if (dim1 == 3 or dim2 == 3)
                elemat2(node1*numdf+dim1,node2*numdf+dim2) = 0.0;
              else
              {
                for (int iterdim=0;iterdim<3;++iterdim)
                  elemat2(node1*numdf+dim1,node2*numdf+dim2) +=
                    u(iterdim) * normalderiv2(iterdim,node1*3+dim1,node2*3+dim2) * fac;
              }
            }
          }
        }
      }

      // consideration of grid velocity
      if (parent_->IsAle())
      {
        // get time step size
        const double dt = params.get<double>("dt", -1.0);
        if (dt < 0.) dserror("invalid time step size");

        if (gridvel == 1)
        {
          for (int node1=0;node1<iel;++node1)
          {
            for (int dim1=0;dim1<3;++dim1)
            {
              for (int node2=0;node2<iel;++node2)
              {
                for (int dim2=0;dim2<3;++dim2)
                {
                  elemat2(node1*numdf+dim1,node2*numdf+dim2) -= (1.0/dt * funct[node1] * normalderiv(dim1, 3*node2+dim2)
                                                                 + 1.0/dt * funct[node2] * normalderiv(dim2, 3*node1+dim1))
                                                                * fac;
                }
              }
            }
          }
        }
        else
          dserror("flowrate calculation: higher order of accuracy of grid velocity not implemented");
      }


      //-------------------------------------------------------------------
    }
  }
}//DRT::ELEMENTS::Fluid3Surface::FlowRateDeriv


 /*----------------------------------------------------------------------*
  |  Impedance related parameters on boundary elements          AC 03/08  |
  *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::ImpedanceIntegration(ParameterList& params,
                  DRT::Discretization&       discretization,
                  vector<int>&               lm,
                  Epetra_SerialDenseVector&  elevec1)
{
  const int iel   = NumNode();
  const DiscretizationType distype = Shape();
  const int numdf = 4;
  const double thsl = params.get("thsl",0.0);

  double invdensity=0.0; // inverse density of my parent element

  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = parent_->Material();

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else
    dserror("Fluid material expected but got type %d", mat->MaterialType());


  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector funct(iel);
  LINALG::SerialDenseMatrix deriv(2,iel);

  // node coordinates
  LINALG::SerialDenseMatrix  x(iel,3);

  // pressure from time integration
  double pressure = params.get<double>("ConvolutedPressure");

  GaussRule2D  gaussrule = intrule2D_undefined;
  switch(distype)
  {
    case quad4:
      gaussrule = intrule_quad_4point;
    break;
    case quad8: case quad9:
      gaussrule = intrule_quad_9point;
    break;
    case tri3 :
      gaussrule = intrule_tri_3point;
    break;
    case tri6:
      gaussrule = intrule_tri_6point;
    break;
    default:
      dserror("shape type unknown!\n");
  }

  for(int i=0; i<iel; ++i)
  {
    x(i,0) = Nodes()[i]->X()[0];
    x(i,1) = Nodes()[i]->X()[1];
    x(i,2) = Nodes()[i]->X()[2];
  }

  const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct,e0,e1,distype);
    shape_function_2D_deriv1(deriv,e0,e1,distype);

    vector<double> normal(3);
    LINALG::SerialDenseMatrix dxyzdrs(2,3);
    dxyzdrs.Multiply('N','N',1.0,deriv,x,0.0);
    normal[0] = dxyzdrs(0,1) * dxyzdrs(1,2) - dxyzdrs(0,2) * dxyzdrs(1,1);
    normal[1] = dxyzdrs(0,2) * dxyzdrs(1,0) - dxyzdrs(0,0) * dxyzdrs(1,2);
    normal[2] = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(0,1) * dxyzdrs(1,0);
    // detA is equal to length of cross product
    const double length = sqrt( normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] );
    // here we need an inward normal!!!
    for (int i=0; i<3; ++i) normal[i] = -normal[i] / length;

    const double fac = intpoints.qwgt[gpid] * length * thsl * pressure * invdensity;
    for (int node=0;node<iel;++node)
      for(int dim=0;dim<3;++dim)
	elevec1[node*numdf+dim] += funct[node] * fac * normal[dim];
  }
  return;
} //DRT::ELEMENTS::Fluid3Surface::ImpedanceIntegration

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
