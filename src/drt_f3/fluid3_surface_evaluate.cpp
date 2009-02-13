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
#include "../drt_mat/sutherland_fluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"

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
    else if (action == "Outlet impedance")
        act = Fluid3Surface::Outletimpedance;
    else if (action == "calc_node_normal")
        act = Fluid3Surface::calc_node_normal;
    else if (action == "calc_surface_tension")
        act = Fluid3Surface::calc_surface_tension;
    else if (action == "enforce_weak_dbc")
        act = Fluid3Surface::enforce_weak_dbc;
    else if (action == "conservative_outflow_bc")
        act = Fluid3Surface::conservative_outflow_bc;
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
        FlowRateParameterCaculation(params,discretization,lm);
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
    case calc_surface_tension:
    {

// 2D or TET: To possibilities work. FSTENS1 is a direct implementation of Wall et
// al. eq. (25) with node normals obtained by weighted assembly of element
// normals. Because geometric considerations are used to find the normals of
// our flat (!) surface elements no second derivatives appear. FSTENS2 employs the
// divergence theorem acc. to Saksono eq. (24).

#define FSTENS2
#undef FSTENS1


      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }

      vector<double> mynormals;
#ifdef FSTENS1
      RefCountPtr<const Epetra_Vector> normals;

      normals = discretization.GetState("normals");
      if (normals!=null)
      {
        mynormals.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*normals,mynormals,lm);
      }
#endif

      ElementSurfaceTension(params,discretization,lm,elevec1,mydispnp,mynormals);
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
                                           Epetra_SerialDenseVector& elevec1)
{
  // there are 3 velocities and 1 pressure
  const int numdf = 4;

  // get time parameter
  const double thsl = params.get("thsl",0.0);

  // get constant density (only relevant for incompressible flow)
  //const double inc_dens = params.get("inc_density",0.0);

  // get flag whether outflow stabilization or not
  string outflowstabstr = params.get("outflow stabilization","no_outstab");
  bool outflowstab = false;
  if(outflowstabstr =="yes_outstab") outflowstab = true;

  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel   = this->NumNode();

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

  // get velocity/density vector
  RCP<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  if (vedenp==null) dserror("Cannot get state vector 'vedenp'");

  // extract local values from global vector
  vector<double> myvedenp(lm.size());
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);

  // create blitz object for density array
  blitz::Array<double, 1> edensnp(iel);

  // insert density into element array
  for (int i=0;i<iel;++i)
  {
    edensnp(i) = myvedenp[3+(i*4)];
  }

  // this part will be run when an outflow stabilization term is required
  if (outflowstab)
  {
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

    length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

    // outward pointing normal of length 1
    for (int i=0; i<3; i++) normal[i] = normal[i] / length;

    RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
    if (velnp==null) dserror("Cannot get state vector 'velnp'");

    // extract local values from the global vectors
    vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

    // create blitz object for element array
    blitz::Array<double, 2> evelnp(3,iel,blitz::ColumnMajorArray<2>());

    // insert velocity into element array
    for (int i=0;i<iel;++i)
    {
      evelnp(0,i) = myvelnp[0+(i*4)];
      evelnp(1,i) = myvelnp[1+(i*4)];
      evelnp(2,i) = myvelnp[2+(i*4)];
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

      std::vector<double> vel(3);
      double normvel=0.0;
      for(int dim=0;dim<3;dim++)
      {
        vel[dim] = 0.0;
        for (int node=0;node<iel;++node)
        {
          vel[dim]+= funct[node] * evelnp(dim,node);
        }
        normvel += vel[dim]*normal[dim];
      }

      if (normvel<-0.0001)
      {
        // compute measure tensor for surface element and the infinitesimal
        // area element drs for the integration
        shape_function_2D_deriv1(deriv, e0, e1, distype);
        DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

        const double fac = intpoints.qwgt[gpid] * drs * thsl * normvel;

        for (int node=0;node<iel;++node)
        {
          for(int dim=0;dim<3;dim++)
          {
            elevec1[node*numdf+dim] += funct[node] * edensnp(node) * evelnp(dim,node) * fac;
          }
        }
      }
    } /* end of loop over integration points gpid */
  }
  else
  {
    // find out whether we will use a time curve
    bool usetime = true;
    const double time = params.get("total time",-1.0);
    if (time<0.0) usetime = false;

    // find out whether we will use a time curve and get the factor
    const vector<int>* curve  = condition.Get<vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

    // get values, switches and spatial functions from the condition
    // (assumed to be constant on element boundary)
    const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
    const vector<double>* val   = condition.Get<vector<double> >("val"  );
    const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

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
      DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

      // values are multiplied by the product from inf. area element,
      // the gauss weight, the timecurve factor and the constant
      // belonging to the time integration algorithm (theta*dt for
      // one step theta, 2/3 for bdf with dt const.)
      // Furthermore, there may be a divison by the constant scalar density,
      // only relevant (i.e., it may be unequal 1.0) in incompressible flow case
      const double fac = intpoints.qwgt[gpid] * drs * curvefac * thsl;

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
              functfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dim,coordgpref);
            }
            else
              functfac = 1.0;
          }

          elevec1[node*numdf+dim]+= edensnp(node)*funct[node]*(*onoff)[dim]*(*val)[dim]*fac*functfac;
        }
      }
    } /* end of loop over integration points gpid */
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
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);

    // values are multiplied by the product from inf. area element,
    // and the gauss weight
    const double fac = intpoints.qwgt[gpid] * drs;

    // compute this element's normal vector scaled by infinitesimal area
    // element and gaussweight
    double length = 0.0;
    norm(0) = (xyze(1,1)-xyze(1,0))*(xyze(2,2)-xyze(2,0))-(xyze(2,1)-xyze(2,0))*(xyze(1,2)-xyze(1,0));
    norm(1) = (xyze(2,1)-xyze(2,0))*(xyze(0,2)-xyze(0,0))-(xyze(0,1)-xyze(0,0))*(xyze(2,2)-xyze(2,0));
    norm(2) = (xyze(0,1)-xyze(0,0))*(xyze(1,2)-xyze(1,0))-(xyze(1,1)-xyze(1,0))*(xyze(0,2)-xyze(0,0));

    length = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));

    norm(0) = fac*(1.0/length)*norm(0);
    norm(1) = fac*(1.0/length)*norm(1);
    norm(2) = fac*(1.0/length)*norm(2);

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

    // rescaled flux (accoriding to time integration)
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

  //this element's normal vector
  Epetra_SerialDenseVector   norm(3);
  double length = 0.0;
  norm[0] = (xyze(1,1)-xyze(1,0))*(xyze(2,2)-xyze(2,0))-(xyze(2,1)-xyze(2,0))*(xyze(1,2)-xyze(1,0));
  norm[1] = (xyze(2,1)-xyze(2,0))*(xyze(0,2)-xyze(0,0))-(xyze(0,1)-xyze(0,0))*(xyze(2,2)-xyze(2,0));
  norm[2] = (xyze(0,1)-xyze(0,0))*(xyze(1,2)-xyze(1,0))-(xyze(1,1)-xyze(1,0))*(xyze(0,2)-xyze(0,0));

  length = sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2]);

  norm[0] = (1.0/length)*norm[0];
  norm[1] = (1.0/length)*norm[1];
  norm[2] = (1.0/length)*norm[2];

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

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        elevec1[node*numdf+dim]+= norm[dim] * funct[node] * fac;
      }
      elevec1[node*numdf+3] = 0.0;
    }
  } /* end of loop over integration points gpid */
} // DRT::ELEMENTS::Fluid3Surface::ElementNodeNormal


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::ElementSurfaceTension(ParameterList& params,
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
  else if (mat->MaterialType()==INPAR::MAT::m_sutherland_fluid)
  {
    // no Gamma available
  }
  else
    dserror("newtonian or sutherland fluid material expected but got type %d", mat->MaterialType());

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

#ifdef FSTENS1
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
#endif


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
    const double fac = intpoints.qwgt[gpid] * sqrtdetg;

    Epetra_SerialDenseMatrix dxyzdrs(2,3);
    dxyzdrs.Multiply('N','T',1.0,deriv,xyze,0.0);


#ifdef FSTENS1

    // calculate normal vector at integration point
    Epetra_SerialDenseVector  norm(3);
    for (int dim=0;dim<3;dim++)
    {
      for (int node=0;node<iel;++node)
      {
        norm[dim] += funct[node] * norm_elem(dim,node);
      }
    }
    // set length to 1.0
    double length = 0.0;
    for (int dim=0;dim<3;dim++)
    {
      length += norm[dim] * norm[dim];
    }
    length = sqrt(length);
    for (int dim=0;dim<3;dim++)
    {
      norm[dim] = (1.0/length) * norm[dim];
    }

    // calculate double mean curvature 2*H at integration point.
    double twoH = 0.0;
    Epetra_SerialDenseMatrix dn123drs(2,3);

    for (int i=0;i<2;i++)
    {
      for (int dim=0;dim<3;dim++)
      {
        for (int node=0;node<iel;node++)
        {
          dn123drs(i,dim) = deriv(i,node) * norm_elem(dim,node);
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
    twoH = (metrictensor(0,0)*N - twoM*metrictensor(0,1)
    + metrictensor(1,1)*L)/(sqrtdetg*sqrtdetg);


//     //Acc. to Saksono eq. (4): 2H = - Surface_gradient*norm
//     double abs_dxyzdr = 0.0;
//     double abs_dxyzds = 0.0;
//     double pointproduct = 0.0;

//     for (int dim=0;dim<3;dim++)
//     {
//       abs_dxyzdr += dxyzdrs(0,dim) * dxyzdrs(0,dim);
//       abs_dxyzds += dxyzdrs(1,dim) * dxyzdrs(1,dim);
//       pointproduct += dxyzdrs(0,dim) * dxyzdrs(1,dim);
//     }
//     abs_dxyzdr = sqrt(abs_dxyzdr);
//     abs_dxyzds = sqrt(abs_dxyzds);


//TODO: for 2H = (L+N)/E (as below) 0<E=G (!) and F=0 are required!


//     for(int dim=0;dim<3;dim++)
//     {
//        twoH += dn123drs(1, dim) * dxyzdrs(1, dim)
//         - dn123drs(0, dim) * pointproduct / abs_dxyzdr / abs_dxyzds * dxyzdrs(1, dim)
//         + dn123drs(0, dim) * pointproduct * pointproduct / abs_dxyzdr / abs_dxyzds * dxyzdrs(0, dim)
//         - dn123drs(1, dim) * pointproduct * dxyzdrs(0, dim)
//         + dn123drs(0, dim) * dxyzdrs(0, dim);
//     }
//     twoH = (-1.0) * twoH;


     for (int node=0;node<iel;++node)
     {
       for(int dim=0;dim<3;dim++)
       {
        // according to Saksono (23)
         elevec1[node*numdf+dim]+= SFgamma *
                                   twoH * norm[dim] * funct[node]
                                   * fac;
       }
      elevec1[node*numdf+3] = 0.0;
     }

#else  //----> FSTENS2

    double abs_dxyzdr = 0.0;
    double abs_dxyzds = 0.0;
    double pointproduct = 0.0;

    for (int dim=0;dim<3;dim++)
    {
      abs_dxyzdr += dxyzdrs(0,dim) * dxyzdrs(0,dim);
      abs_dxyzds += dxyzdrs(1,dim) * dxyzdrs(1,dim);
      pointproduct += dxyzdrs(0,dim) * dxyzdrs(1,dim);
    }
    abs_dxyzdr = sqrt(abs_dxyzdr);
    abs_dxyzds = sqrt(abs_dxyzds);

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
        // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
        // function) on Gamma_FS
        // See Saksono eq. (26)
        // discretized as surface gradient * ( Shapefunction-Matrix
        // transformed )

//TODO: This works well. But is E=G required?

        elevec1[node*numdf+dim] += SFgamma *
                                   (-1.0) * (
                                     deriv(1, node) * dxyzdrs(1, dim)
                                     - deriv(0, node) * pointproduct / abs_dxyzdr / abs_dxyzds * dxyzdrs(1, dim)
                                     + deriv(0, node) * pointproduct * pointproduct / abs_dxyzdr / abs_dxyzds * dxyzdrs(0, dim)
                                     - deriv(1, node) * pointproduct * dxyzdrs(0, dim)
                                     + deriv(0, node) * dxyzdrs(0, dim)
                                     )
                                   * fac;
      }
      elevec1[node*numdf+3] = 0.0;
    }

#endif

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
  RefCountPtr<MAT::Material> mat = parent_->Material();
  double density=0.0, viscosity=0.0;

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_sutherland_fluid
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    density = actmat->Density();
    viscosity =  actmat->Viscosity();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_sutherland_fluid)
  {
    //const MAT::SutherlandFluid* actmat = static_cast<const MAT::SutherlandFluid*>(mat.get());
    dserror("How to extract viscosity from Sutherland law material for artery tree??");
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
    dserror("fluid material expected but got type %d", mat->MaterialType());

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
void DRT::ELEMENTS::Fluid3Surface::FlowRateParameterCaculation(ParameterList& params,
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

  RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");

  if (velnp==null)
    dserror("Cannot get state vector 'velnp'");


  // extract local values from the global vectors
  vector<double> myvelnp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

  double flowrate    = params.get<double>("Outlet flowrate");

  // create blitz objects for element arrays
  const int numnode = NumNode();
  blitz::Array<double, 1> eprenp(numnode);
  blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
  blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());

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
}//DRT::ELEMENTS::Fluid3Surface::FlowRateParameterCaculation





  /*----------------------------------------------------------------------*
  |  Impedance related parameters on boundary elements          AC 03/08  |
  *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Surface::ImpedanceIntegration(ParameterList& params,
                  DRT::Discretization&       discretization,
                  vector<int>&               lm,
                  Epetra_SerialDenseVector&  elevec1)
{
  const int iel   = this->NumNode();
  const DiscretizationType distype = this->Shape();
  const int numdf = 4;
  const double thsl = params.get("thsl",0.0);

  double invdensity=0.0; // inverse density of my parent element

  // get material of volume element this surface belongs to
  RefCountPtr<MAT::Material> mat = parent_->Material();

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_sutherland_fluid
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  if(mat->MaterialType()== INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    invdensity = 1.0/actmat->Density();
  }
  else if(mat->MaterialType()== INPAR::MAT::m_sutherland_fluid)
  {
    //const MAT::SutherlandFluid* actmat = static_cast<const MAT::SutherlandFluid*>(mat.get());
    dserror("You really want to compute this with Sutherland fluid???");
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
    dserror("fluid material expected but got type %d", mat->MaterialType());


  // allocate vector for shape functions and matrix for derivatives
  LINALG::SerialDenseVector  	funct       (iel);
  LINALG::SerialDenseMatrix  	deriv       (2,iel);

  // node coordinates
  LINALG::SerialDenseMatrix  	xyze (3,iel);

  // the metric tensor and the area of an infintesimal surface element
  LINALG::SerialDenseMatrix 	metrictensor (2,2);
  double                drs;

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

  // here we need an inward normal!!!
  for (int i=0; i<3; i++)
    normal[i] = -normal[i] / length;

  const IntegrationPoints2D  intpoints(gaussrule);
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e0 = intpoints.qxg[gpid][0];
    const double e1 = intpoints.qxg[gpid][1];

    // get shape functions and derivatives in the plane of the element
    shape_function_2D(funct, e0, e1, distype);
    shape_function_2D_deriv1(deriv, e0, e1, distype);

    // Calculate infinitesimal area of element (drs)
    DRT::UTILS::ComputeMetricTensorForSurface(xyze,deriv,metrictensor,&drs);


    const double fac = intpoints.qwgt[gpid] * drs * thsl * pressure * invdensity;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
	elevec1[node*numdf+dim] += funct[node] * fac * normal[dim];
      }
    }
  }
  return;
}//DRT::ELEMENTS::Fluid3Surface::ImpedanceIntegration

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID3
