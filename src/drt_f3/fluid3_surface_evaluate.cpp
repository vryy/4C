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
    else if (action == "enforce_weak_dbc")
        act = Fluid3Surface::enforce_weak_dbc;
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
      xyze(0,i) += edispnp[4*i];
      xyze(1,i) += edispnp[4*i+1];
      xyze(2,i) += edispnp[4*i+2];
    }
  }

  //this element's normal vector
  Epetra_SerialDenseVector   norm(numdf);
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
        elevec1[node*numdf+dim]+=
          funct[node] * fac * norm[dim];
      }
    }
  } /* end of loop over integration points gpid */
} // DRT::ELEMENTS::Fluid3Surface::ElementNodeNormal



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

  if( mat->MaterialType()    != m_carreauyasuda
      && mat->MaterialType() != m_modpowerlaw
      && mat->MaterialType() != m_fluid)
          dserror("Material law is not a fluid");

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_fluid)
  {
    actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();
    density = actmat->m.fluid->density;
    viscosity =  actmat->m.fluid->viscosity;
  }
  else if(mat->MaterialType()== m_carreauyasuda)
  {
    actmat = static_cast<MAT::CarreauYasuda*>(mat.get())->MaterialData();
    density = actmat->m.carreauyasuda->density;
    dserror("How to extract viscosity from Carreau Yasuda material law for artery tree??");
  }
  else if(mat->MaterialType()== m_modpowerlaw)
  {
    actmat = static_cast<MAT::ModPowerLaw*>(mat.get())->MaterialData();
    density = actmat->m.modpowerlaw->density;
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

  if( mat->MaterialType()    != m_carreauyasuda
      && mat->MaterialType() != m_modpowerlaw
      && mat->MaterialType() != m_fluid)
          dserror("Material law is not a fluid");

  MATERIAL* actmat = NULL;

  if(mat->MaterialType()== m_fluid)
  {
    actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();
    invdensity = 1.0/actmat->m.fluid->density;
  }
  else if(mat->MaterialType()== m_carreauyasuda)
  {
    actmat = static_cast<MAT::CarreauYasuda*>(mat.get())->MaterialData();
    invdensity = 1.0/actmat->m.carreauyasuda->density;
  }
  else if(mat->MaterialType()== m_modpowerlaw)
  {
    actmat = static_cast<MAT::ModPowerLaw*>(mat.get())->MaterialData();
    invdensity = 1.0/actmat->m.modpowerlaw->density;
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
