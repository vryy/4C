/*!----------------------------------------------------------------------
\file fluid2_line_evaluate.cpp
\brief

<pre>
Maintainer: Peter Gmanitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_nurbs_discret.H"
#include "../drt_lib/drt_control_point.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/modpowerlaw.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid2Line::Evaluate(        ParameterList&            params,
                                                DRT::Discretization&      discretization,
                                                vector<int>&              lm,
                                                Epetra_SerialDenseMatrix& elemat1,
                                                Epetra_SerialDenseMatrix& elemat2,
                                                Epetra_SerialDenseVector& elevec1,
                                                Epetra_SerialDenseVector& elevec2,
                                                Epetra_SerialDenseVector& elevec3)
{
    DRT::ELEMENTS::Fluid2Line::ActionType act = Fluid2Line::none;
    string action = params.get<string>("action","none");

    if (action == "none") dserror("No action supplied");
    else if (action == "integrate_Shapefunction")
        act = Fluid2Line::integrate_Shapefunction;
    else if (action == "calc_node_normal")
        act = Fluid2Line::calc_node_normal;
    else if (action == "calc_surface_tension")
        act = Fluid2Line::calc_surface_tension;
    else dserror("Unknown type of action for Fluid2_Line");

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
    case calc_surface_tension:
    {
      RefCountPtr<const Epetra_Vector> dispnp;
      vector<double> mydispnp;

      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
      RefCountPtr<const Epetra_Vector> normals;
      vector<double> mynormals;

      normals = discretization.GetState("normals");
      if (normals!=null)
      {
        mynormals.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*normals,mynormals,lm);
      }
      ElementSurfaceTension(params,discretization,lm,elemat1,elevec1,mydispnp,mynormals);
      break;
    }
    default:
        dserror("Unknown type of action for Fluid2Line");
    } // end of switch(act)

    return 0;

} // DRT::ELEMENTS::Fluid2Line::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid2Line::EvaluateNeumann(
    ParameterList& params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{

  // there are 2 velocities and 1 pressure
  const int numdf = 3;

  const double thsl = params.get("thsl",0.0);

  // we need the material here to denormalise the body force making sure that 
  // we can cope with physical parameters (force per volume) in the input file
  // check here, if we really have a fluid !!
  // get density
  double invdensity=0.0;
  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = parent_->Material();

  if( mat->MaterialType()    != m_carreauyasuda
      && mat->MaterialType() != m_modpowerlaw
      && mat->MaterialType() != m_fluid)
          dserror("Material law is not a fluid");

  MATERIAL* actmat = NULL;

  // we shall need the inverse of rho
  if(mat->MaterialType()== m_fluid)
  {
    actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();
    invdensity = 1.0 / actmat->m.fluid->density;
  }
  else if(mat->MaterialType()== m_carreauyasuda)
  {
    actmat = static_cast<MAT::CarreauYasuda*>(mat.get())->MaterialData();
    invdensity = 1.0 / actmat->m.carreauyasuda->density;
  }
  else if(mat->MaterialType()== m_modpowerlaw)
  {
    actmat = static_cast<MAT::ModPowerLaw*>(mat.get())->MaterialData();
    invdensity = 1.0 / actmat->m.modpowerlaw->density;
  }
  else
    dserror("fluid material expected but got type %d", mat->MaterialType());

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

  // get values and switches from the condition
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff     = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val       = condition.Get<vector<double> >("val"  );
  const vector<int>*    functions = condition.Get<vector<int> >   ("funct");

    // set number of nodes
  const int iel   = this->NumNode();

  const DiscretizationType distype = this->Shape();

  // gaussian points
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);

  // --------------------------------------------------
  // Now initialise the nurbs specific stuff
  blitz::Array<double,1> myknots;

  blitz::Array<double,1> weights(iel);
  
  // for isogeometric elements
  if(this->Shape()==nurbs2 || this->Shape()==nurbs3)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));
    
    std::vector<blitz::Array<double,1> > surfaceknots(2);
    (*((*nurbsdis).GetKnotVector())).GetEleKnots(surfaceknots,parent_->Id());
    
    //           6    7    8
    //            +---+---+
    //            *       |
    //            *       |
    //            +   +   +
    //           3*   4   |5
    //            *       |
    //            +---+---+
    //           0    1    2  
    //
    // Example: Id()=3 -> (Id())%2=1, i.e. direction v (OK)
    myknots.resize((surfaceknots[(Id())%2]).extent(blitz::firstDim));
    myknots=surfaceknots[(Id())%2];

    DRT::Node** nodes = Nodes();

    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp 
	= 
	dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights(inode) = cp->W();
    }
  }

  // allocate vector for shape functions and for derivatives
  blitz::Array<double,1> funct(iel);
  blitz::Array<double,1> deriv(iel);

  // node coordinates
  blitz::Array<double,2> xye(2,iel);

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }


  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid];
    // get shape functions and derivatives in the line
    if(!(distype == DRT::Element::nurbs2 
	 || 
	 distype == DRT::Element::nurbs3))
    {
      DRT::UTILS::shape_function_1D(funct,
				    e1,
				    distype);
      DRT::UTILS::shape_function_1D_deriv1(deriv,
					   e1,
					   distype);
    }
    else
    {
      DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv(funct  ,
						  deriv  ,
						  e1     ,
						  myknots,
						  weights,
						  distype);
    }

    // The Jacobian is computed using the formula
    //
    //            +-----
    //   dx_j(u)   \      dN_k(r)
    //   -------  = +     ------- * (x_j)_k
    //     du      /       du         |
    //            +-----    |         |
    //            node k    |         |
    //                  derivative    |
    //                   of shape     |
    //                   function     |
    //                           component of
    //                          node coordinate
    //

    blitz::firstIndex  i;   // Placeholder for the first index
    blitz::secondIndex j;   // Placeholder for the second index

    blitz::Array<double,1> der_par(2);
    der_par = blitz::sum(deriv(j)*xye(i,j),j);

    // compute infinitesimal line element dr for integration along the line
    const double dr = sqrt(der_par(0)*der_par(0)+der_par(1)*der_par(1));

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)

    const double fac = intpoints.qwgt[gpid] *dr* curvefac * thsl;

    // factor given by spatial function
    double functionfac = 1.0;
    // determine coordinates of current Gauss point
    double coordgp[2];
    coordgp[0]=0.0;
    coordgp[0]=0.0;
    for (int i = 0; i< iel; i++)
    {
      coordgp[0]+=xye(0,i)*funct(i);
      coordgp[1]+=xye(1,i)*funct(i);
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;dim++)
      {
	// factor given by spatial function
	if (functions) functnum = (*functions)[dim];
	{
	  if (functnum>0)
	  {
	    // evaluate function at current gauss point
	    functionfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dim,coordgpref);
	  }
	  else
	  {
	    functionfac = 1.0;
	  }
	}
	
	elevec1[node*numdf+dim]+=
          funct(node)*(*onoff)[dim]*(*val)[dim]*fac*functionfac*invdensity;
      }
    }
  } //end of loop over integrationen points

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
GaussRule1D DRT::ELEMENTS::Fluid2Line::getOptimalGaussrule(const DiscretizationType& distype)
{
  GaussRule1D rule = intrule1D_undefined;
  switch (distype)
    {
        case line2: case nurbs2:
          rule = intrule_line_2point;
          break;
        case line3: case nurbs3:
          rule = intrule_line_3point;
          break;
    default:
    dserror("unknown number of nodes for gaussrule initialization");
    }
  return rule;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double  DRT::ELEMENTS::Fluid2Line::f2_substitution(
  const Epetra_SerialDenseMatrix  xye,
  const Epetra_SerialDenseMatrix  deriv,
  const int iel)
{
  // compute derivative of parametrization
  double dr = 0.0;
  Epetra_SerialDenseVector der_par (iel);
  der_par.Multiply('N','T',1.0,xye,deriv,0.0);
  dr=der_par.Norm2();
  return dr;
}

/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over line (public)              g.bau 07/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::IntegrateShapeFunction(ParameterList& params,
                                                       DRT::Discretization&       discretization,
                                                       vector<int>&               lm,
                                                       Epetra_SerialDenseVector&  elevec1,
                                                       const std::vector<double>& edispnp)
{
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

//  const double thsl = params.get("thsl",1.0);

/*
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;
*/

  // set number of nodes
  const int iel   = this->NumNode();

  // gaussian points
  const DiscretizationType distype = this->Shape();
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);

  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix   deriv(1,iel);

  // node coordinates
  Epetra_SerialDenseMatrix 	xye(2,iel);

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int i=0;i<iel;i++)
    {
      xye(0,i) += edispnp[3*i];
      xye(1,i) += edispnp[3*i+1];
    }
  }

  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid];
    // get shape functions and derivatives in the line
    shape_function_1D(funct,e1,distype);
    shape_function_1D_deriv1(deriv,e1,distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye,deriv,iel);

    // values are multiplied by the product from inf. area element,
    // the gauss weight, the timecurve factor and the constant
    // belonging to the time integration algorithm (theta*dt for
    // one step theta, 2/3 for bdf with dt const.)

    //double fac = intpoints.qwgt[gpid] *dr * thsl;
    const double fac = intpoints.qwgt[gpid] *dr;

    for (int node=0;node<iel;++node)
      {
       for(int dim=0;dim<numdf;dim++)
        {
          elevec1[node*numdf+dim]+=funct[node] * fac;
        }
      }
  } //end of loop over integrationen points

return;
} // DRT::ELEMENTS::Fluid2Line::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::ElementNodeNormal(ParameterList& params,
                                                  DRT::Discretization&       discretization,
                                                  vector<int>&               lm,
                                                  Epetra_SerialDenseVector&  elevec1,
                                                  const std::vector<double>& edispnp)
{
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

  // set number of nodes
  const int iel   = this->NumNode();

  // gaussian points
  const DiscretizationType distype = this->Shape();
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);

  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix   deriv(1,iel);

  // node coordinates
  Epetra_SerialDenseMatrix 	xye(2,iel);

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int i=0;i<iel;i++)
    {
      xye(0,i) += edispnp[3*i];
      xye(1,i) += edispnp[3*i+1];
    }
  }

  //this element's normal vector
  Epetra_SerialDenseVector   norm(numdf);
  double length = 0.0;
  norm[0] = xye(1,1) - xye(1,0);
  norm[1] = (-1.0)*(xye(0,1) - xye(0,0));

  length = sqrt(norm[0]*norm[0]+norm[1]*norm[1]);

  norm[0] = (1.0/length)*norm[0];
  norm[1] = (1.0/length)*norm[1];

  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid];
    // get shape functions and derivatives in the line
    shape_function_1D(funct,e1,distype);
    shape_function_1D_deriv1(deriv,e1,distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye,deriv,iel);

    // values are multiplied by the product from inf. area element and the gauss weight
    const double fac = dr * intpoints.qwgt[gpid];

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<numdf;dim++)
      {
        elevec1[node*numdf+dim]+=funct[node] * fac * norm[dim];
      }
    }
  } //end of loop over integration points
} // DRT::ELEMENTS::Fluid2Line::ElementNodeNormal

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::ElementSurfaceTension(ParameterList& params,
                                                      DRT::Discretization& discretization,
                                                      vector<int>& lm,
                                                      Epetra_SerialDenseMatrix& elemat1,
                                                      Epetra_SerialDenseVector& elevec1,
                                                      const std::vector<double>& edispnp,
                                                      std::vector<double>& enormals)
{
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

  // set number of nodes
  const int iel   = this->NumNode();

  const double thsl = params.get("thsl",0.0);
  if (thsl <= 0.0) dserror("recieved illegal time integration value");

  const double dta = params.get("dta",0.0);
  if (dta <= 0.0) dserror("recieved illegal time step value");

  // get material data
  RCP<MAT::Material> mat = parent_->Material();
  if (mat==null) dserror("no mat from parent!");
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  MATERIAL* actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

  const double density  = actmat->m.fluid->density;
  // isotropic and isothermal surface tension coefficient
  const double SFgamma = actmat->m.fluid->gamma;

  // gaussian points
  const DiscretizationType distype = this->Shape();
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);

  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector   funct(iel);
  Epetra_SerialDenseMatrix   deriv(1,iel);

  // node coordinates
  Epetra_SerialDenseMatrix 	xye(2,iel);

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }

  if (parent_->IsAle())
  {
    dsassert(edispnp.size()!=0,"paranoid");

    for (int i=0;i<iel;i++)
    {
      xye(0,i) += edispnp[3*i];
      xye(1,i) += edispnp[3*i+1];
    }
  }

  //set normal vectors to length = 1.0
  for (int node=0;node<iel;++node)
  {
    double length = 0.0;
    for (int dim=0;dim<numdf;dim++)
      length += enormals[numdf*node+dim] * enormals[numdf*node+dim];
    length = sqrt(length);
    for (int dim=0;dim<numdf;dim++)
    {
      enormals[numdf*node+dim] = (1.0/length) * enormals[numdf*node+dim];
    }
  }

  //This element's surface force vector, obtained by integrating the surface
  //stress h over the element's surface. See Wall et al. 2.3 (14)
  Epetra_SerialDenseVector   H(numdf);

  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid];
    // get shape functions and derivatives in the line
    shape_function_1D(funct,e1,distype);
    shape_function_1D_deriv1(deriv,e1,distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye,deriv,iel);

    // Values are multiplied by the product from inf. area element,
    // the gauss weight, the constant belonging to the time integration
    // algorithm (theta*dt for one step theta, 2/3 for bdf with dt const.)
    // and the inverse of the fluid's density
    const double fac = dr * intpoints.qwgt[gpid] * thsl * (1.0/density);

    //------ calculate surface tension force
    // Determinant of the metric-tensor. metric-tensor for surface element (a line) is 1x1
    double A = 0.0;
    Epetra_SerialDenseVector  a(numdf);
    for (int node=0;node<iel;++node)
    {
      for (int dim=0;dim<numdf-1;dim++)
        a[dim] += xye(dim,node)*deriv(0,node);
    }
    for (int dim=0;dim<numdf-1;dim++)
      A += a[dim]*a[dim];

    // calculate normal vector at integration point
    Epetra_SerialDenseVector  norm(numdf);
    for (int dim=0;dim<numdf;dim++)
    {
      for (int node=0;node<iel;++node)
      {
        norm[dim] += funct[node] * enormals[numdf*node+dim];
      }
    }

    // calculate surface stress at integration point
    double h = 0.0, left = 0.0, right = 0.0;
    for(int dim=0;dim<numdf;dim++)
    {
      left = 0.0;
      right = 0.0;
      for (int node=0;node<iel;++node)
      {
        left += enormals[numdf*node+dim] * deriv(0,node);
        right += xye(dim,node) * deriv(0,node);
      }
      h += left * right;
    }
    h = SFgamma * (-1.0/A) * h;

    // Gauss-Integration
    for (int dim=0;dim<numdf;dim++)
        {
          H[dim] += fac * h  * norm[dim];
        }
    //END--- calculate surface tension force

    //------ calculate element matrix contribution. Slikkerveer, Van Lohuizen, O'Brien (eq.29)
    for (int nodei = 0; nodei < iel; nodei++)
      for (int nodej = nodei; nodej < iel; nodej++)
        for (int l = 0; l < 2; l++)
          for (int p = 0; p < 2; p++)
            // it turns out that the factor 0.06 reduces calculation-time by
            // some percent compared to 0.0. much slower convergence with 1.0
            // ->don't use dta for whole timestep. dr?
            elemat1(nodei*3+l,nodej*3+p) += 0.06 * SFgamma * dta
                                            * dr * intpoints.qwgt[gpid]
                                            * deriv(0,nodei) * deriv(0,nodej)
                                            * norm[p] * norm[l];
  } //end of loop over integration points

  for (int node=0;node<iel;++node)
  {
    for(int dim=0;dim<numdf;dim++)
    {
      elevec1[numdf*node+dim] = H[dim] * (1.0/iel);
    }
  }
} // DRT::ELEMENTS::Fluid2Line::ElementSurfaceTension


#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
