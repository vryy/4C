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
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_control_point.H"
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
    else if (action == "enforce_weak_dbc")
        act = Fluid2Line::enforce_weak_dbc;
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
    case enforce_weak_dbc:
    {
      //--------------------------------------------------
      // get the condition information
      RefCountPtr<DRT::Condition> wdbc_cond 
	= 
	params.get<RefCountPtr<DRT::Condition> >("condition");

      // default is adjoint consistent
      const string* consistency  
	=
	(*wdbc_cond).Get<string>("Choice of gamma parameter");

      double wd_gamma=0.0;
      if(*consistency=="adjoint-consistent")
      {
	wd_gamma = 1.0;
      }
      else if(*consistency=="diffusive-optimal")
      {
	wd_gamma =-1.0;
      }
      else
      {
	dserror("Unknown type of definition for gamma parameter: %s",(*consistency).c_str());
      }


      // find out whether we will use a time curve
      bool usetime = true;
      const double time = params.get("total time",-1.0);
      if (time<0.0) usetime = false;

      // find out whether we will use a time curve and get the factor
      const vector<int>* curve  = (*wdbc_cond).Get<vector<int> >("curve");
      int curvenum = -1;
      if (curve) curvenum = (*curve)[0];
      double curvefac = 1.0;
      if (curvenum>=0 && usetime)
      curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      
      // get values and switches from the condition
      // (assumed to be constant on element boundary)
      const vector<int>* functions = (*wdbc_cond).Get<vector<int> >   ("funct");

      // I hope we have a linear element.
      // Ciarlet PG. The ﬁnite element method for elliptic 
      // problems. Amsterdam: North-Holland; 1978.
      if(parent_->Shape()!=Fluid2::quad4)
      {
	dserror("Cb up to now only implemented for (bi)linears");
      }
      double Cb       = 4.0;

      // get value for boundary condition
      const vector<double>* val = (*wdbc_cond).Get<vector<double> >("val");

      // get time integration parameter
      double afgdt = params.get<double>("afgdt");

      //--------------------------------------------------
      // get parent elements location vector and ownerships

      // the vectors have been allocated outside in
      // EvaluateConditionUsingParentData
      RefCountPtr<vector<int> > plm      
	=
	params.get<RefCountPtr<vector<int> > >("plm");
      RefCountPtr<vector<int> > plmowner 
	=
	params.get<RefCountPtr<vector<int> > >("plmowner");
      
      parent_->LocationVector(discretization,*plm,*plmowner);

      // Reshape element matrices and vectors and init to zero
      const int eledim = (int)(*plm).size();
      elemat1.Shape(eledim,eledim);
      elemat2.Shape(eledim,eledim);
      elevec1.Size(eledim);
      elevec2.Size(eledim);
      elevec3.Size(eledim);

      // --------------------------------------------------
      // extract velocities from global distributed vectors

      // velocities (intermediate time step, n+alpha_F)
      RefCountPtr<const Epetra_Vector> velaf 
	=
	discretization.GetState("u and p (n+alpha_F,trial)");
      if (velaf==null)
      {
        dserror("Cannot get state vector 'velaf'");
      }
      vector<double> mypvelaf((*plm).size());
      DRT::UTILS::ExtractMyValues(*velaf,mypvelaf,*plm);

      // create blitz matrix object for convenience
      const int numnode = parent_->NumNode();
      blitz::Array<double, 2> pevelaf(2,
				      numnode,
				      blitz::ColumnMajorArray<2>());

      // extract velocities
      for (int i=0;i<numnode;++i)
      {
	pevelaf(0,i) = mypvelaf[0+(i*3)];
	pevelaf(1,i) = mypvelaf[1+(i*3)];
      }

      // call the special evaluation method
      EvaluateWeakDirichlet(lm            ,
			    *plm          ,
			    elemat1       ,
			    elevec1       ,
			    pevelaf       ,
			    *val          ,
			    functions     ,
			    curvefac      ,
			    Cb            ,
			    wd_gamma      ,
			    afgdt         );

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

  // get time parameter
  const double thsl = params.get("thsl",0.0);

  // get constant density (only relevant for incompressible flow)
  //const double inc_dens = params.get("inc_density",0.0);

  // get flag whether outflow stabilization or not
  string outflow_stab = params.get("outflow stabilization","no_outstab");

  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel   = this->NumNode();

  // Gaussian points
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
      DRT::NURBS::ControlPoint* cp =
         dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights(inode) = cp->W();
    }
  }

  // allocate vector for shape functions and matrix for derivatives
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
    edensnp(i) = myvedenp[2+(i*3)];
  }

  // this part will be run when outflow stabilization is required
  if(outflow_stab == "yes_outstab")
  {
    // Determine normal to this element
    std::vector<double> normal(2);
    double length;

    normal[0] = xye(1,1) - xye(1,0);
    normal[1] = (-1.0)*(xye(0,1) - xye(0,0));

    length = sqrt(normal[0]*normal[0]+normal[1]*normal[1]);

    // outward pointing normal of length 1
    for (int i=0; i<2; i++) normal[i] = normal[i] / length;

    RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
    if (velnp==null) dserror("Cannot get state vector 'velnp'");

    // extract local values from global vector
    vector<double> myvelnp(lm.size());
    DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

    // create blitz object for element array
    blitz::Array<double, 2> evelnp(2,iel,blitz::ColumnMajorArray<2>());

    // insert velocity into element array
    for (int i=0;i<iel;++i)
    {
      evelnp(0,i) = myvelnp[0+(i*3)];
      evelnp(1,i) = myvelnp[1+(i*3)];
    }

    /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
    *----------------------------------------------------------------------*/
    for (int gpid=0; gpid<intpoints.nquad; gpid++)
    {
      const double e1 = intpoints.qxg[gpid];

      // get shape functions and derivatives for linear element
      DRT::UTILS::shape_function_1D(funct,e1,distype);
      // explicit determination of first derivatives, since respective function
      //DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);
      // provides mismatch due to matrix/vector inconsistency
      switch (distype)
      {
        case DRT::Element::point1:
        {
          deriv(0) = 0.0;
          break;
        }
        case DRT::Element::line2:
        {
          deriv(0)= -0.5;
          deriv(1)= 0.5;
          break;
        }
        case DRT::Element::line3:
        {
          deriv(0)= e1 - 0.5;
          deriv(1)= e1 + 0.5;
          deriv(2)= -2.0*e1;
          break;
        }
        default:
           dserror("distype unknown\n");
      }

      std::vector<double> vel(2);
      double normvel=0;
      for(int dim=0;dim<2;dim++)
      {
        vel[dim] = 0;
        for (int node=0;node<iel;++node)
        {
          vel[dim]+= funct(node) * evelnp(dim,node);
        }
        normvel += vel[dim]*normal[dim];
      }

      if(normvel<-0.0001)
      {
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

        const double fac = intpoints.qwgt[gpid] * dr * thsl * normvel;

        for (int node=0;node<iel;++node)
        {
          for(int dim=0;dim<2;dim++)
          {
            elevec1[node*numdf+dim]+= funct(node) * edensnp(node) * evelnp(dim,node) * fac;
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

    // loop over integration points
    for (int gpid=0;gpid<intpoints.nquad;gpid++)
    {
      const double e1 = intpoints.qxg[gpid];

      // get shape functions and derivatives for line element
      if(!(distype == DRT::Element::nurbs2 || distype == DRT::Element::nurbs3))
      {
        DRT::UTILS::shape_function_1D(funct,e1,distype);
        // explicit determination of first derivatives, since respective function
        //DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);
        // provides mismatch due to matrix/vector inconsistency
        switch (distype)
        {
        case DRT::Element::point1:
        {
          deriv(0) = 0.0;
          break;
        }
        case DRT::Element::line2:
        {
          deriv(0)= -0.5;
          deriv(1)= 0.5;
          break;
        }
        case DRT::Element::line3:
        {
          deriv(0)= e1 - 0.5;
          deriv(1)= e1 + 0.5;
          deriv(2)= -2.0*e1;
          break;
        }
        default:
           dserror("distype unknown\n");
        }
      }
      else
        DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv(funct,deriv,e1,myknots,weights,distype);

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
      // Furthermore, there may be a divison by the constant scalar density,
      // only relevant (i.e., it may be unequal 1.0) in incompressible flow case
      const double fac = intpoints.qwgt[gpid] *dr* curvefac * thsl;

      // factor given by spatial function
      double functfac = 1.0;
      // determine coordinates of current Gauss point
      double coordgp[2];
      coordgp[0]=0.0;
      coordgp[1]=0.0;
      for (int i = 0; i< iel; i++)
      {
        coordgp[0]+=xye(0,i)*funct(i);
        coordgp[1]+=xye(1,i)*funct(i);
      }

      int functnum = -1;
      const double* coordgpref = &coordgp[0]; // needed for function evaluation

      for (int node=0;node<iel;++node)
      {
        for(int dim=0;dim<2;dim++)
        {
          // factor given by spatial function
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

          elevec1[node*numdf+dim]+= edensnp(node)*funct(node)*(*onoff)[dim]*(*val)[dim]*fac*functfac;
        }
      }
    } //end of loop over integrationen points
  }

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
    // the gauss weight and the constant belonging to the time integration
    // algorithm (theta*dt for one step theta, 2/3 for bdf with dt const.)
    const double fac = dr * intpoints.qwgt[gpid] * thsl;

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


/*--------------------------------------------------------------------
  Evaluate weakly imposed Dirichlet conditions
  --------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::EvaluateWeakDirichlet(
  vector<int>&               lm            ,
  vector<int>&               plm           ,
  Epetra_SerialDenseMatrix&  elemat        ,
  Epetra_SerialDenseVector&  elevec        ,
  blitz::Array<double,2>&    pevelaf       ,
  vector<double>             val           ,
  const vector<int>*         functions     ,
  double                     curvefac      ,
  double                     Cb            ,
  double                     wd_gamma      ,
  double                     afgdt
  )
{
  //--------------------------------------------------
  // gausspoint quantites
  blitz::Array<double,1> velintaf(2);
  blitz::Array<double,2> vderxyaf(2,2,blitz::ColumnMajorArray<2>());

  blitz::Array<double,2> pxji(2,2,blitz::ColumnMajorArray<2>());
  blitz::Array<double,2> pxjm(2,2,blitz::ColumnMajorArray<2>());


  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------
  const int piel = parent_->NumNode();

  //--------------------------------------------------
  // get node coordinates
  blitz::Array<double,2> pxye(2,piel);

  for(int i=0;i<piel;++i)
  {
    pxye(0,i)=parent_->Nodes()[i]->X()[0];
    pxye(1,i)=parent_->Nodes()[i]->X()[1];
  }

  //--------------------------------------------------
  // allocate arrays for shapefunctions and derivatives
  blitz::Array<double,1> pfunct(  parent_->NumNode());
  blitz::Array<double,2> pderiv(2,parent_->NumNode());

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = parent_->Material();

  if( mat->MaterialType()    != m_carreauyasuda
      && mat->MaterialType() != m_modpowerlaw
      && mat->MaterialType() != m_fluid)
          dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(mat->MaterialType() == m_fluid)
  {
    MATERIAL* actmat 
      =  
      static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

    visc = actmat->m.fluid->viscosity;
  }
  else
  {
    dserror("up to now I expect a constant viscosity to inforce weak DBCs\n");
  }
  
  //--------------------------------------------------
  //          GET BOUNDARY ELEMENT DATA
  //--------------------------------------------------

  //--------------------------------------------------
  // get node coordinates
  blitz::Array<double,2> xye(2,NumNode());

  for(int i=0;i<NumNode();++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }

  //--------------------------------------------------
  // allocate arrays for shapefunctions and derivatives
  blitz::Array<double,1> funct(  NumNode());
  blitz::Array<double,2> deriv(1,NumNode());

  //--------------------------------------------------
  // array for an outward normal
  blitz::Array<double,1> n(2);

  //--------------------------------------------------
  // get gausspoints to integrate over boundary element

  // get gauss rule
  GaussRule1D intrule=DRT::UTILS::intrule1D_undefined;
  switch (Shape())
  {
  case line2:
  {
    intrule=DRT::UTILS::intrule_line_2point;
    break;
  }
  default:
    dserror("invalid discretization type for fluid2line weak DBC evaluation");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(intrule);

  //--------------------------------------------------
  // the gausspoints above have to be mapped to the 
  // parent element to be able to evaluate one sided
  // derivatives on the boundary
  
  // check if parent is a quad4
  if(parent_->Shape()!=Fluid2::quad4)
  {
    dserror("Dont know how to determine gausspoint in parent element for other elements than quad4\n");
  }
  
  // startnode is the first node of the edge --- we use it to determine the gausspoints
  // coordinates in the global element
  int startnode=0;
  for(unsigned rr=0;rr<plm.size();++rr)
  {
    if(plm[rr]==lm[0])
    {
      break;
    }
    ++startnode;
  }

  if(startnode%3!=0)
  {
    dserror("problem with dof numbering when trying to determine my start node");
  }
  else
  {
    startnode=startnode/3;
  }
  // or maybe use a sweet blitz array/Epetra vector here? ;-) greetings Axel
  vector<vector<double> > intpointinparent(intpoints.nquad, std::vector<double>(2,0.0));
  //double intpointinparent[intpoints.nquad][2];
  switch(startnode)
  {
  case 0:
  {
    /*                s|
                       | 

             3                   2
              +-----------------+
              |                 |
              |                 |
	      |                 |
	      |        |        |             r
              |        +--      |         -----
              |                 |
              |                 |
              |                 |
              |                 |
              +-----------*-----+
             0                   1
                    -->|gp|<--               */

    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      intpointinparent[iquad][0]=intpoints.qxg [iquad];
      intpointinparent[iquad][1]=-1.0;
    }
    break;
  }
  case 1:
  {
    /*                s|
                       | 

             3                   2
              +-----------------+
              |                 | |
              |                 | v
	      |                 *---
	      |        |        | gp          r
              |        +--      |---      -----
              |                 | ^
              |                 | |
              |                 |
              |                 |
              +-----------------+
             0                   1
                                             */

    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      intpointinparent[iquad][0]= 1.0;
      intpointinparent[iquad][1]=intpoints.qxg [iquad];
    }
    break;
  }  
  case 2:
  {
    /*                s|
                       | 

             3   -->|gp|<--       
              +-----*-----------+
              |                 |  
              |                 |  
	      |                 |   
	      |        |        |             r
              |        +--      |         -----
              |                 |  
              |                 |  
              |                 |
              |                 |
              +-----------------+
             0                   1
                                             */

    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      intpointinparent[iquad][0]=-intpoints.qxg [iquad];
      intpointinparent[iquad][1]= 1.0;
    }
    break;
  }  
  case 3:
  {
    /*                s|
                       | 

             3                    
              +-----*-----------+
              |                 |  
              |                 |  
	    | |                 |   
	    v |        |        |             r
           ---|        +--      |         -----
            gp|                 |  
           ---*                 |  
            ^ |                 |
            | |                 |
              +-----------------+
             0                   1
                                             */

    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      intpointinparent[iquad][0]=-1.0;
      intpointinparent[iquad][1]=-intpoints.qxg [iquad];
    }
    break;
  }
  default:
  {
    dserror("unable to determine intpoint in parent");
  }
  }

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (int iquad=0;iquad<intpoints.nquad;++iquad)
  {
    // gaussian weight
    const double wquad = intpoints.qwgt[iquad];

    // gaussian point in boundary elements local coordinates
    const double gp    = intpoints.qxg [iquad];

    // gaussian point in parent elements local coordinates
    const double r     = intpointinparent[iquad][0];
    const double s     = intpointinparent[iquad][1];

    // shape function derivs of boundary element at gausspoint
    DRT::UTILS::shape_function_1D       (funct,gp,Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gp,Shape());

    // compute infinitesimal line element dr
    /*
    //      /     \    +---+          /     \
    //  d  |  x(r) |    \     dN_i   |  x_i  |
    // --- |       | =   +    ---- * |       |
    // dr  |  y(r) |    /      dr    |  y_i  |
    //      \     /    +---+          \     / 
    //                   i
    */
    blitz::Array<double,1> dxydr(2);
    dxydr=0;
    
    for(int i=0;i<NumNode();++i)
    {
      dxydr(0)+=xye(0,i)*deriv(0,i);
      dxydr(1)+=xye(1,i)*deriv(0,i);
    }

    const double dr = sqrt(dxydr(0)*dxydr(0)+dxydr(1)*dxydr(1));

    if (dr<0.0)
    {
      dserror("negative line element\n");
    }

    // compute normal
    /*
      //      /    \           /   dy  \
      //     |  n   |         |    --   |
      //     |   1  |   1.0   |    dr   |
      //     |      | = --- * |         |
      //     |  n   |    dr   |    dx   |
      //     |   2  |         |  - --   |
      //      \    /           \   dr  /
      // 
    */

    // compute normal
    n(0)= dxydr(1)/dr;
    n(1)=-dxydr(0)/dr;

    // factor given by spatial function
    double functionfac[2];
    functionfac[0]= 1.0;
    functionfac[1]= 1.0;

    // determine coordinates of current Gauss point
    double coordgp[2];
    coordgp[0]=0.0;
    coordgp[1]=0.0;
    for (int i = 0; i< NumNode(); i++)
    {
      coordgp[0]+=xye(0,i)*funct(i);
      coordgp[1]+=xye(1,i)*funct(i);
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation

    for (int node=0;node<NumNode();++node)
    {
      for(int dim=0;dim<2;dim++)
      {
	// factor given by spatial function
	if (functions) 
	{
	  functnum = (*functions)[dim];
	  if (functnum>0)
	  {
	    // evaluate function at current gauss point
	    functionfac[dim] = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dim,coordgpref);
	  }
	  else
	  {
	    functionfac[dim] = 1.0;
	  }
	}
      }
    }

    // shape functions and derivs of corresponding parent at gausspoint
    DRT::UTILS::shape_function_2D       (pfunct,r,s,parent_->Shape());
    DRT::UTILS::shape_function_2D_deriv1(pderiv,r,s,parent_->Shape());

    //-----------------------------------------------------
    //
    //                       +-------------+       
    //                      / /  T       \ |
    //           h = 2 * \ / |  n * G * n |
    //            b       +   \          /
    //              

    // get Jacobian matrix and determinant
    pxjm=0;

    for(int i=0;i<piel;++i)
    {
      for(int rr=0;rr<2;++rr)
      {
	for(int mm=0;mm<2;++mm)
	{
	  pxjm(rr,mm)+=pderiv(rr,i)*pxye(mm,i);
	}
      }
    }
    
    const double pdet = pxjm(0,0)*pxjm(1,1) - pxjm(0,1)*pxjm(1,0);

    // check for degenerated elements
    if (pdet < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", parent_->Id(), pdet);
    }

    /*
      Use the Jacobian and the known derivatives in element coordinate
      directions on the right hand side to compute the derivatives in
      global coordinate directions

          +-          -+     +-    -+      +-    -+
          |  dx    dy  |     | dN_k |      | dN_k |
          |  --    --  |     | ---- |      | ---- |
          |  dr    dr  |     |  dx  |      |  dr  |
          |            |  *  |      |   =  |      | for all k
          |  dx    dy  |     | dN_k |      | dN_k |
          |  --    --  |     | ---- |      | ---- |
          |  ds    ds  |     |  dy  |      |  ds  |
          +-          -+     +-    -+      +-    -+

          Matrix is inverted analytically
    */
    // inverse of jacobian
    pxji(0,0) = ( pxjm(1,1))/pdet;
    pxji(0,1) = (-pxjm(0,1))/pdet;
    pxji(1,0) = (-pxjm(1,0))/pdet;
    pxji(1,1) = ( pxjm(0,0))/pdet;

    
    /*          +-           -+   +-           -+   +-           -+
		|             |   |             |   |             |
		|  dr    dr   |   |  ds    ds   |   |  dt    dt   |
          G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
           ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
	        |    i     j  |   |    i     j  |   |    i     j  |
		+-           -+   +-           -+   +-           -+
    */
    blitz::Array<double,2> G(2,2,blitz::ColumnMajorArray<2>());

    for (int nn=0;nn<2;++nn)
    {
      for (int rr=0;rr<2;++rr)
      {
	G(nn,rr) = pxji(nn,0)*pxji(rr,0)+pxji(nn,1)*pxji(rr,1);
      }
    }

    //
    //                           2.0
    //             h  = ---------------------
    //              b        +-------------+       
    //                      / /  T       \ |
    //                   \ / |  n * G * n |
    //                    +   \          /
    //              

    const double h =2.0/sqrt(n(0)*G(0,0)*n(0)
			     +
			     n(0)*G(0,1)*n(1)
			     +
			     n(1)*G(1,0)*n(0)
			     +
			     n(1)*G(1,1)*n(1));


    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----
    //       n+af       \                  n+af
    //    vel    (x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    velintaf=0;
    for(int j=0;j<piel;++j)
    {
      velintaf(0)+=pfunct(j)*pevelaf(0,j);
      velintaf(1)+=pfunct(j)*pevelaf(1,j);
    }

    // get velocity (n+alpha_F,i) derivatives at integration point
    //
    //       n+af      +-----  dN (x)
    //   dvel    (x)    \        k         n+af
    //   ----------- =   +     ------ * vel
    //       dx         /        dx        k
    //         j       +-----      j
    //                 node k
    //
    // j : direction of derivative x/y/z
    //
    vderxyaf=0;
    for(int k=0;k<piel;++k)
    {
      for (int j=0;j<2;++j)
      {
	vderxyaf(0,j)+=pderiv(j,k)*pevelaf(0,k);
	vderxyaf(1,j)+=pderiv(j,k)*pevelaf(1,k);
      }
    }

    //--------------------------------------------------
    // partially integrated viscous term
    /*
    // factor: 2*nu*afgdt
    //
    //    /                     \
    //   |           s           |
    // - |  v , nabla  Dacc * n  |
    //   |                       |
    //    \                     / boundaryele
    //
    */
    for (int ui=0; ui<piel; ++ui) 
    {
      for (int vi=0; vi<piel; ++vi) 
      {
	elemat[vi*3    ][ui*3    ] -= dr*wquad*2.0*visc*afgdt*pfunct(vi)*
	  (pderiv(0,ui)*n(0)+0.5*pderiv(1,ui)*n(1));
        elemat[vi*3    ][ui*3 + 1] -= dr*wquad*2.0*visc*afgdt*pfunct(vi)*
	  (0.5*pderiv(0,ui)*n(1));
	elemat[vi*3 + 1][ui*3    ] -= dr*wquad*2.0*visc*afgdt*pfunct(vi)*
	  (0.5*pderiv(1,ui)*n(0));
	elemat[vi*3 + 1][ui*3 + 1] -= dr*wquad*2.0*visc*afgdt*pfunct(vi)*
	  (0.5*pderiv(0,ui)*n(0)+pderiv(1,ui)*n(1));
      }
    }
    /*
    // factor: 2*nu
    //
    //    /                     \
    //   |           s  n+af     |
    // + |  v , nabla  u    * n  |
    //   |                       |
    //    \                     / boundaryele
    //
    */
    for (int vi=0; vi<piel; ++vi) 
    {
      elevec[vi*3    ] += dr*wquad*2.0*visc*pfunct(vi)*
	(vderxyaf(0,0)*n(0)
	 +
	 0.5*(vderxyaf(0,1)+vderxyaf(1,0))*n(1));
      elevec[vi*3 + 1] += dr*wquad*2.0*visc*pfunct(vi)*
	(0.5*(vderxyaf(0,1)+vderxyaf(1,0))*n(0)
	 +
	 vderxyaf(1,1)*n(1));
    }

    //--------------------------------------------------
    // (adjoint) consistency term, viscous part

    /*
    // factor: 2*nu*gamma_wd*afgdt
    //
    //    /                     \
    //   |        s              |
    // - |   nabla  w * n , Dacc |
    //   |                       |
    //    \                     / boundaryele
    //
    */
    for (int ui=0; ui<piel; ++ui) 
    {
      for (int vi=0; vi<piel; ++vi) 
      {
	elemat[vi*3    ][ui*3    ] -= dr*wquad*2.0*visc*wquad*wd_gamma*afgdt*
	  (pderiv(0,vi)*n(0)+0.5*pderiv(1,vi)*n(1))*pfunct(ui);
        elemat[vi*3    ][ui*3 + 1] -= dr*wquad*2.0*visc*wquad*wd_gamma*afgdt*
	  (0.5*pderiv(1,vi)*n(0)                  )*pfunct(ui);
	elemat[vi*3 + 1][ui*3    ] -= dr*wquad*2.0*visc*wquad*wd_gamma*afgdt*
	  (0.5*pderiv(0,vi)*n(1)                  )*pfunct(ui);
	elemat[vi*3 + 1][ui*3 + 1] -= dr*wquad*2.0*visc*wquad*wd_gamma*afgdt*
	  (0.5*pderiv(0,vi)*n(0)+pderiv(1,vi)*n(1))*pfunct(ui);
      }
    }
    /*
    // factor: 2*nu*gamma_wd
    //
    //    /                           \
    //   |        s          n+af      |
    // + |   nabla  w * n , u    - u   |
    //   |                          b  |
    //    \                           / boundaryele
    //
    */
    for (int vi=0; vi<piel; ++vi) 
    {
      elevec[vi*3    ] += dr*wquad*2.0*visc*wquad*wd_gamma*
	((pderiv(0,vi)*n(0)+0.5*pderiv(1,vi)*n(1))*(velintaf(0)-val[0]*functionfac[0]*curvefac)
	 +
	 (0.5*pderiv(1,vi)*n(0)                  )*(velintaf(1)-val[1]*functionfac[1]*curvefac));

      elevec[vi*3 + 1] += dr*wquad*2.0*visc*wquad*wd_gamma*
	((0.5*pderiv(0,vi)*n(1)                  )*(velintaf(0)-val[0]*functionfac[0]*curvefac)
	 +
	 (0.5*pderiv(0,vi)*n(0)+pderiv(1,vi)*n(1))*(velintaf(1)-val[1]*functionfac[1]*curvefac));
    }

    //--------------------------------------------------
    // penalty term

    /*
    // factor: nu*Cb/h*afgdt
    //
    //    /          \
    //   |            |
    // + |  w , Dacc  |
    //   |            |
    //    \          / boundaryele
    //
    */

    for (int ui=0; ui<piel; ++ui) 
    {
      for (int vi=0; vi<piel; ++vi) 
      {
	elemat[vi*3    ][ui*3    ] += afgdt*dr*wquad*Cb*visc/h*pfunct(ui)*pfunct(vi);
	elemat[vi*3 + 1][ui*3 + 1] += afgdt*dr*wquad*Cb*visc/h*pfunct(ui)*pfunct(vi);
      }
    }

    /*
    // factor: nu*Cb/h
    //
    //    /                \
    //   |        n+af      |
    // + |   w , u    - u   |
    //   |               b  |
    //    \                / boundaryele
    //
    */

    for (int vi=0; vi<piel; ++vi) 
    {
      elevec[vi*3    ] -= dr*wquad*Cb*visc/h*pfunct(vi)*(velintaf(0)-val[0]*functionfac[0]*curvefac);
      elevec[vi*3 + 1] -= dr*wquad*Cb*visc/h*pfunct(vi)*(velintaf(1)-val[1]*functionfac[1]*curvefac);
    }

    //--------------------------------------------------
    // adjoint consistency term, convective part 
    
    double flux=velintaf(0)*n(0)+velintaf(1)*n(1);

    if(flux>0)
    {
      flux=0;
    }
    else
    {
      /*
	// This linearisation has only to be included if
	// u*n is negative --- otherwise it's nonesense
	//
	// factor: afgdt
	//
	//    /                             \
	//   |    /        \       n+af      |
	// - |   | Dacc * n | w , u    - u   |
	//   |    \        /              b  |
	//    \                             / boundaryele, inflow
	//               
      */

      for (int ui=0; ui<piel; ++ui) 
      {
	for (int vi=0; vi<piel; ++vi) 
	{
	  elemat[vi*3    ][ui*3    ] -= afgdt*dr*wquad*(velintaf(0)-val[0]*functionfac[0]*curvefac)*pfunct(ui)*n(0)*pfunct(vi);
	  elemat[vi*3    ][ui*3 + 1] -= afgdt*dr*wquad*(velintaf(0)-val[0]*functionfac[0]*curvefac)*pfunct(ui)*n(1)*pfunct(vi);
	  elemat[vi*3 + 1][ui*3    ] -= afgdt*dr*wquad*(velintaf(1)-val[1]*functionfac[1]*curvefac)*pfunct(ui)*n(0)*pfunct(vi);
	  elemat[vi*3 + 1][ui*3 + 1] -= afgdt*dr*wquad*(velintaf(1)-val[1]*functionfac[1]*curvefac)*pfunct(ui)*n(1)*pfunct(vi);
	}
      }
    }

    /*
    // factor: afgdt
    //
    //    /                       \
    //   |    / n+af   \           |
    // - |   | u    * n | w , Dacc |
    //   |    \        /           |
    //    \  |          |         / boundaryele, inflow
    //       +----------+
    //           <0
    */


    for (int ui=0; ui<piel; ++ui) 
    {
      for (int vi=0; vi<piel; ++vi) 
      {
	elemat[vi*3    ][ui*3    ] -= afgdt*dr*wquad*flux*pfunct(ui)*pfunct(vi);
	elemat[vi*3 + 1][ui*3 + 1] -= afgdt*dr*wquad*flux*pfunct(ui)*pfunct(vi);
      }
    }

    /*
    // factor: 1
    //
    //    /                             \
    //   |    / n+af   \       n+af      |
    // - |   | u    * n | w , u    - u   |
    //   |    \        /              b  |
    //    \  |          |               / boundaryele, inflow
    //       +----------+
    //           <0
    */

    for (int vi=0; vi<piel; ++vi) 
    {
      elevec[vi*3    ] += dr*wquad*pfunct(vi)*flux*(velintaf(0)-val[0]*functionfac[0]*curvefac);
      elevec[vi*3 + 1] += dr*wquad*pfunct(vi)*flux*(velintaf(1)-val[1]*functionfac[1]*curvefac);
    }

  } // end the gausspointloop

  return;
} // end EvaluateWeakDirichlet

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
