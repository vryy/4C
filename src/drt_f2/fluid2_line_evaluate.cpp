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
#include "fluid2_weak_dbc.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_control_point.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/sutherland_fluid.H"
#include "../drt_mat/carreauyasuda.H"
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
    else if (action == "calc_Neumann_inflow")
        act = Fluid2Line::calc_Neumann_inflow;
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
    case enforce_weak_dbc:
    {
      return DRT::ELEMENTS::Fluid2LineWeakDBCInterface::Impl(this)->EvaluateWeakDBC(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
      break;

#if 0
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
      
      // there is only one definition of tauB available in 2d
      const string* deftauB
        =
        (*wdbc_cond).Get<string>("Definition of penalty parameter");
      
      if(*deftauB!="constant")
      {
        dserror("Definition of penalty parameter tauB not supported in 2d: %s",(*deftauB).c_str());
      }

      // linearisation of adjoint convective flux
      const string* linearisation_approach
        =
        (*wdbc_cond).Get<string>("Linearisation");
      
      bool complete_linearisation=false;
      
      if(*linearisation_approach=="lin_all")
      {
        complete_linearisation=true;
      }
      else if(*linearisation_approach=="no_lin_conv_inflow")
      {
        complete_linearisation=false;
      }
      else
      {
        dserror("Unknown definition of linearisation approach: %s",(*linearisation_approach).c_str());
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
      // Ciarlet PG. The ï¬nite element method for elliptic
      // problems. Amsterdam: North-Holland; 1978.
      if(parent_->Shape()!=Fluid2::quad4)
      {
	dserror("Cb up to now only implemented for (bi)linears");
      }
      double Cb       = 4.0;

      // get value for boundary condition
      const vector<double>* val = (*wdbc_cond).Get<vector<double> >("val");

      // get time integration parameters
      const double afgdt = params.get<double>("afgdt");

      const double gdt   = params.get<double>("gdt");

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
      RefCountPtr<const Epetra_Vector> velnp 
	=
	discretization.GetState("u and p (n+1      ,trial)");

      if (velaf==null || velnp==null)
      {
        dserror("Cannot get state vector 'velaf', 'velnp'");
      }
      vector<double> mypvelaf((*plm).size());
      DRT::UTILS::ExtractMyValues(*velaf,mypvelaf,*plm);

      vector<double> mypvelnp((*plm).size());
      DRT::UTILS::ExtractMyValues(*velnp,mypvelnp,*plm);

      // create  matrix object for convenience
      const int numnode = parent_->NumNode();
      Epetra_SerialDenseMatrix pevelaf(2,numnode);
      Epetra_SerialDenseMatrix pevelnp(2,numnode);
      Epetra_SerialDenseVector peprenp(  numnode);

      // extract velocities
      for (int i=0;i<numnode;++i)
      {
	pevelaf(0,i) = mypvelaf[0+(i*3)];
	pevelaf(1,i) = mypvelaf[1+(i*3)];

	pevelnp(0,i) = mypvelnp[0+(i*3)];
	pevelnp(1,i) = mypvelnp[1+(i*3)];
        peprenp(i)   = mypvelnp[2+(i*3)];
      }

      // call the special evaluation method
      EvaluateWeakDirichlet(lm                    ,
			    *plm                  ,
			    elemat1               ,
			    elevec1               ,
			    pevelaf               ,
			    pevelnp               ,
                            peprenp               ,
			    *val                  ,
			    functions             ,
			    curvefac              ,
			    Cb                    ,
			    wd_gamma              ,
			    afgdt                 ,
                            gdt                   ,
                            complete_linearisation
        );

#endif
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
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

  // get time parameter
  const double thsl = params.get("thsl",0.0);

  // get constant density (only relevant for incompressible flow)
  //const double inc_dens = params.get("inc_density",0.0);

  // get discretization type
  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel = this->NumNode();

  // Gaussian points
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints(gaussrule);

  // --------------------------------------------------
  // Now initialise the nurbs specific stuff
  Epetra_SerialDenseVector myknots;

  Epetra_SerialDenseVector weights(iel);

  // for isogeometric elements
  if(this->Shape()==nurbs2 || this->Shape()==nurbs3)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    std::vector<Epetra_SerialDenseVector> surfaceknots(2);
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
    myknots.Size((surfaceknots[(Id())%2]).Length());
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
  Epetra_SerialDenseVector funct(iel);
  Epetra_SerialDenseVector deriv(iel);

  // node coordinates
  Epetra_SerialDenseMatrix xye(2,iel);

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

  // create object for density array
  Epetra_SerialDenseVector edensnp(iel);

  // insert density into element array
  for (int i=0;i<iel;++i)
  {
    edensnp(i) = myvedenp[2+(i*3)];
  }

  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid][0];

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
    Epetra_SerialDenseVector der_par(2);
    for(int rr=0;rr<2;++rr)
    {
      der_par(rr)=0;
      for(int mm=0;mm<iel;++mm)
      {
        der_par(rr) += deriv(mm)*xye(rr,mm);
      }
    }

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
            functfac = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dim,coordgpref,0.0,NULL);
          }
          else
            functfac = 1.0;
        }

        elevec1[node*numdf+dim]+= edensnp(node)*funct(node)*(*onoff)[dim]*(*val)[dim]*fac*functfac;
      }
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | compute potential Neumann inflow                            vg 03/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::NeumannInflow(
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

  // get discretization type
  const DiscretizationType distype = this->Shape();

  // set number of nodes
  const int iel = this->NumNode();

  // Gaussian points
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints(gaussrule);

  // (density-weighted) shape functions and first derivatives
  Epetra_SerialDenseVector funct(iel);
  Epetra_SerialDenseVector densfunct(iel);
  Epetra_SerialDenseVector deriv(iel);

  // node coordinates
  Epetra_SerialDenseMatrix xye(2,iel);

  // the element's normal vector
  Epetra_SerialDenseVector normal(2);

  // velocity and momentum at gausspoint
  Epetra_SerialDenseVector velint(2);
  Epetra_SerialDenseVector momint(2);

  // Jacobian
  Epetra_SerialDenseVector der_par(2);

  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=this->Nodes()[i]->X()[0];
    xye(1,i)=this->Nodes()[i]->X()[1];
  }

  // determine outward-pointing normal to this element
  normal(0) = xye(1,1)-xye(1,0);
  normal(1) = (-1.0)*(xye(0,1)-xye(0,0));

  // length of normal
  double length = 0.0;
  length = sqrt(normal(0)*normal(0)+normal(1)*normal(1));

  // outward-pointing normal of unit length
  for(int inode=0;inode<2;inode++)
  {
    normal(inode) = normal(inode)/length;
  }

  // get velocity and density vector
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  RCP<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  if (velnp==null or vedenp==null)
    dserror("Cannot get state vector 'velnp' and/or 'vedenp'");

  // extract local values from global vector
  vector<double> myvelnp(lm.size());
  vector<double> myvedenp(lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);

  // create blitz object for density array
  Epetra_SerialDenseMatrix evelnp(2,iel);
  Epetra_SerialDenseVector edensnp(iel);

  // insert velocity and density into element array
  for (int i=0;i<iel;++i)
  {
    evelnp(0,i) = myvelnp[0+(i*3)];
    evelnp(1,i) = myvelnp[1+(i*3)];

    edensnp(i) = myvedenp[2+(i*3)];
  }

  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.nquad; gpid++)
  {
    const double e = intpoints.qxg[gpid][0];

    // get shape functions
    DRT::UTILS::shape_function_1D(funct,e,distype);

    // compute momentum (i.e., density times velocity) and normal momentum
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double normmom = 0.0;
    for (int j=0;j<2;++j)
    {
      velint(j) = 0.0;
      momint(j) = 0.0;
      for (int i=0;i<iel;++i)
      {
        velint(j) += funct(i)*evelnp(j,i);
        momint(j) += edensnp(i)*velint(j);
      }
      normmom += momint(j)*normal(j);
    }

    // computation only required for negative normal momentum
    if (normmom<-0.0001)
    {
      // first derivatives
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
          deriv(0)= e - 0.5;
          deriv(1)= e + 0.5;
          deriv(2)= -2.0*e;
          break;
        }
        default:
           dserror("distype unknown\n");
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
      for(int rr=0;rr<2;++rr)
      {
        der_par(rr)=0;
        for(int mm=0;mm<iel;++mm)
        {
          der_par(rr) += deriv(mm)*xye(rr,mm);
        }
      }

      // compute infinitesimal line element dr for integration along the line
      const double dr = sqrt(der_par(0)*der_par(0)+der_par(1)*der_par(1));

      // integration factor
      const double fac = intpoints.qwgt[gpid] * dr;

      // integration factor for left- and right-hand side
      const double lhsfac = normmom*timefac*fac;
      double rhsfac = normmom*fac;
      if (not is_genalpha) rhsfac *= timefac;

      // matrix
      for (int vi=0; vi<iel; ++vi)
      {
        const double vlhs = lhsfac*funct(vi);

        const int fvi   = 3*vi;
        const int fvip  = fvi+1;

        for (int ui=0; ui<iel; ++ui)
        {
          const int fui   = 3*ui;
          const int fuip  = fui+1;

          elemat1(fvi  ,fui  ) -= vlhs*funct(ui);
          elemat1(fvip ,fuip ) -= vlhs*funct(ui);
        }
      }

      // rhs
      const double vrhs0 = rhsfac*velint(0);
      const double vrhs1 = rhsfac*velint(1);
      for (int vi=0; vi<iel; ++vi)
      {
        const int fvi   =3*vi;
        const int fvip  =fvi+1;

        elevec1(fvi  ) += funct(vi)*vrhs0;
        elevec1(fvip ) += funct(vi)*vrhs1;
      }
    }
  }

  return;
}// DRT::ELEMENTS::Fluid2Line::NeumannInflow


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
  const IntegrationPoints1D  intpoints(gaussrule);

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
    const double e1 = intpoints.qxg[gpid][0];
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
  const IntegrationPoints1D  intpoints(gaussrule);

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
      xye(0,i) += edispnp[numdf*i];
      xye(1,i) += edispnp[numdf*i+1];
    }
  }

  //this element's normal vector
  Epetra_SerialDenseVector   norm(2);
  double length = 0.0;
  norm[0] = xye(1,1) - xye(1,0);
  norm[1] = (-1.0)*(xye(0,1) - xye(0,0));

  length = sqrt(norm[0]*norm[0]+norm[1]*norm[1]);

  norm[0] = (1.0/length)*norm[0];
  norm[1] = (1.0/length)*norm[1];

  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid][0];
    // get shape functions and derivatives in the line
    shape_function_1D(funct,e1,distype);
    shape_function_1D_deriv1(deriv,e1,distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye,deriv,iel);

    // values are multiplied by the product from inf. area element and the gauss weight
    const double fac = dr * intpoints.qwgt[gpid];

    for (int node=0;node<iel;++node)
    {
       for(int dim=0;dim<2;dim++)
      {
        elevec1[node*numdf+dim]+=funct[node] * fac * norm[dim];
      }
       elevec1[node*numdf+2] = 0.0;
    }

  } //end of loop over integration points
} // DRT::ELEMENTS::Fluid2Line::ElementNodeNormal

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Line::ElementSurfaceTension(ParameterList& params,
                                                      DRT::Discretization& discretization,
                                                      vector<int>& lm,
                                                      Epetra_SerialDenseVector& elevec1,
                                                      const std::vector<double>& edispnp,
                                                      std::vector<double>& enormals)
{
  // there are 2 velocities and 1 pressure
  const int numdf = 3;

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
    MAT::NewtonianFluid* actmat = static_cast<MAT::NewtonianFluid*>(mat.get());
    SFgamma = actmat->Gamma();
  }
  else if (mat->MaterialType()==INPAR::MAT::m_sutherland_fluid)
    //MAT::SutherlandFluid* actmat = static_cast<MAT::SutherlandFluid*>(mat.get());
    dserror("no gamma for Sutherland's fluid");
  else
    dserror("newtonian or sutherland fluid material expected but got type %d", mat->MaterialType());

  // gaussian points
  const DiscretizationType distype = this->Shape();
  const GaussRule1D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints1D  intpoints(gaussrule);

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
      xye(0,i) += edispnp[numdf*i];
      xye(1,i) += edispnp[numdf*i+1];
    }
  }

#ifdef FSTENS1
  //set normal vectors to length = 1.0
  for (int node=0;node<iel;++node)
  {
    double length = 0.0;
    for (int dim=0;dim<2;dim++)
    {
      length += enormals[numdf*node+dim] * enormals[numdf*node+dim];
    }
    length = sqrt(length);
    for (int dim=0;dim<2;dim++)
    {
      enormals[numdf*node+dim] = (1.0/length) * enormals[numdf*node+dim];
    }
  }
#endif

  // loop over integration points
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {
    const double e1 = intpoints.qxg[gpid][0];
    // get shape functions and derivatives in the line
    shape_function_1D(funct,e1,distype);
    shape_function_1D_deriv1(deriv,e1,distype);

    // compute infinitesimal line element dr for integration along the line
    const double dr = f2_substitution(xye,deriv,iel);

    // values are multiplied by the product from inf. area element and the
    // the gauss weight
    const double fac = intpoints.qwgt[gpid] *dr;

#ifdef FSTENS1

    // Metric-tensor for surface element (a line) is 1x1
    double metrictensor = 0.0;
    Epetra_SerialDenseVector  dxyds(2);
    for (int node=0;node<iel;++node)
    {
      for (int dim=0;dim<2;dim++)
        dxyds[dim] += xye(dim,node)*deriv(0,node);
    }
    for (int dim=0;dim<2;dim++)
      metrictensor += dxyds[dim]*dxyds[dim];

    // calculate normal vector at integration point
    Epetra_SerialDenseVector  norm(2);
    for (int dim=0;dim<2;dim++)
    {
      for (int node=0;node<iel;++node)
      {
        norm[dim] += funct[node] * enormals[numdf*node+dim];
      }
    }
    // set length to 1.0
    double length = 0.0;
    for (int dim=0;dim<2;dim++)
    {
      length += norm[dim] * norm[dim];
    }
    length = sqrt(length);
    for (int dim=0;dim<2;dim++)
    {
      norm[dim] = (1.0/length) * norm[dim];
    }

    // calculate double mean curvature 2*H at integration point.
    double twoH = 0.0;
    Epetra_SerialDenseVector dn12dr(2);

    for(int dim=0;dim<2;dim++)
    {
      for (int node=0;node<iel;++node)
      {
        dn12dr[dim] += enormals[numdf*node+dim] * deriv(0,node);
      }
    }

    //Acc. to Saksono eq. (4): 2H = - Surface_gradient*norm
    //see also Bronstein ..."mittlere Kruemmung" (2D: with H instead of 2H,
    //E=G, F=0,N=L)
    twoH = (-1.0) * dn12dr[0] * dxyds[0] + (-1.0) * dn12dr[1] * dxyds[1];


    for (int node=0;node<iel;++node)
        {
       for(int dim=0;dim<2;dim++)
        {
          // according to Saksono (23)
          elevec1[node*numdf+dim]+= SFgamma *
                                    twoH * norm[dim] * funct[node]
                                    * fac;
        }
        }

#else //----> FSTENS2
      // works fine

    Epetra_SerialDenseVector  dxyds(2);
  for (int node=0;node<iel;++node)
  {
      for (int dim=0;dim<2;dim++)
        dxyds[dim] += xye(dim,node)*deriv(0,node);
    }

    for (int node=0;node<iel;++node)
      {
       for(int dim=0;dim<2;dim++)
    {
          // Right hand side Integral (SFgamma * -Surface_Gradient, weighting
          // function) on Gamma_FS
          // See Saksono eq. (26)
          // discretized as surface gradient * ( Shapefunction-Matrix
          // transformed )
          elevec1[node*numdf+dim]+= SFgamma *
                                    (-1.0) * deriv(0, node) * dxyds[dim]
                                    * fac;
    }
  }
#endif

  } //end of loop over integration points

} // DRT::ELEMENTS::Fluid2Line::ElementSurfaceTension

#endif  // #ifdef CCADISCRET
#endif // #ifdef D_FLUID2
