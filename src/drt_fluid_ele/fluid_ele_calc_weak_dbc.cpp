
#include "fluid_ele_calc_weak_dbc.H"

#include "fluid_ele.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

// TODO: remove after Nurbs functions are changed
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"


//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                        INTERFACE CLASS
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//   Allocate one static instance of the internal implementation
//   class for weak dirichlet condition and return pointer to it
//-----------------------------------------------------------------
DRT::ELEMENTS::FluidBoundaryWeakDBCInterface* DRT::ELEMENTS::FluidBoundaryWeakDBCInterface::Impl(
  DRT::ELEMENTS::FluidBoundary* f3bdry
  )
{
  switch (f3bdry->Shape())
  {
  // 3D:
  case DRT::Element::quad4:
  {
    static FluidSurfaceWeakDBC<DRT::Element::quad4,DRT::Element::hex8>* fsurfq4;

    if(f3bdry->ParentElement()->Shape()==DRT::Element::hex8)
    {
      if (fsurfq4==NULL)
        fsurfq4 = new FluidSurfaceWeakDBC<DRT::Element::quad4,DRT::Element::hex8>();
    }
    else
    {
      dserror("expected combination quad4/hex8 for surface/parent pair");
    }
    return fsurfq4;
  }
  case DRT::Element::nurbs9:
  {
    static FluidSurfaceWeakDBC<DRT::Element::nurbs9,DRT::Element::nurbs27>* fsurfn9;

    if(f3bdry->ParentElement()->Shape()==DRT::Element::nurbs27)
    {
      if (fsurfn9==NULL)
        fsurfn9 = new FluidSurfaceWeakDBC<DRT::Element::nurbs9,DRT::Element::nurbs27>();
    }
    else
    {
      dserror("expected combination quad4/hex8 for surface/parent pair");
    }
    return fsurfn9;
  }
  // 2D:
  case DRT::Element::line2:
  {
    static FluidLineWeakDBC<DRT::Element::line2,DRT::Element::quad4>* fline2;

    if(f3bdry->ParentElement()->Shape()==DRT::Element::quad4)
    {
      if (fline2==NULL)
      {
        fline2 = new FluidLineWeakDBC<DRT::Element::line2,DRT::Element::quad4>();
      }
    }
    else
    {
      dserror("expected combination line2/quad4 for line/parent pair");
    }

    return fline2;
  }
  case DRT::Element::nurbs3:
  {
    static FluidLineWeakDBC<DRT::Element::nurbs3,DRT::Element::nurbs9>* flinen3;

    if(f3bdry->ParentElement()->Shape()==DRT::Element::nurbs9)
    {
      if (flinen3==NULL)
      {
        flinen3 = new FluidLineWeakDBC<DRT::Element::nurbs3,DRT::Element::nurbs9>();
      }
    }
    else
    {
      dserror("expected combination nurbs3/nurbs9 for line/parent pair");
    }

    return flinen3;
  }
  default:
    dserror("shape %d (%d nodes) not supported by weak DBC", f3bdry->Shape(), f3bdry->NumNode());
  }

  return NULL;
}


//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                        IMPLEMENTATION
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//                       empty constructor
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype>
DRT::ELEMENTS::FluidSurfaceWeakDBC<distype,pdistype>::FluidSurfaceWeakDBC()
{
  // pointer to class FluidImplParameter (access to the general parameter)
  fldpara_ = DRT::ELEMENTS::FluidEleParameter::Instance();

  return;
}

//-----------------------------------------------------------------
//             evaluate implementation for weak dbcs
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype>
int DRT::ELEMENTS::FluidSurfaceWeakDBC<distype,pdistype>::EvaluateWeakDBC(
  FluidBoundary*             surfele       ,
  ParameterList&             params        ,
  DRT::Discretization&       discretization,
  vector<int>&               lm            ,
  Epetra_SerialDenseMatrix&  elemat_epetra ,
  Epetra_SerialDenseVector&  elevec_epetra )
{

  //--------------------------------------------------
  // get the condition information
  RCP<DRT::Condition> wdbc_cond
    =
    params.get<RCP<DRT::Condition> >("condition");

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

  // initialise Spaldings law with parameters chi=0.4 and B=5.5
  FluidSurfaceWeakDBCSpaldingsLaw SpaldingsLaw(0.4,5.5);

  // decide whether to use it or not
  const string* deftauB
    =
    (*wdbc_cond).Get<string>("Definition of penalty parameter");

  bool spalding=false;

  if(*deftauB=="Spalding")
  {
    spalding=true;
  }
  else if(*deftauB=="constant")
  {
    spalding=false;
  }
  else
  {
    dserror("Unknown definition of penalty parameter tauB: %s",(*deftauB).c_str());
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
  // apply time curve at (n+1) for all time integration schemes (dirichlet condition)
  // time_ in fluid3_parameter is time(n+alphaF) in the case of genalpha
  // therefore, it need to be reset to time(n+1)
  const double time = fldpara_->Time()+(1-fldpara_->AlphaF())*fldpara_->Dt();
  if (time<0.0) usetime = false;


  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = (*wdbc_cond).Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
  curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and switches from the condition
  // (assumed to be constant on element boundary)
  const vector<int>* functions = (*wdbc_cond).Get<vector<int> >   ("funct");

  // I hope we have a linear element.
  // Ciarlet PG. The Finite element method for elliptic
  // problems. Amsterdam: North-Holland; 1978.

  // Bazilevs Michler etal use 4.0 for quadratic nurbs as well
  // (in combination with a dynamic computation of tau, so lets
  // try it as well)
  if(surfele->ParentElement()->Shape()!=Fluid::hex8 && surfele->ParentElement()->Shape()!=Fluid::nurbs27)
  {
    dserror("Cb up to now only implemented for trilinear orb nurbs elements");
  }

  // find out whether to apply weak dbcs only in normal direction
  bool onlynormal = false;

  const string* active_components
    =
    (*wdbc_cond).Get<string>("Directions to apply weak dbc");

  if(*active_components=="all_directions")
  {
     onlynormal=false;
  }
  else if(*active_components=="only_in_normal_direction")
  {
     onlynormal=true;
  }
  else
  {
    dserror("Unknown definition of active components: %s",(*active_components).c_str());
  }

  // optional scaling of penalty parameter
  const double scaling
    =
    (*wdbc_cond).GetDouble("TauBscaling");

  if(spalding && fabs(scaling-1.0)>1e-9)
  {
    dserror("tauB computed according to Spaldings law. Do not apply scaling!=1.0 !!!\n");
  }

  const double Cb = 4.0*scaling;

  // get value for boundary condition
  const vector<double>* val = (*wdbc_cond).Get<vector<double> >("val");

  if(spalding)
  {
    for(int i=0;i<3;++i)
    {
      if((*val)[i]*(*val)[i]>1e-9)
      {
        dserror("Applying Spaldings law to a wall with nonzero velocity\n");
      }
    }
  }

  // get time integration parameter
  //const double afgdt = params.get<double>("afgdt");
  //const double gdt   = params.get<double>("gdt");

  const double timefac = fldpara_->TimeFac();
  const double timefacpre   = fldpara_->TimeFacPre();
  const double timefacrhs   = fldpara_->TimeFacRhs();

  //--------------------------------------------------
  // get parent elements location vector and ownerships

  // the vectors have been allocated outside in
  // EvaluateConditionUsingParentData
  RCP<vector<int> > plm
    =
    params.get<RCP<vector<int> > >("plm");
  RCP<vector<int> > plmowner
    =
    params.get<RCP<vector<int> > >("plmowner");
  RCP<vector<int> > plmstride
    =
    params.get<RCP<vector<int> > >("plmstride");

  surfele->ParentElement()->LocationVector(discretization,*plm,*plmowner,*plmstride);

  // --------------------------------------------------
  // Reshape element matrices and vectors and init to zero, construct views
  const int eledim = 4*piel;
  elemat_epetra.Shape(eledim,eledim);
  elevec_epetra.Size (eledim);
  //
  LINALG::Matrix<eledim,eledim> elemat(elemat_epetra.A(),true);
  LINALG::Matrix<eledim,     1> elevec(elevec_epetra.A(),true);


  // --------------------------------------------------
  // extract velocities from global distributed vectors

  // velocities (intermediate time step, n+alpha_F)
  RCP<const Epetra_Vector> velaf
    =
    discretization.GetState("velaf");
  if (velaf==null)
    dserror("Cannot get state vector 'velaf'");

  vector<double> mypvelaf((*plm).size());
  DRT::UTILS::ExtractMyValues(*velaf,mypvelaf,*plm);

  // velocities n+1
  vector<double> mypvelnp((*plm).size());

  if((fldpara_->TimeAlgo()==INPAR::FLUID::timeint_gen_alpha) or
      (fldpara_->TimeAlgo()==INPAR::FLUID::timeint_npgenalpha))
  {
    // velocities (intermediate time step, n+1)
    RCP<const Epetra_Vector> velnp
      =
      discretization.GetState("velnp");
    if (velnp==null)
      dserror("Cannot get state vector 'velnp'");

    DRT::UTILS::ExtractMyValues(*velnp,mypvelnp,*plm);
  }
  // mypvelnp = mypvelaf
  else
    DRT::UTILS::ExtractMyValues(*velaf,mypvelnp,*plm);

  vector<double> myedispnp ((lm  ).size());
  vector<double> mypedispnp((*plm).size());
  if (surfele->ParentElement()->IsAle())
  {
    // mesh displacements, new time step, n+1
    RCP<const Epetra_Vector> dispnp
      =
      discretization.GetState("dispnp");
    if (dispnp==null)
    {
      dserror("Cannot get state vector 'dispnp'");
    }

    DRT::UTILS::ExtractMyValues(*dispnp,myedispnp ,lm  );
    DRT::UTILS::ExtractMyValues(*dispnp,mypedispnp,*plm);
  }

  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------

  // extract intermediate velocities
  for(int i=0;i<piel;++i)
  {
    const int fi=4*i;

    pevelaf_(0,i) = mypvelaf[  fi];
    pevelaf_(1,i) = mypvelaf[1+fi];
    pevelaf_(2,i) = mypvelaf[2+fi];
  }

  // extract current velocities and pressure
  for(int i=0;i<piel;++i)
  {
    const int fi=4*i;

    pevelnp_(0,i) = mypvelnp[  fi];
    pevelnp_(1,i) = mypvelnp[1+fi];
    pevelnp_(2,i) = mypvelnp[2+fi];

    peprenp_(  i) = mypvelnp[3+fi];
  }

  if (surfele->ParentElement()->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      const int fi=4*i;

      pedispnp_(0,i) = mypedispnp[  fi];
      pedispnp_(1,i) = mypedispnp[1+fi];
      pedispnp_(2,i) = mypedispnp[2+fi];
    }

    for (int i=0;i<iel;++i)
    {
      const int fi=4*i;

      edispnp_(0,i) = myedispnp[  fi];
      edispnp_(1,i) = myedispnp[1+fi];
      edispnp_(2,i) = myedispnp[2+fi];
    }
  }

  // extract node coords
  for(int i=0;i<piel;++i)
  {
    pxyze_(0,i)=surfele->ParentElement()->Nodes()[i]->X()[0];
    pxyze_(1,i)=surfele->ParentElement()->Nodes()[i]->X()[1];
    pxyze_(2,i)=surfele->ParentElement()->Nodes()[i]->X()[2];
  }

  if (surfele->ParentElement()->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      pxyze_(0,i) += pedispnp_(0,i);
      pxyze_(1,i) += pedispnp_(1,i);
      pxyze_(2,i) += pedispnp_(2,i);
    }
  }

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = surfele->ParentElement()->Material();

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(mat->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    // we need the kinematic viscosity here
    visc = actmat->Viscosity()/actmat->Density();
    if (actmat->Density() != 1.0)
      dserror("density 1.0 expected: the density need to be included in the linearization terms");
  }
  else
  {
    dserror("up to now I expect a constant viscosity to inforce weak DBCs\n");
  }

  //--------------------------------------------------
  //          GET BOUNDARY ELEMENT DATA
  //--------------------------------------------------

  // local surface id
  int surfaceid =surfele->SurfaceNumber();

  // extract node coords
  for(int i=0;i<iel;++i)
  {
    xyze_(0,i)=surfele->Nodes()[i]->X()[0];
    xyze_(1,i)=surfele->Nodes()[i]->X()[1];
    xyze_(2,i)=surfele->Nodes()[i]->X()[2];
  }

  if (surfele->ParentElement()->IsAle())
  {
    for (int i=0;i<iel;++i)
    {
      xyze_(0,i) += edispnp_(0,i);
      xyze_(1,i) += edispnp_(1,i);
      xyze_(2,i) += edispnp_(2,i);
    }
  }

  //--------------------------------------------------
  // get gausspoints to integrate over boundary element

  // get gauss rule
  DRT::UTILS::GaussRule2D gaussrule=DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
  case DRT::Element::quad4:
  {
    gaussrule = DRT::UTILS::intrule_quad_4point;
    break;
  }
  case DRT::Element::nurbs9:
  {
    gaussrule = DRT::UTILS::intrule_quad_9point;
    break;
  }
  default:
    dserror("invalid discretization type for fluid3surface weak DBC evaluation");
  }

  // gaussian points on surface
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //--------------------------------------------------
  // the gausspoints above have to be mapped to the
  // parent element to be able to evaluate one sided
  // derivatives on the boundary
  //
  // in addition, get information on the orientation of the
  // outward normal

  Epetra_SerialDenseMatrix pqxg(intpoints.nquad,3);

  SurfaceGPToParentGP(pqxg     ,
                      intpoints,
                      pdistype ,
                      distype  ,
                      surfaceid);

  // --------------------------------------------------
  // Now do the nurbs specific stuff

  std::vector<Epetra_SerialDenseVector> mypknots(3);
  std::vector<Epetra_SerialDenseVector> myknots (2);

  Epetra_SerialDenseVector weights(iel);
  LINALG::Matrix<piel,1>   pweights;

  // orientation of outward normal
  double normalfac=0.0;

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(surfele->Shape()==Fluid::nurbs4 || surfele->Shape()==Fluid::nurbs9)
  {
    // --------------------------------------------------
    // get knotvector

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

    bool zero_sized_parent
      =knots->GetBoundaryEleAndParentKnots(mypknots              ,
                                           myknots               ,
                                           normalfac             ,
                                           surfele->ParentElement()->Id(),
                                           surfaceid             );

    if(zero_sized_parent)
    {
      return 0;
    }

    // --------------------------------------------------
    // get node weights for nurbs elements
    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (surfele->Nodes()[inode]);

      weights(inode) = cp->W();
    }

    // extract node coords
    for(int i=0;i<piel;++i)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (surfele->ParentElement()->Nodes()[i]);

      pweights(i) = cp->W();
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
    const double xi    = intpoints.qxg [iquad][0];
    const double eta   = intpoints.qxg [iquad][1];

    // gaussian point in parent elements local coordinates
    const double r     = pqxg(iquad,0);
    const double s     = pqxg(iquad,1);
    const double t     = pqxg(iquad,2);

    if(!(distype == DRT::Element::nurbs9))
    {
      // ------------------------------------------------
      // shape function derivs of boundary element at gausspoint
      DRT::UTILS::shape_function_2D       (funct_,xi,eta,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv_,xi,eta,distype);
    }
    else
    {
      // this is just a temporary work-around
      Epetra_SerialDenseVector gp(2);
      gp(0)=xi;
      gp(1)=eta;

      Epetra_SerialDenseVector tempfunct(iel);
      Epetra_SerialDenseMatrix tempderiv(2,iel);

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (tempfunct,
         tempderiv,
         gp       ,
         myknots  ,
         weights  ,
         distype  );

      for(int rr=0;rr<iel;++rr)
      {
        funct_(  rr)=tempfunct(  rr);
        deriv_(0,rr)=tempderiv(0,rr);
        deriv_(1,rr)=tempderiv(1,rr);
      }
    }

    // ------------------------------------------------
    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration

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
    dxyzdrs_.MultiplyNT(deriv_,xyze_);
    /*
      |
      |      +-           -+    +-            -+   +-            -+ T
      |      |             |    | dx   dy   dz |   | dx   dy   dz |
      |      |  g11   g12  |    | --   --   -- |   | --   --   -- |
      |      |             |    | dr   dr   dr |   | dr   dr   dr |
      |      |             |  = |              | * |              |
      |      |             |    | dx   dy   dz |   | dx   dy   dz |
      |      |  g21   g22  |    | --   --   -- |   | --   --   -- |
      |      |             |    | ds   ds   ds |   | ds   ds   ds |
      |      +-           -+    +-            -+   +-            -+
      |
      | the calculation of g21 is redundant since g21=g12
    */
    metrictensor_.MultiplyNT(dxyzdrs_,dxyzdrs_);

    /*
                          +--------------+
                         /               |
           sqrtdetg =   /  g11*g22-g12^2
                      \/
    */

    drs_= sqrt(metrictensor_(0,0)*metrictensor_(1,1)
               -
               metrictensor_(0,1)*metrictensor_(1,0));


    // total integration factor
    const double fac = drs_*wquad;

    // ------------------------------------------------
    // compute normal
    if(distype!=DRT::Element::nurbs9)
    {
      double length = 0.0;
      n_(0) = (xyze_(1,1)-xyze_(1,0))*(xyze_(2,2)-xyze_(2,0))
        -
        (xyze_(2,1)-xyze_(2,0))*(xyze_(1,2)-xyze_(1,0));
      n_(1) = (xyze_(2,1)-xyze_(2,0))*(xyze_(0,2)-xyze_(0,0))
        -
        (xyze_(0,1)-xyze_(0,0))*(xyze_(2,2)-xyze_(2,0));
      n_(2) = (xyze_(0,1)-xyze_(0,0))*(xyze_(1,2)-xyze_(1,0))
        -
        (xyze_(1,1)-xyze_(1,0))*(xyze_(0,2)-xyze_(0,0));

      length = n_.Norm2();

      for(int i=0;i<3;++i)
      {
        n_(i)/=length;
      }
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
      n_(0) = dxyzdrs_(0,1)*dxyzdrs_(1,2)-dxyzdrs_(1,1)*dxyzdrs_(0,2);
      n_(1) = dxyzdrs_(0,2)*dxyzdrs_(1,0)-dxyzdrs_(1,2)*dxyzdrs_(0,0);
      n_(2) = dxyzdrs_(0,0)*dxyzdrs_(1,1)-dxyzdrs_(1,0)*dxyzdrs_(0,1);

      const double length = n_.Norm2()*normalfac;

      for(int i=0;i<3;++i)
      {
        n_(i)/=length;
      }
    }

    // ------------------------------------------------
    // factor given by spatial function
    LINALG::Matrix<3,1> functionfac;
    for(int i=0;i<3;++i)
    {
      functionfac(i)= 1.0;
    }

    // determine coordinates of current Gauss point
    LINALG::Matrix<3,1> coordgp;
    for(int i=0;i<3;++i)
    {
      coordgp(i)= 0.0;
    }
    for (int i=0;i<iel;++i)
    {
      for(int j=0;j<3;++j)
      {
        coordgp(j)+=xyze_(j,i)*funct_(i);
      }
    }

    int functnum = -1;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<3;++dim)
      {
	// factor given by spatial function
	if (functions)
	{
	  functnum = (*functions)[dim];
	  if (functnum>0)
	  {
	    // evaluate function at current gauss point (important: requires 3D position vector)
	    functionfac(dim) = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgp.A(),time,NULL);
	  }
	  else
	  {
	    functionfac(dim) = 1.0;
	  }
	}
      }
    }

    // ------------------------------------------------
    // shape functions and derivs of corresponding parent at gausspoint
    if(!(pdistype == DRT::Element::nurbs27))
    {
      DRT::UTILS::shape_function_3D       (pfunct_,r,s,t,pdistype);
      DRT::UTILS::shape_function_3D_deriv1(pderiv_,r,s,t,pdistype);
    }
    else
    {
      // set gauss point coordinates
      LINALG::Matrix<3,1> gp;

      gp(0)=r;
      gp(1)=s;
      gp(2)=t;

      DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
        (pfunct_ ,
         pderiv_ ,
         gp      ,
         mypknots,
         pweights,
         pdistype);
    }

    //-----------------------------------------------------
    //
    //                       +-------------+
    //                      / /  T       \ |
    //           h = 2 * \ / |  n * G * n |
    //            b       +   \          /
    //

    // get Jacobian matrix and determinant
    pxjm_=0;

    for(int i=0;i<piel;++i)
    {
      for(int rr=0;rr<3;++rr)
      {
	for(int mm=0;mm<3;++mm)
	{
	  pxjm_(rr,mm)+=pderiv_(rr,i)*pxyze_(mm,i);
	}
      }
    }

    const double pdet =
      pxjm_(0,0)*pxjm_(1,1)*pxjm_(2,2)+
      pxjm_(0,1)*pxjm_(1,2)*pxjm_(2,0)+
      pxjm_(0,2)*pxjm_(1,0)*pxjm_(2,1)-
      pxjm_(0,2)*pxjm_(1,1)*pxjm_(2,0)-
      pxjm_(0,0)*pxjm_(1,2)*pxjm_(2,1)-
      pxjm_(0,1)*pxjm_(1,0)*pxjm_(2,2);

    // check for degenerated elements
    if (pdet < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",
              surfele->ParentElement()->Id(),
              pdet);
    }

    //-----------------------------------------------------
    //
    //             compute global first derivates
    //
    /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

    */

    // inverse of jacobian (transposed)
    /*
          +-                 -+     +-                 -+ -1
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dx    dx    dx   |     |  dr    dr    dr   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |  =  |  --    --    --   |
          |  dy    dy    dy   |     |  ds    ds    ds   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dz    dz    dz   |     |  dt    dt    dt   |
          +-                 -+     +-                 -+

    */
    pxji_(0,0) = (  pxjm_(1,1)*pxjm_(2,2) - pxjm_(2,1)*pxjm_(1,2))/pdet;
    pxji_(1,0) = (- pxjm_(1,0)*pxjm_(2,2) + pxjm_(2,0)*pxjm_(1,2))/pdet;
    pxji_(2,0) = (  pxjm_(1,0)*pxjm_(2,1) - pxjm_(2,0)*pxjm_(1,1))/pdet;
    pxji_(0,1) = (- pxjm_(0,1)*pxjm_(2,2) + pxjm_(2,1)*pxjm_(0,2))/pdet;
    pxji_(1,1) = (  pxjm_(0,0)*pxjm_(2,2) - pxjm_(2,0)*pxjm_(0,2))/pdet;
    pxji_(2,1) = (- pxjm_(0,0)*pxjm_(2,1) + pxjm_(2,0)*pxjm_(0,1))/pdet;
    pxji_(0,2) = (  pxjm_(0,1)*pxjm_(1,2) - pxjm_(1,1)*pxjm_(0,2))/pdet;
    pxji_(1,2) = (- pxjm_(0,0)*pxjm_(1,2) + pxjm_(1,0)*pxjm_(0,2))/pdet;
    pxji_(2,2) = (  pxjm_(0,0)*pxjm_(1,1) - pxjm_(1,0)*pxjm_(0,1))/pdet;


    //-----------------------------------------------------
    /*            +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
    */
    LINALG::Matrix<3,3> G;

    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        G(nn,rr) = pxji_(nn,0)*pxji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          G(nn,rr) += pxji_(nn,mm)*pxji_(rr,mm);
        }
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

    double nGn=0;
    for (int nn=0;nn<3;++nn)
    {
      for (int rr=0;rr<3;++rr)
      {
        nGn+=n_(rr)*G(rr,nn)*n_(nn);
      }
    }
    const double h =2.0/sqrt(nGn);

    //-----------------------------------------------------
    // compute global derivates at integration point
    //
    //   dN    +-----  dN (xi)    dxi
    //     i    \        i           k
    //   --- =   +     ------- * -----
    //   dx     /        dxi      dx
    //     j   +-----       k       j
    //         node k
    //
    // j : direction of derivative x/y/z
    //
    for(int nn=0;nn<piel;++nn)
    {
      for(int rr=0;rr<3;++rr)
      {
        pderxy_(rr,nn)=pxji_(rr,0)*pderiv_(0,nn);

        for(int mm=1;mm<3;++mm)
        {
          pderxy_(rr,nn)+=pxji_(rr,mm)*pderiv_(mm,nn);
        }
      }
    }

    //-----------------------------------------------------
    // get velocities (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    vel   (x) =   +      N (x) * vel
    //                 /        j         j
    //                +-----
    //                node j
    //
    for(int rr=0;rr<3;++rr)
    {
      velintnp_(rr)=pfunct_(0)*pevelnp_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        velintnp_(rr)+=pfunct_(nn)*pevelnp_(rr,nn);
      }
    }

    //-----------------------------------------------------
    // get pressure (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    pre   (x) =   +      N (x) * pre
    //                 /        i         i
    //                +-----
    //                node i
    //
    prenp_=pfunct_(0)*peprenp_(0);
    for(int nn=1;nn<piel;++nn)
    {
      prenp_+=pfunct_(nn)*peprenp_(nn);
    }

    //-----------------------------------------------------
    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----
    //       n+af       \                  n+af
    //    vel    (x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    for(int rr=0;rr<3;++rr)
    {
      velintaf_(rr)=pfunct_(0)*pevelaf_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        velintaf_(rr)+=pfunct_(nn)*pevelaf_(rr,nn);
      }
    }

    //-----------------------------------------------------
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
    for(int rr=0;rr<3;++rr)
    {
      for(int mm=0;mm<3;++mm)
      {
        vderxyaf_(rr,mm)=pderxy_(mm,0)*pevelaf_(rr,0);
        for(int nn=1;nn<piel;++nn)
        {
          vderxyaf_(rr,mm)+=pderxy_(mm,nn)*pevelaf_(rr,nn);
        }
      }
    }

    /*
    // ---------------------------------------------------
    // define (initial) penalty parameter
    //
    //
    //                         /    \
    //                        |  nu  |
    //         tau      = C * | ---- |
    //            B(,0)    b  |  h   |
    //                         \  b /
    */
    double tau_B = Cb*visc/h;

    // ---------------------------------------------------
    // update penalty parameter for Spalding law
    if(spalding)
    {
      //                             +--------------------+
      //                            /  +---+              |
      //       ||  n+af ||         /    \    n+af    n+af
      //       || u     ||  =     /      +  u     * u
      //                         /      /    j       j
      //                      \ /      +---+
      //                       +       dim j

      const double normu = velintaf_.Norm2();

      SpaldingsLaw.ComputeTauBUsingSpaldingsLaw(tau_B,
                                                normu,
                                                h    ,
                                                Cb   ,
                                                visc );
    }

    // ---------------------------------------------------
    // velocities on boundary
    //
    //         n+af
    //        u    - u
    //                b
    //

    LINALG::Matrix<3,1> bvres;

    bvres(0)=(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);
    bvres(1)=(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);
    bvres(2)=(velintaf_(2)-(*val)[2]*functionfac(2)*curvefac);


    // ---------------------------------------------------
    // momentum flux over boundary

    const double flux=velintaf_(0)*n_(0)+velintaf_(1)*n_(1)+velintaf_(2)*n_(2);

    //--------------------------------------------------
    // partially integrated pressure term, rescaled by gamma*dt
    /*
    // factor: 1.0
    //
    //             /            \
    //            |              |
    //          + |  v , Dp * n  |
    //            |              |
    //             \            / boundaryele
    //
    */
    for (int ui=0; ui<piel; ++ui)
    {
      for (int vi=0; vi<piel; ++vi)
      {
        elemat(vi*4  ,ui*4+3) += fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*4+1,ui*4+3) += fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(1);
        elemat(vi*4+2,ui*4+3) += fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(2);
      }
    }

    for (int vi=0; vi<piel; ++vi)
    {
      elevec(vi*4    ) -= fac*timefacrhs*pfunct_(vi)*n_(0)*prenp_;
      elevec(vi*4 + 1) -= fac*timefacrhs*pfunct_(vi)*n_(1)*prenp_;
      elevec(vi*4 + 2) -= fac*timefacrhs*pfunct_(vi)*n_(2)*prenp_;
    }

    //--------------------------------------------------
    // adjoint consistency term, pressure/continuity part
    /*
    // factor: gdt
    //
    //             /              \
    //            |                |
    //          - |  q , Dacc * n  |
    //            |                |
    //             \              / boundaryele
    //
    */
    for (int ui=0; ui<piel; ++ui)
    {
      for (int vi=0; vi<piel; ++vi)
      {
        elemat(vi*4+3,ui*4  ) -= fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*4+3,ui*4+1) -= fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(1);
        elemat(vi*4+3,ui*4+2) -= fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(2);
      }
    }

    /*
    // factor: 1.0
    //
    //             /                       \
    //            |       / n+1     \       |
    //          + |  q , | u   - u   | * n  |
    //            |       \ (i)   B /       |
    //             \                       / boundaryele
    //
    */
    for (int vi=0; vi<piel; ++vi)
    {
      elevec(vi*4+3) += fac*timefacrhs*pfunct_(vi)*
        ((velintnp_(0)-(*val)[0]*functionfac(0)*curvefac)*n_(0)
         +
         (velintnp_(1)-(*val)[1]*functionfac(1)*curvefac)*n_(1)
         +
         (velintnp_(2)-(*val)[2]*functionfac(2)*curvefac)*n_(2));
    }

    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    //           weak boundary conditions in all directions
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    if(!onlynormal)
    {
      const double timefacnu=fac*2.0*visc*timefac;

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
        double nabla_u_o_n_lin[3][3];

        nabla_u_o_n_lin[0][0]=timefacnu*(    pderxy_(0,ui)*n_(0)
                                         +
                                         0.5*pderxy_(1,ui)*n_(1)
                                         +
                                         0.5*pderxy_(2,ui)*n_(2));
        nabla_u_o_n_lin[0][1]=timefacnu*(0.5*pderxy_(0,ui)*n_(1));
        nabla_u_o_n_lin[0][2]=timefacnu*(0.5*pderxy_(0,ui)*n_(2));

        nabla_u_o_n_lin[1][0]=timefacnu*(0.5*pderxy_(1,ui)*n_(0));
        nabla_u_o_n_lin[1][1]=timefacnu*(0.5*pderxy_(0,ui)*n_(0)
                                         +
                                         pderxy_(1,ui)*n_(1)
                                         +
                                         0.5*pderxy_(2,ui)*n_(2));
        nabla_u_o_n_lin[1][2]=timefacnu*(0.5*pderxy_(1,ui)*n_(2));

        nabla_u_o_n_lin[2][0]=timefacnu*(0.5*pderxy_(2,ui)*n_(0));
        nabla_u_o_n_lin[2][1]=timefacnu*(0.5*pderxy_(2,ui)*n_(1));
        nabla_u_o_n_lin[2][2]=timefacnu*(0.5*pderxy_(0,ui)*n_(0)
                                         +
                                         0.5*pderxy_(1,ui)*n_(1)
                                         +
                                         pderxy_(2,ui)*n_(2));


        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4    ,ui*4    ) -= pfunct_(vi)*nabla_u_o_n_lin[0][0];
          elemat(vi*4    ,ui*4 + 1) -= pfunct_(vi)*nabla_u_o_n_lin[0][1];
          elemat(vi*4    ,ui*4 + 2) -= pfunct_(vi)*nabla_u_o_n_lin[0][2];

          elemat(vi*4 + 1,ui*4    ) -= pfunct_(vi)*nabla_u_o_n_lin[1][0];
          elemat(vi*4 + 1,ui*4 + 1) -= pfunct_(vi)*nabla_u_o_n_lin[1][1];
          elemat(vi*4 + 1,ui*4 + 2) -= pfunct_(vi)*nabla_u_o_n_lin[1][2];

          elemat(vi*4 + 2,ui*4    ) -= pfunct_(vi)*nabla_u_o_n_lin[2][0];
          elemat(vi*4 + 2,ui*4 + 1) -= pfunct_(vi)*nabla_u_o_n_lin[2][1];
          elemat(vi*4 + 2,ui*4 + 2) -= pfunct_(vi)*nabla_u_o_n_lin[2][2];
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

      {
        double nabla_u_o_n[3];
        nabla_u_o_n[0]=fac*timefacrhs*2.0*visc*
          (                     vderxyaf_(0,0) *n_(0)
           +0.5*(vderxyaf_(0,1)+vderxyaf_(1,0))*n_(1)
           +0.5*(vderxyaf_(0,2)+vderxyaf_(2,0))*n_(2));
        nabla_u_o_n[1]=fac*2.0*visc*
          ( 0.5*(vderxyaf_(1,0)+vderxyaf_(0,1))*n_(0)
            +                    vderxyaf_(1,1) *n_(1)
            +0.5*(vderxyaf_(1,2)+vderxyaf_(2,1))*n_(2));
        nabla_u_o_n[2]=fac*2.0*visc*
          ( 0.5*(vderxyaf_(2,0)+vderxyaf_(0,2))*n_(0)
           +0.5*(vderxyaf_(2,1)+vderxyaf_(1,2))*n_(1)
           +                    vderxyaf_(2,2) *n_(2));

        for (int vi=0; vi<piel; ++vi)
        {
          elevec(vi*4    ) += pfunct_(vi)*nabla_u_o_n[0];
          elevec(vi*4 + 1) += pfunct_(vi)*nabla_u_o_n[1];
          elevec(vi*4 + 2) += pfunct_(vi)*nabla_u_o_n[2];
        }
      }

      //--------------------------------------------------
      // (adjoint) consistency term, viscous part

      const double consistencytimefac=fac*2.0*visc*wd_gamma*timefac;


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

      LINALG::Matrix<3,3> nabla_s_w_o_n;
      for (int vi=0; vi<piel; ++vi)
      {
        nabla_s_w_o_n(0,0)=consistencytimefac*(n_(0)*(    pderxy_(0,vi))+n_(1)*(0.5*pderxy_(1,vi))+n_(2)*(0.5*pderxy_(2,vi)));
        nabla_s_w_o_n(0,1)=consistencytimefac*(n_(0)*(0.5*pderxy_(1,vi)));
        nabla_s_w_o_n(0,2)=consistencytimefac*(n_(0)*(0.5*pderxy_(2,vi)));

        nabla_s_w_o_n(1,0)=consistencytimefac*(n_(1)*(0.5*pderxy_(0,vi)));
        nabla_s_w_o_n(1,1)=consistencytimefac*(n_(0)*(0.5*pderxy_(0,vi))+n_(1)*(    pderxy_(1,vi))+n_(2)*(0.5*pderxy_(2,vi)));
        nabla_s_w_o_n(1,2)=consistencytimefac*(n_(1)*(0.5*pderxy_(2,vi)));

        nabla_s_w_o_n(2,0)=consistencytimefac*(n_(2)*(0.5*pderxy_(0,vi)));
        nabla_s_w_o_n(2,1)=consistencytimefac*(n_(2)*(0.5*pderxy_(1,vi)));
        nabla_s_w_o_n(2,2)=consistencytimefac*(n_(0)*(0.5*pderxy_(0,vi))+n_(1)*(0.5*pderxy_(1,vi))+n_(2)*(    pderxy_(2,vi)));


        for (int ui=0; ui<piel; ++ui)
        {
          elemat(vi*4    ,ui*4    ) -= pfunct_(ui)*nabla_s_w_o_n(0,0);
          elemat(vi*4    ,ui*4 + 1) -= pfunct_(ui)*nabla_s_w_o_n(0,1);
          elemat(vi*4    ,ui*4 + 2) -= pfunct_(ui)*nabla_s_w_o_n(0,2);

          elemat(vi*4 + 1,ui*4    ) -= pfunct_(ui)*nabla_s_w_o_n(1,0);
          elemat(vi*4 + 1,ui*4 + 1) -= pfunct_(ui)*nabla_s_w_o_n(1,1);
          elemat(vi*4 + 1,ui*4 + 2) -= pfunct_(ui)*nabla_s_w_o_n(1,2);

          elemat(vi*4 + 2,ui*4    ) -= pfunct_(ui)*nabla_s_w_o_n(2,0);
          elemat(vi*4 + 2,ui*4 + 1) -= pfunct_(ui)*nabla_s_w_o_n(2,1);
          elemat(vi*4 + 2,ui*4 + 2) -= pfunct_(ui)*nabla_s_w_o_n(2,2);
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
        elevec(vi*4    ) += fac*timefacrhs*2.0*visc*wd_gamma*(
          n_(0)*bvres(0)*(    pderxy_(0,vi))
          +
          n_(1)*bvres(0)*(0.5*pderxy_(1,vi))
          +
          n_(2)*bvres(0)*(0.5*pderxy_(2,vi))
          +
          n_(0)*bvres(1)*(0.5*pderxy_(1,vi))
          +
          n_(0)*bvres(2)*(0.5*pderxy_(2,vi)));

        elevec(vi*4 + 1) += fac*timefacrhs*2.0*visc*wd_gamma*(
          n_(1)*bvres(0)*(0.5*pderxy_(0,vi))
          +
          n_(0)*bvres(1)*(0.5*pderxy_(0,vi))
          +
          n_(1)*bvres(1)*(    pderxy_(1,vi))
          +
          n_(2)*bvres(1)*(0.5*pderxy_(2,vi))
          +
          n_(1)*bvres(2)*(0.5*pderxy_(2,vi)));

        elevec(vi*4 + 2) += fac*timefacrhs*2.0*visc*wd_gamma*(
          n_(2)*bvres(0)*(0.5*pderxy_(0,vi))
          +
          n_(2)*bvres(1)*(0.5*pderxy_(1,vi))
          +
          n_(0)*bvres(2)*(0.5*pderxy_(0,vi))
          +
          n_(1)*bvres(2)*(0.5*pderxy_(1,vi))
          +
          n_(2)*bvres(2)*(    pderxy_(2,vi)));
      }

      //--------------------------------------------------
      // adjoint consistency term, convective part (only on inflow)

      if(flux<0)
      {
        if(complete_linearisation)
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
              elemat(vi*4    ,ui*4    ) -= fac*timefac*pfunct_(vi)*n_(0)*bvres(0);
              elemat(vi*4    ,ui*4 + 1) -= fac*timefac*pfunct_(vi)*n_(1)*bvres(0);
              elemat(vi*4    ,ui*4 + 2) -= fac*timefac*pfunct_(vi)*n_(2)*bvres(0);

              elemat(vi*4 + 1,ui*4    ) -= fac*timefac*pfunct_(vi)*n_(0)*bvres(1);
              elemat(vi*4 + 1,ui*4 + 1) -= fac*timefac*pfunct_(vi)*n_(1)*bvres(1);
              elemat(vi*4 + 1,ui*4 + 2) -= fac*timefac*pfunct_(vi)*n_(2)*bvres(1);

              elemat(vi*4 + 2,ui*4    ) -= fac*timefac*pfunct_(vi)*n_(0)*bvres(2);
              elemat(vi*4 + 2,ui*4 + 1) -= fac*timefac*pfunct_(vi)*n_(1)*bvres(2);
              elemat(vi*4 + 2,ui*4 + 2) -= fac*timefac*pfunct_(vi)*n_(2)*bvres(2);
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

          const double fluxtimefac =fac*timefac*flux;

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*4    ,ui*4    ) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
              elemat(vi*4 + 1,ui*4 + 1) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
              elemat(vi*4 + 2,ui*4 + 2) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
            }
          }
        } // end if full_linearisation

        const double fluxfac =fac*flux*timefacrhs;

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
          elevec(vi*4    ) += fluxfac*pfunct_(vi)*bvres(0);
          elevec(vi*4 + 1) += fluxfac*pfunct_(vi)*bvres(1);
          elevec(vi*4 + 2) += fluxfac*pfunct_(vi)*bvres(2);
        }
      } // end if flux<0, i.e. boundary is an inflow boundary

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

      const double penaltytimefac=timefac*tau_B*fac;

      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          const double temp=penaltytimefac*pfunct_(ui)*pfunct_(vi);

          elemat(vi*4    ,ui*4    ) += temp;
          elemat(vi*4 + 1,ui*4 + 1) += temp;
          elemat(vi*4 + 2,ui*4 + 2) += temp;
        }
      }

      const double penaltyfac=tau_B*fac*timefacrhs;

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
        elevec(vi*4    ) -= penaltyfac*pfunct_(vi)*bvres(0);
        elevec(vi*4 + 1) -= penaltyfac*pfunct_(vi)*bvres(1);
        elevec(vi*4 + 2) -= penaltyfac*pfunct_(vi)*bvres(2);
      }
    } // !onlynormal
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    //          weak boundary conditions only in normal direction
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    else
    {
      //--------------------------------------------------
      // partially integrated viscous term

      const double timefacnu=fac*2.0*visc*timefac;

      /*
      // factor: 2*nu
      //
      //    /                             \
      //   |                   s           |
      // + |  v * n , n * nabla  Dacc * n  |
      //   |                               |
      //    \                             / boundaryele
      //
      */
      for (int ui=0; ui<piel; ++ui)
      {

        const double aux=timefacnu*(pderxy_(0,ui)*n_(0)+pderxy_(1,ui)*n_(1)+pderxy_(2,ui)*n_(2));
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4    ,ui*4    ) -= pfunct_(vi)*n_(0)*n_(0)*aux;
          elemat(vi*4    ,ui*4 + 1) -= pfunct_(vi)*n_(0)*n_(1)*aux;
          elemat(vi*4    ,ui*4 + 2) -= pfunct_(vi)*n_(0)*n_(2)*aux;

          elemat(vi*4 + 1,ui*4    ) -= pfunct_(vi)*n_(1)*n_(0)*aux;
          elemat(vi*4 + 1,ui*4 + 1) -= pfunct_(vi)*n_(1)*n_(1)*aux;
          elemat(vi*4 + 1,ui*4 + 2) -= pfunct_(vi)*n_(1)*n_(2)*aux;

          elemat(vi*4 + 2,ui*4    ) -= pfunct_(vi)*n_(2)*n_(0)*aux;
          elemat(vi*4 + 2,ui*4 + 1) -= pfunct_(vi)*n_(2)*n_(1)*aux;
          elemat(vi*4 + 2,ui*4 + 2) -= pfunct_(vi)*n_(2)*n_(2)*aux;
        }
      }

      /*
      // factor: 2*nu
      //
      //    /                             \
      //   |                   s  n+af     |
      // + |  v * n , n * nabla  u    * n  |
      //   |                               |
      //    \                             / boundaryele
      //
      */
      double n_o_nabla_u_o_n =
        vderxyaf_(0,0)*n_(0)*n_(0)
        +
        vderxyaf_(1,1)*n_(1)*n_(1)
        +
        vderxyaf_(2,2)*n_(2)*n_(2)
        +
        (vderxyaf_(0,1)+vderxyaf_(1,0))*n_(0)*n_(1)
        +
        (vderxyaf_(0,2)+vderxyaf_(2,0))*n_(0)*n_(2)
        +
        (vderxyaf_(1,2)+vderxyaf_(2,1))*n_(2)*n_(1);

      for (int vi=0; vi<piel; ++vi)
      {
        elevec(vi*4    ) += fac*timefacrhs*2.0*visc*wd_gamma*pfunct_(vi)*n_(0)*n_o_nabla_u_o_n;
        elevec(vi*4 + 1) += fac*timefacrhs*2.0*visc*wd_gamma*pfunct_(vi)*n_(1)*n_o_nabla_u_o_n;
        elevec(vi*4 + 2) += fac*timefacrhs*2.0*visc*wd_gamma*pfunct_(vi)*n_(2)*n_o_nabla_u_o_n;
      }

      //--------------------------------------------------
      // (adjoint) consistency term, viscous part

      const double consistencytimefac=fac*2.0*visc*wd_gamma*timefac;

      /*
      // factor: 2*nu*gamma_wd
      //
      //    /                               \
      //   |            s                    |
      // + |   n * nabla  w * n ,  Dacc * n  |
      //   |                                 |
      //    \                               / boundaryele
      //
      */

      for (int vi=0; vi<piel; ++vi)
      {
        for (int ui=0; ui<piel; ++ui)
        {
          const double aux=pderxy_(0,vi)*n_(0)+pderxy_(1,vi)*n_(1)+pderxy_(2,vi)*n_(2);

          elemat(vi*4    ,ui*4    ) -= aux*n_(0)*n_(0)*pfunct_(ui)*consistencytimefac;
          elemat(vi*4    ,ui*4 + 1) -= aux*n_(0)*n_(1)*pfunct_(ui)*consistencytimefac;
          elemat(vi*4    ,ui*4 + 2) -= aux*n_(0)*n_(2)*pfunct_(ui)*consistencytimefac;

          elemat(vi*4 + 1,ui*4    ) -= aux*n_(1)*n_(0)*pfunct_(ui)*consistencytimefac;
          elemat(vi*4 + 1,ui*4 + 1) -= aux*n_(1)*n_(1)*pfunct_(ui)*consistencytimefac;
          elemat(vi*4 + 1,ui*4 + 2) -= aux*n_(1)*n_(2)*pfunct_(ui)*consistencytimefac;

          elemat(vi*4 + 2,ui*4    ) -= aux*n_(2)*n_(0)*pfunct_(ui)*consistencytimefac;
          elemat(vi*4 + 2,ui*4 + 1) -= aux*n_(2)*n_(1)*pfunct_(ui)*consistencytimefac;
          elemat(vi*4 + 2,ui*4 + 2) -= aux*n_(2)*n_(2)*pfunct_(ui)*consistencytimefac;
        }
      }
      /*
      // factor: 2*nu*gamma_wd
      //
      //    /                                        \
      //   |            s          / n+af     \       |
      // + |   n * nabla  w * n , | u    - u   | * n  |
      //   |                       \        b /       |
      //    \                                        / boundaryele
      //
      */

      double bvres_o_n=bvres(0)*n_(0)+bvres(1)*n_(1)+bvres(2)*n_(2);

      for (int vi=0; vi<piel; ++vi)
      {
        double aux=(pderxy_(0,vi)*n_(0)+pderxy_(1,vi)*n_(1)+pderxy_(2,vi)*n_(2));

        elevec(vi*4    ) += fac*timefacrhs*2.0*visc*wd_gamma*aux*n_(0)*bvres_o_n;
        elevec(vi*4 + 1) += fac*timefacrhs*2.0*visc*wd_gamma*aux*n_(1)*bvres_o_n;
        elevec(vi*4 + 2) += fac*timefacrhs*2.0*visc*wd_gamma*aux*n_(2)*bvres_o_n;
      }

      //--------------------------------------------------
      // penalty term

      const double penaltytimefac=timefac*tau_B*fac;

      /*
      // factor: nu*Cb/h*afgdt
      //
      //    /                  \
      //   |                    |
      // + |  w o n , Dacc o n  |
      //   |                    |
      //    \                  / boundaryele
      //
      */
      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*4    ,ui*4    ) += penaltytimefac*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
          elemat(vi*4    ,ui*4 + 1) += penaltytimefac*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);
          elemat(vi*4    ,ui*4 + 2) += penaltytimefac*pfunct_(vi)*n_(0)*n_(2)*pfunct_(ui);

          elemat(vi*4 + 1,ui*4    ) += penaltytimefac*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
          elemat(vi*4 + 1,ui*4 + 1) += penaltytimefac*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
          elemat(vi*4 + 1,ui*4 + 2) += penaltytimefac*pfunct_(vi)*n_(1)*n_(2)*pfunct_(ui);

          elemat(vi*4 + 2,ui*4    ) += penaltytimefac*pfunct_(vi)*n_(2)*n_(0)*pfunct_(ui);
          elemat(vi*4 + 2,ui*4 + 1) += penaltytimefac*pfunct_(vi)*n_(2)*n_(1)*pfunct_(ui);
          elemat(vi*4 + 2,ui*4 + 2) += penaltytimefac*pfunct_(vi)*n_(2)*n_(2)*pfunct_(ui);
        }
      }

      const double penaltyfac=tau_B*fac*timefacrhs;

      /*
      // factor: nu*Cb/h
      //
      //    /                           \
      //   |            / n+af   \       |
      // + |   w o n , | u    - u | o n  |
      //   |            \  b     /       |
      //    \                           / boundaryele
      //
      */

      for (int vi=0; vi<piel; ++vi)
      {
        elevec(vi*4    ) -= penaltyfac*pfunct_(vi)*n_(0)*bvres_o_n;
        elevec(vi*4 + 1) -= penaltyfac*pfunct_(vi)*n_(1)*bvres_o_n;
        elevec(vi*4 + 2) -= penaltyfac*pfunct_(vi)*n_(2)*bvres_o_n;
      }

      //--------------------------------------------------
      // adjoint consistency term, convective part (only on inflow)

      if(flux<0)
      {
        if(complete_linearisation)
        {
          /*
          // These linearisations have only to be included if
          // u*n is negative --- otherwise they're nonesense
          //
          // factor: afgdt
          //
          //    /                                              \
          //   |    /        \    /     \     / n+af     \      |
          // - |   | Dacc * n |  | w o n | , | u    - u   | o n |
          //   |    \        /    \     /     \        b /      |
          //    \                                              / boundaryele, inflow
          //
          */

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*4    ,ui*4    ) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
              elemat(vi*4    ,ui*4 + 1) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);
              elemat(vi*4    ,ui*4 + 2) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(2)*pfunct_(ui);

              elemat(vi*4 + 1,ui*4    ) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
              elemat(vi*4 + 1,ui*4 + 1) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
              elemat(vi*4 + 1,ui*4 + 2) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(2)*pfunct_(ui);

              elemat(vi*4 + 2,ui*4    ) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(2)*n_(0)*pfunct_(ui);
              elemat(vi*4 + 2,ui*4 + 1) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(2)*n_(1)*pfunct_(ui);
              elemat(vi*4 + 2,ui*4 + 2) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(2)*n_(2)*pfunct_(ui);
            }
          }

          const double fluxtimefac =fac*timefac*flux;

          /*
          // factor: afgdt
          //
          //    /                                   \
          //   |    / n+af   \   /     \             |
          // - |   | u    * n | | w o n | , Dacc o n |
          //   |    \        /   \     /             |
          //    \  |          |                     / boundaryele, inflow
          //       +----------+
          //           <0
          */
          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*4    ,ui*4    ) -= fluxtimefac*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
              elemat(vi*4    ,ui*4 + 1) -= fluxtimefac*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);
              elemat(vi*4    ,ui*4 + 2) -= fluxtimefac*pfunct_(vi)*n_(0)*n_(2)*pfunct_(ui);

              elemat(vi*4 + 1,ui*4    ) -= fluxtimefac*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
              elemat(vi*4 + 1,ui*4 + 1) -= fluxtimefac*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
              elemat(vi*4 + 1,ui*4 + 2) -= fluxtimefac*pfunct_(vi)*n_(1)*n_(2)*pfunct_(ui);

              elemat(vi*4 + 2,ui*4    ) -= fluxtimefac*pfunct_(vi)*n_(2)*n_(0)*pfunct_(ui);
              elemat(vi*4 + 2,ui*4 + 1) -= fluxtimefac*pfunct_(vi)*n_(2)*n_(1)*pfunct_(ui);
              elemat(vi*4 + 2,ui*4 + 2) -= fluxtimefac*pfunct_(vi)*n_(2)*n_(2)*pfunct_(ui);
            }
          }

        } // end complete_linearisation

        const double fluxfac =fac*flux*timefacrhs;

        /*
        // factor: 1
        //
        //    /                                            \
        //   |    / n+af   \   /     \    / n+af     \      |
        // - |   | u    * n | | w o n |, | u    - u   | o n |
        //   |    \        /   \     /    \        b /      |
        //    \  |          |                              / boundaryele, inflow
        //       +----------+
        //           <0
        */
        for (int vi=0; vi<piel; ++vi)
        {
          elevec(vi*4    ) += fluxfac*pfunct_(vi)*n_(0)*bvres_o_n;
          elevec(vi*4 + 1) += fluxfac*pfunct_(vi)*n_(1)*bvres_o_n;
          elevec(vi*4 + 2) += fluxfac*pfunct_(vi)*n_(2)*bvres_o_n;
        }
      } // end if flux<0, i.e. boundary is an inflow boundary
    } // onlynormal

  } // end gaussloop

  return 0;
}

//-----------------------------------------------------------------
//-----------------------------------------------------------------
//
//                    SPALDINGS LAW OF THE WALL
//
//-----------------------------------------------------------------
//-----------------------------------------------------------------


//-----------------------------------------------------------------
//                          constructor
//                                            (public) gammi 11/09
//-----------------------------------------------------------------
DRT::ELEMENTS::FluidSurfaceWeakDBCSpaldingsLaw::FluidSurfaceWeakDBCSpaldingsLaw(
  const double chi_in,
  const double B_in  )
  :
  chi_(chi_in),
  B_  (B_in  )
{
  return;
}


//-----------------------------------------------------------------
// dynamic computation of the penalty paramter using Spaldings law
//                                            (public) gammi 11/09
//-----------------------------------------------------------------
void DRT::ELEMENTS::FluidSurfaceWeakDBCSpaldingsLaw::ComputeTauBUsingSpaldingsLaw(
  double&       tau_B,
  const double& normu,
  const double& h    ,
  const double& Cb   ,
  const double& visc )
{

  /*
  // the penalty term could be interpretaed as a traction
  // boundary condition (normally g=0 for wall bounded flows)
  //
  //      /                       \
  //     |             / n+af   \  |
  //     | v , tau  * | u    - g | |
  //     |        B    \        /  |
  //      \                       / boundaryele
  //
  //           |                 |
  //           +-----------------+
  //                t  /   "
  //               " W/ rho
  */

  /*
  // this gives rise to the following definition of the
  // friction velocity
  //
  //                 +---------+     +----------------+
  //      u    =    / t   /      =  / tau  * || n+af||
  //       tau     v   W / rho     v     B   ||u    ||
  */

  /*
  // and hence the following dimensionless value for y
  //
  //                             +-------------+
  //           h           y *  / tau  * ||u||
  //            b      +       v     B
  //      y = ---- -> y  = ----------------------
  //           C                     nu
  //            b
  //                              +
  // note that y is constant but y  is depending on tau !
  //                                                   B
  //
  //              +-------------+       +-------------+
  //        h *  / tau  * ||u||        / tau  * ||u||
  //    +    b  v     B               v     B
  //   y  = ---------------------- = -------------------
  //              C  * nu                  tau
  //               b                          B,0
  */

  /*
  // accordingly, we are able to define the dimensioneless velocity
  //                                            +-------------+
  //                  ||  n+af ||              /  ||  n+af||
  //        +         || u     ||             /   || u    ||
  //       u  = ---------------------- =     /   -------------
  //              +------------------+    \ /        tau
  //             / tau  * ||  n+af ||      v            B
  //            v     B   || u     ||
  */

  // we assume a boundary layer thickness of
  //
  //
  //                    h
  //                     b       nu
  //               y = ---- = --------
  //                    C      tau
  //                     b        B,0
  //
  // (proportional to the grid size normal to the wall)
  const double y =h/Cb;

  // iterate until the residual of the Spalding equation is 0
  double res=SpaldingResidual(y,visc,tau_B,normu);

  int count = 0;

  while(res*res>1e-6)
  {
    const double drdtauB=JacobianSpaldingResidual(y,visc,tau_B,normu);

    if(drdtauB <1e-10)
    {
      dserror("(Nearly) singular Jacobian of Spaldings equation");
    }

    double inc = res/drdtauB;

    // do damping to avoid negative values of tau_B (robustness)
    while(tau_B-inc < 0)
    {
      inc/=2.0;
    }

    // get jacobian, do damped Newton step
    tau_B-=inc;

    // get residual of Spaldings equation (law of the wall)
    res   = SpaldingResidual(y,visc,tau_B,normu);

    ++count;
    if(count>100)
    {
      dserror("no convergence in 100 steps in Newton iteration during solution of Spaldings equation\n");
    }
  }

  return;
}



//-----------------------------------------------------------------
//             evaluate the residual of Spaldings law of the wall
//                                           (private) gammi 11/09
//-----------------------------------------------------------------
double DRT::ELEMENTS::FluidSurfaceWeakDBCSpaldingsLaw::SpaldingResidual(
  const double y     ,
  const double visc  ,
  const double tau_B ,
  const double normu
  )
{
  // get dimensionless velocity
  const double up = Uplus(normu,tau_B);

  //      +
  // get y , a dimensionless boundary layer thickness
  const double yp = Yplus(normu,tau_B,visc,y);

  /*
  // Evaluate Spaldings law of the wall
  //                                /                                                  \
  //                               |                           /     +\ 2    /     +\ 3 |
  //                               |       +                  | chi*u  |    | chi*u  |  |
  //  +     / +\     +    -chi*B   |  chi*u               +    \      /      \      /   |
  // y = f | u  | = u  + e       * | e       - 1.0 - chi*u  - ----------- - ----------- |
  //        \  /                   |                              2.0           6.0     |
  //                               |                                                    |
  //                                \                                                  /
  */
  return(yp-(up+exp(-chi_*B_)*(exp(chi_*up)-1.0-chi_*up*(1+chi_*up/2.0*(1.0+chi_*up/3.0)))));
}

//-----------------------------------------------------------------
//     evaluate the residual of Spaldings law of the wall
//                                           (private) gammi 11/09
//-----------------------------------------------------------------
double DRT::ELEMENTS::FluidSurfaceWeakDBCSpaldingsLaw::JacobianSpaldingResidual(
  const double y    ,
  const double visc ,
  const double tau_B,
  const double normu)
{
  // get dimensionless velocity
  const double up = Uplus(normu,tau_B);

  // compute the derivative of the Spalding residual w.r.t. tau_B
  double drdtauB = y/(2.0*visc*sqrt(tau_B))*sqrt(normu);

  drdtauB += (1+chi_*exp(-chi_*B_)*(exp(chi_*up)-1.0-chi_*up*(1.0+0.5*chi_*up)))
             *0.5*sqrt(normu)/(sqrt(tau_B)*sqrt(tau_B)*sqrt(tau_B));

  return(drdtauB);
}


//-----------------------------------------------------------------
// compute dimensionless velocity u+
//                                           (private) gammi 11/09
//-----------------------------------------------------------------
double DRT::ELEMENTS::FluidSurfaceWeakDBCSpaldingsLaw::Uplus(
  const double normu,
  const double tau_B)
{
  /*
  // define dimensionless velocity
  //                                            +-------------+
  //                  ||  n+af ||              /  ||  n+af||
  //        +         || u     ||             /   || u    ||
  //       u  = ---------------------- =     /   -------------
  //              +------------------+    \ /        tau
  //             / tau  * ||  n+af ||      v            B
  //            v     B   || u     ||
  */
  return(sqrt(normu/tau_B));
}

//-----------------------------------------------------------------
// compute dimensionless thickness of modeled layer y+
//                                           (private) gammi 11/09
//-----------------------------------------------------------------
double DRT::ELEMENTS::FluidSurfaceWeakDBCSpaldingsLaw::Yplus(
  const double normu,
  const double tau_B,
  const double visc ,
  const double y    )
{
  /*
  //  +
  // y  is a dimensionless boundary layer thickness (law of the
  // wall is a model for the flow between this point and 0)
  //
  //            +-------------+              +-------------+
  //           / tau  * ||u||          h *  / tau  * ||u||
  //    +     v     B                   b  v     B
  //   y  =  ------------------- * y = ---------------------- =
  //                 nu                      C  * nu
  //                                          b
  //           +-------------+
  //          / tau  * ||u||
  //         v     B
  //      = -------------------
  //              tau
  //                 B,0
  //
  //
  // note that this means that in the first iteration we initialise
  // y+ and u+ as the near wall limit of Spaldings equation
  //
  //                     +-------------+
  //                    /  ||  n+af||
  //     +    +        /   || u    ||
  //    y  = u  =     /   -------------
  //               \ /        tau
  //                v            B
  */

  return((sqrt(tau_B*normu)/visc)*y);
}


//-----------------------------------------------------------------
//                       FluidLineWeakDBC
//-----------------------------------------------------------------


//-----------------------------------------------------------------
//                       empty constructor
//-----------------------------------------------------------------

template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype>
DRT::ELEMENTS::FluidLineWeakDBC<distype,pdistype>::FluidLineWeakDBC()
{
  // pointer to class FluidEleParameter (access to the general parameter)
  fldpara_ = DRT::ELEMENTS::FluidEleParameter::Instance();

  return;
}

//-----------------------------------------------------------------
//             evaluate implementation for weak dbcs
//-----------------------------------------------------------------

template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype>
int DRT::ELEMENTS::FluidLineWeakDBC<distype,pdistype>::EvaluateWeakDBC(
  FluidBoundary*                lineele       ,
  ParameterList&             params        ,
  DRT::Discretization&       discretization,
  vector<int>&               lm            ,
  Epetra_SerialDenseMatrix&  elemat_epetra ,
  Epetra_SerialDenseVector&  elevec_epetra )
{
  //--------------------------------------------------
  // get the condition information
  RCP<DRT::Condition> wdbc_cond
    =
    params.get<RCP<DRT::Condition> >("condition");

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

  // sanity check
  const string* deftauB
    =
    (*wdbc_cond).Get<string>("Definition of penalty parameter");

  if(*deftauB!="constant")
  {
    dserror("you cannot apply Spaldings law etc in 2d!\n");
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
  // apply time curve at (n+1) for all time integration schemes (dirichlet condition)
  // time_ in fluid3_parameter is time(n+alphaF) in the case of genalpha
  // therefore, it need to be reset to time(n+1)
  const double time = fldpara_->Time()+(1-fldpara_->AlphaF())*fldpara_->Dt();
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = (*wdbc_cond).Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
  curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and switches from the condition
  // (assumed to be constant on element boundary)
  const vector<int>* functions = (*wdbc_cond).Get<vector<int> >   ("funct");

  // I hope we have a linear element.
  // Ciarlet PG. The Finite element method for elliptic
  // problems. Amsterdam: North-Holland; 1978.

  // Bazilevs Michler etal use 4.0 for quadratic nurbs as well
  // (in combination with a dynamic computation of tau, so lets
  // try it as well)
  if(lineele->ParentElement()->Shape()!=DRT::Element::quad4 && lineele->ParentElement()->Shape()!=DRT::Element::nurbs9)
  {
    dserror("Cb up to now only implemented for trilinear or nurbs elements");
  }

  // find out whether to apply weak dbcs only in normal direction
  bool onlynormal = false;

  const string* active_components
    =
    (*wdbc_cond).Get<string>("Directions to apply weak dbc");

  if(*active_components=="all_directions")
  {
     onlynormal=false;
  }
  else if(*active_components=="only_in_normal_direction")
  {
     onlynormal=true;
  }
  else
  {
    dserror("Unknown definition of active components: %s",(*active_components).c_str());
  }

  // optional scaling of penalty parameter
  const double scaling
    =
    (*wdbc_cond).GetDouble("TauBscaling");

  const double Cb = 4.0*scaling;

  // get value for boundary condition
  const vector<double>* val = (*wdbc_cond).Get<vector<double> >("val");

  // get time integration parameter
  //const double afgdt = params.get<double>("afgdt");
  //const double gdt   = params.get<double>("gdt");

  const double timefac = fldpara_->TimeFac();
  const double timefacpre   = fldpara_->TimeFacPre();
  const double timefacrhs   = fldpara_->TimeFacRhs();

  //--------------------------------------------------
  // get parent elements location vector and ownerships

  // the vectors have been allocated outside in
  // EvaluateConditionUsingParentData
  RCP<vector<int> > plm
    =
    params.get<RCP<vector<int> > >("plm");
  RCP<vector<int> > plmowner
    =
    params.get<RCP<vector<int> > >("plmowner");
  RCP<vector<int> > plmstride
    =
    params.get<RCP<vector<int> > >("plmstride");

  lineele->ParentElement()->LocationVector(discretization,*plm,*plmowner,*plmstride);

  // --------------------------------------------------
  // Reshape element matrices and vectors and init to zero, construct views
  const int eledim = 3*piel;
  elemat_epetra.Shape(eledim,eledim);
  elevec_epetra.Size (eledim);
  //
  LINALG::Matrix<eledim,eledim> elemat(elemat_epetra.A(),true);
  LINALG::Matrix<eledim,     1> elevec(elevec_epetra.A(),true);


  // --------------------------------------------------
  // extract velocities from global distributed vectors

  // velocities (intermediate time step, n+alpha_F)
  RCP<const Epetra_Vector> velaf
    =
    discretization.GetState("velaf");
  if (velaf==null)
    dserror("Cannot get state vector 'velaf'");

  vector<double> mypvelaf((*plm).size());
  DRT::UTILS::ExtractMyValues(*velaf,mypvelaf,*plm);

  // velocities n+1
  vector<double> mypvelnp((*plm).size());

  if((fldpara_->TimeAlgo()==INPAR::FLUID::timeint_gen_alpha) or
      (fldpara_->TimeAlgo()==INPAR::FLUID::timeint_npgenalpha))
  {
    // velocities (intermediate time step, n+1)
    RCP<const Epetra_Vector> velnp
      =
      discretization.GetState("velnp");
    if (velnp==null)
      dserror("Cannot get state vector 'velnp'");

    DRT::UTILS::ExtractMyValues(*velnp,mypvelnp,*plm);
  }
  // mypvelnp = mypvelaf
  else
    DRT::UTILS::ExtractMyValues(*velaf,mypvelnp,*plm);

  vector<double> myedispnp ((lm  ).size());
  vector<double> mypedispnp((*plm).size());
  if (lineele->ParentElement()->IsAle())
  {
    // mesh displacements, new time step, n+1
    RCP<const Epetra_Vector> dispnp
      =
      discretization.GetState("dispnp");
    if (dispnp==null)
    {
      dserror("Cannot get state vector 'dispnp'");
    }

    DRT::UTILS::ExtractMyValues(*dispnp,myedispnp ,lm  );
    DRT::UTILS::ExtractMyValues(*dispnp,mypedispnp,*plm);
  }

  //--------------------------------------------------
  //                GET PARENT DATA
  //--------------------------------------------------

  // extract intermediate velocities
  for(int i=0;i<piel;++i)
  {
    const int ti=3*i;

    pevelaf_(0,i) = mypvelaf[  ti];
    pevelaf_(1,i) = mypvelaf[1+ti];
  }

  // extract current velocities and pressure
  for(int i=0;i<piel;++i)
  {
    const int ti=3*i;

    pevelnp_(0,i) = mypvelnp[  ti];
    pevelnp_(1,i) = mypvelnp[1+ti];

    peprenp_(  i) = mypvelnp[2+ti];
  }

  if (lineele->ParentElement()->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      const int ti=3*i;

      pedispnp_(0,i) = mypedispnp[  ti];
      pedispnp_(1,i) = mypedispnp[1+ti];
    }

    for (int i=0;i<piel;++i)
    {
      const int ti=3*i;

      edispnp_(0,i) = myedispnp[  ti];
      edispnp_(1,i) = myedispnp[1+ti];
    }
  }

  // extract node coords
  for(int i=0;i<piel;++i)
  {
    pxyze_(0,i)=lineele->ParentElement()->Nodes()[i]->X()[0];
    pxyze_(1,i)=lineele->ParentElement()->Nodes()[i]->X()[1];
  }

  if (lineele->ParentElement()->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      pxyze_(0,i) += pedispnp_(0,i);
      pxyze_(1,i) += pedispnp_(1,i);
    }
  }

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = lineele->ParentElement()->Material();

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(mat->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    // we need the kinematic viscosity here
    visc = actmat->Viscosity()/actmat->Density();
    if (actmat->Density() != 1.0)
      dserror("density 1.0 expected: the density need to be included in the linearization terms");
  }
  else
  {
    dserror("up to now I expect a constant viscosity to inforce weak DBCs\n");
  }

  //--------------------------------------------------
  //          GET BOUNDARY ELEMENT DATA
  //--------------------------------------------------

  // local surface id
  int lineid =lineele->SurfaceNumber();

  // extract node coords
  for(int i=0;i<iel;++i)
  {
    xyze_(0,i)=lineele->Nodes()[i]->X()[0];
    xyze_(1,i)=lineele->Nodes()[i]->X()[1];
  }

  if (lineele->ParentElement()->IsAle())
  {
    for (int i=0;i<iel;++i)
    {
      xyze_(0,i) += edispnp_(0,i);
      xyze_(1,i) += edispnp_(1,i);
    }
  }

  //--------------------------------------------------
  // get gausspoints to integrate over boundary element

  // get gauss rule
  DRT::UTILS::GaussRule1D gaussrule=DRT::UTILS::intrule1D_undefined;
  switch (distype)
  {
  case DRT::Element::line2:
  {
    gaussrule = DRT::UTILS::intrule_line_2point;
    break;
  }
  case DRT::Element::nurbs3:
  {
    gaussrule = DRT::UTILS::intrule_line_3point;
    break;
  }
  default:
    dserror("invalid discretization type for fluid2line weak DBC evaluation");
  }

  // gaussian points on surface
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule);

  //--------------------------------------------------
  // the gausspoints above have to be mapped to the
  // parent element to be able to evaluate one sided
  // derivatives on the boundary
  //
  // in addition, get information on the orientation of the
  // outward normal

  Epetra_SerialDenseMatrix pqxg(intpoints.nquad,2);

  DRT::UTILS::LineGPToParentGP(pqxg     ,
                               intpoints,
                               pdistype ,
                               distype  ,
                               lineid);


  // --------------------------------------------------
  // Now do the nurbs specific stuff

  std::vector<Epetra_SerialDenseVector> mypknots(2);
  std::vector<Epetra_SerialDenseVector> myknots (1);

  Epetra_SerialDenseVector weights(iel);
  LINALG::Matrix<piel,1>   pweights;

  // orientation of outward normal
  double normalfac=1.0;

  // for isogeometric elements --- get knotvectors for parent
  // element and surface element, get weights
  if(lineele->Shape()==DRT::Element::nurbs3)
  {
    // --------------------------------------------------
    // get knotvector

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    RCP<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

    bool zero_sized_parent
      =knots->GetBoundaryEleAndParentKnots(mypknots              ,
                                           myknots               ,
                                           normalfac             ,
                                           lineele->ParentElement()->Id(),
                                           lineid                );

    if(zero_sized_parent)
    {
      return 0;
    }

    // --------------------------------------------------
    // get node weights for nurbs elements
    for (int inode=0; inode<iel; inode++)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (lineele->Nodes()[inode]);

      weights(inode) = cp->W();
    }

    // extract node coords
    for(int i=0;i<piel;++i)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (lineele->ParentElement()->Nodes()[i]);

      pweights(i) = cp->W();
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
    const double xi    = intpoints.qxg [iquad][0];

    // gaussian point in parent elements local coordinates
    const double r     = pqxg(iquad,0);
    const double s     = pqxg(iquad,1);

    if(!(distype == DRT::Element::nurbs3))
    {
      // ------------------------------------------------
      // shape function derivs of boundary element at gausspoint
      DRT::UTILS::shape_function_1D       (funct_,xi,distype);

      LINALG::Matrix<1,2> deriv1;

      DRT::UTILS::shape_function_1D_deriv1(deriv1,xi,distype);

      deriv_(0)=deriv1(0,0);
      deriv_(1)=deriv1(0,1);

    }
    else
    {
      DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv
        (funct_    ,
         deriv_    ,
         xi        ,
         myknots[0],
         weights   ,
         distype   );
    }

    // ------------------------------------------------
    // compute measure tensor for surface element and the infinitesimal
    // area element drs for the integration

    /*
      |
      |                                        0 1
      |                                       +-+-+
      |                                       | | | 0
      |         0 1             0...iel-1     +-+-+
      |        +-+-+            +-+-+-+-+     | | | .
      |        | | |      =     | | | | |  *  +-+-+ .
      |        +-+-+            +-+-+-+-+     | | | .
      |                                       +-+-+
      |                                       | | | iel-1
      |                                       +-+-+
      |
      |       dxydr               deriv       xyze^T
      |
      |
      |
      |
      |
      |                                 +-        -+
      |                                 | dx   dy  |
      |     yields             dxydr =  | --   --  |
      |                                 | dr   dr  |
      |                                 +-        -+
      |
      |
      |
    */
    dxydr_.MultiplyNN(xyze_,deriv_);
    /*
      |
      |         +-       -+   +-       -+ T
      |         | dx   dy |   | dx   dy |
      |    g =  | --   -- | * | --   -- |
      |         | dr   dr |   | dr   dr |
      |         +-       -+   +-       -+
      |
    */
    const double g = dxydr_(0)*dxydr_(0)+dxydr_(1)*dxydr_(1);

    /*
                         +-+
              sqrtg =  \/ g

    */

    dr_= sqrt(g);


    // total integration factor
    const double fac = dr_*wquad;

    // ------------------------------------------------
    // compute normal
    {
      /*
      |
      |                      +-  -+     +-  -+
      |                      | dx |     |    |
      |                      | -- |     |  0 |
      |                      | dr |     |    |
      |                      |    |     |    |
      |             1.0      | dy |     |    |
      |    n  =  --------- * | -- |  X  |  0 |
      |                      | dr |     |    |
      |          ||.....||   |    |     |    |
      |                      |    |     |    |
      |                      |  0 |     |  1 |
      |                      |    |     |    |
      |                      +-  -+     +-  -+
      |
    */
      n_(0) =  dxydr_(1);
      n_(1) = -dxydr_(0);

      const double length = n_.Norm2()*normalfac;

      for(int i=0;i<2;++i)
      {
        n_(i)/=length;
      }
    }

    // ------------------------------------------------
    // factor given by spatial function
    LINALG::Matrix<2,1> functionfac;
    for(int i=0;i<2;++i)
    {
      functionfac(i)= 1.0;
    }

    // determine coordinates of current Gauss point
    LINALG::Matrix<3,1> coordgp;
    for(int i=0;i<3;++i)
    {
      coordgp(i)= 0.0;
    }
    for (int i=0;i<iel;++i)
    {
      for(int j=0;j<2;++j)
      {
        coordgp(j)+=xyze_(j,i)*funct_(i);
      }
    }

    int functnum = -1;

    for (int node=0;node<iel;++node)
    {
      for(int dim=0;dim<2;++dim)
      {
  // factor given by spatial function
  if (functions)
  {
    functnum = (*functions)[dim];
    if (functnum>0)
    {
      // evaluate function at current gauss point (important: requires 3D position vector)
      functionfac(dim) = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dim,coordgp.A(),time,NULL);
    }
    else
    {
      functionfac(dim) = 1.0;
    }
  }
      }
    }

    // ------------------------------------------------
    // shape functions and derivs of corresponding parent at gausspoint
    if(!(pdistype == DRT::Element::nurbs9))
    {
      DRT::UTILS::shape_function_2D       (pfunct_,r,s,pdistype);
      DRT::UTILS::shape_function_2D_deriv1(pderiv_,r,s,pdistype);
    }
    else
    {
      // set gauss point coordinates
      LINALG::Matrix<2,1> gp;

      gp(0)=r;
      gp(1)=s;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
        (pfunct_ ,
         pderiv_ ,
         gp      ,
         mypknots,
         pweights,
         pdistype);
    }

    //-----------------------------------------------------
    //
    //                       +-------------+
    //                      / /  T       \ |
    //           h = 2 * \ / |  n * G * n |
    //            b       +   \          /
    //

    // get Jacobian matrix and determinant
    pxjm_=0;

    for(int i=0;i<piel;++i)
    {
      for(int rr=0;rr<2;++rr)
      {
  for(int mm=0;mm<2;++mm)
  {
    pxjm_(rr,mm)+=pderiv_(rr,i)*pxyze_(mm,i);
  }
      }
    }

    const double pdet = pxjm_(0,0)*pxjm_(1,1)-pxjm_(0,1)*pxjm_(1,0);

    // check for degenerated elements
    if (pdet < 0.0)
    {
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",
              lineele->ParentElement()->Id(),
              pdet);
    }

    //-----------------------------------------------------
    //
    //             compute global first derivates
    //
    /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

    */

    // inverse of jacobian (transposed)
    /*
          +-                 -+     +-                 -+ -1
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dx    dx    dx   |     |  dr    dr    dr   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |  =  |  --    --    --   |
          |  dy    dy    dy   |     |  ds    ds    ds   |
          |                   |     |                   |
          |  dr    ds    dt   |     |  dx    dy    dz   |
          |  --    --    --   |     |  --    --    --   |
          |  dz    dz    dz   |     |  dt    dt    dt   |
          +-                 -+     +-                 -+

    */
    pxji_(0,0) = (  pxjm_(1,1) )/pdet;
    pxji_(1,0) = (- pxjm_(1,0) )/pdet;
    pxji_(0,1) = (- pxjm_(0,1) )/pdet;
    pxji_(1,1) = (  pxjm_(0,0) )/pdet;


    //-----------------------------------------------------
    /*            +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
    */
    LINALG::Matrix<2,2> G;

    for (int nn=0;nn<2;++nn)
    {
      for (int rr=0;rr<2;++rr)
      {
        G(nn,rr) = pxji_(nn,0)*pxji_(rr,0)+pxji_(nn,1)*pxji_(rr,1);
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

    double nGn=0;
    for (int nn=0;nn<2;++nn)
    {
      for (int rr=0;rr<2;++rr)
      {
        nGn+=n_(rr)*G(rr,nn)*n_(nn);
      }
    }
    const double h =2.0/sqrt(nGn);

    //-----------------------------------------------------
    // compute global derivates at integration point
    //
    //   dN    +-----  dN (xi)    dxi
    //     i    \        i           k
    //   --- =   +     ------- * -----
    //   dx     /        dxi      dx
    //     j   +-----       k       j
    //         node k
    //
    // j : direction of derivative x/y/z
    //
    for(int nn=0;nn<piel;++nn)
    {
      for(int rr=0;rr<2;++rr)
      {
        pderxy_(rr,nn)=pxji_(rr,0)*pderiv_(0,nn)+pxji_(rr,1)*pderiv_(1,nn);
      }
    }

    //-----------------------------------------------------
    // get velocities (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    vel   (x) =   +      N (x) * vel
    //                 /        j         j
    //                +-----
    //                node j
    //
    for(int rr=0;rr<2;++rr)
    {
      velintnp_(rr)=pfunct_(0)*pevelnp_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        velintnp_(rr)+=pfunct_(nn)*pevelnp_(rr,nn);
      }
    }

    //-----------------------------------------------------
    // get pressure (n+1,i) at integration point
    //
    //                +-----
    //       n+1       \                  n+1
    //    pre   (x) =   +      N (x) * pre
    //                 /        i         i
    //                +-----
    //                node i
    //
    prenp_=pfunct_(0)*peprenp_(0);
    for(int nn=1;nn<piel;++nn)
    {
      prenp_+=pfunct_(nn)*peprenp_(nn);
    }

    //-----------------------------------------------------
    // get velocities (n+alpha_F,i) at integration point
    //
    //                 +-----
    //       n+af       \                  n+af
    //    vel    (x) =   +      N (x) * vel
    //                  /        j         j
    //                 +-----
    //                 node j
    //
    for(int rr=0;rr<2;++rr)
    {
      velintaf_(rr)=pfunct_(0)*pevelaf_(rr,0);
      for(int nn=1;nn<piel;++nn)
      {
        velintaf_(rr)+=pfunct_(nn)*pevelaf_(rr,nn);
      }
    }

    //-----------------------------------------------------
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
    for(int rr=0;rr<2;++rr)
    {
      for(int mm=0;mm<2;++mm)
      {
        vderxyaf_(rr,mm)=pderxy_(mm,0)*pevelaf_(rr,0);
        for(int nn=1;nn<piel;++nn)
        {
          vderxyaf_(rr,mm)+=pderxy_(mm,nn)*pevelaf_(rr,nn);
        }
      }
    }

    /*
    // ---------------------------------------------------
    // define (initial) penalty parameter
    //
    //
    //                         /    \
    //                        |  nu  |
    //         tau      = C * | ---- |
    //            B(,0)    b  |  h   |
    //                         \  b /
    */
    double tau_B = Cb*visc/h;

    // ---------------------------------------------------
    // velocities on boundary
    //
    //         n+af
    //        u    - u
    //                b
    //

    LINALG::Matrix<2,1> bvres;

    bvres(0)=(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);
    bvres(1)=(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);


    // ---------------------------------------------------
    // momentum flux over boundary

    const double flux=velintaf_(0)*n_(0)+velintaf_(1)*n_(1);

    //--------------------------------------------------
    // partially integrated pressure term, rescaled by gamma*dt
    /*
    // factor: 1.0
    //
    //             /            \
    //            |              |
    //          + |  v , Dp * n  |
    //            |              |
    //             \            / boundaryele
    //
    */
    for (int ui=0; ui<piel; ++ui)
    {
      for (int vi=0; vi<piel; ++vi)
      {
        elemat(vi*3  ,ui*3+2) += fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*3+1,ui*3+2) += fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(1);
      }
    }

    for (int vi=0; vi<piel; ++vi)
    {
      elevec(vi*3    ) -= fac*timefacrhs*pfunct_(vi)*n_(0)*prenp_;
      elevec(vi*3 + 1) -= fac*timefacrhs*pfunct_(vi)*n_(1)*prenp_;
    }

    //--------------------------------------------------
    // adjoint consistency term, pressure/continuity part
    /*
    // factor: gdt
    //
    //             /              \
    //            |                |
    //          - |  q , Dacc * n  |
    //            |                |
    //             \              / boundaryele
    //
    */
    for (int ui=0; ui<piel; ++ui)
    {
      for (int vi=0; vi<piel; ++vi)
      {
        elemat(vi*3+2,ui*3  ) -= fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*3+2,ui*3+1) -= fac*timefacpre*pfunct_(vi)*pfunct_(ui)*n_(1);
      }
    }

    /*
    // factor: 1.0
    //
    //             /                       \
    //            |       / n+1     \       |
    //          + |  q , | u   - u   | * n  |
    //            |       \ (i)   B /       |
    //             \                       / boundaryele
    //
    */
    for (int vi=0; vi<piel; ++vi)
    {
      elevec(vi*3+2) += fac*timefacrhs*pfunct_(vi)*
        ((velintnp_(0)-(*val)[0]*functionfac(0)*curvefac)*n_(0)
         +
         (velintnp_(1)-(*val)[1]*functionfac(1)*curvefac)*n_(1));
    }

    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    //           weak boundary conditions in all directions
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    if(!onlynormal)
    {
      const double timefacnu=fac*2.0*visc*timefac;

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
        double nabla_u_o_n_lin[2][2];

        nabla_u_o_n_lin[0][0]=timefacnu*(    pderxy_(0,ui)*n_(0)
                                         +
                                         0.5*pderxy_(1,ui)*n_(1));
        nabla_u_o_n_lin[0][1]=timefacnu*(0.5*pderxy_(0,ui)*n_(1));

        nabla_u_o_n_lin[1][0]=timefacnu*(0.5*pderxy_(1,ui)*n_(0));
        nabla_u_o_n_lin[1][1]=timefacnu*(0.5*pderxy_(0,ui)*n_(0)
                                         +
                                         pderxy_(1,ui)*n_(1));


        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*3    ,ui*3    ) -= pfunct_(vi)*nabla_u_o_n_lin[0][0];
          elemat(vi*3    ,ui*3 + 1) -= pfunct_(vi)*nabla_u_o_n_lin[0][1];

          elemat(vi*3 + 1,ui*3    ) -= pfunct_(vi)*nabla_u_o_n_lin[1][0];
          elemat(vi*3 + 1,ui*3 + 1) -= pfunct_(vi)*nabla_u_o_n_lin[1][1];
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

      {
        double nabla_u_o_n[2];
        nabla_u_o_n[0]=fac*timefacrhs*2.0*visc*(vderxyaf_(0,0)*n_(0)+0.5*(vderxyaf_(0,1)+vderxyaf_(1,0))*n_(1));
        nabla_u_o_n[1]=fac*timefacrhs*2.0*visc*(0.5*(vderxyaf_(1,0)+vderxyaf_(0,1))*n_(0)+vderxyaf_(1,1) *n_(1));

        for (int vi=0; vi<piel; ++vi)
        {
          elevec(vi*3    ) += pfunct_(vi)*nabla_u_o_n[0];
          elevec(vi*3 + 1) += pfunct_(vi)*nabla_u_o_n[1];
        }
      }

      //--------------------------------------------------
      // (adjoint) consistency term, viscous part

      const double consistencytimefac=fac*2.0*visc*wd_gamma*timefac;


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

      LINALG::Matrix<2,2> nabla_s_w_o_n;
      for (int vi=0; vi<piel; ++vi)
      {
        nabla_s_w_o_n(0,0)=consistencytimefac*(n_(0)*(    pderxy_(0,vi))+n_(1)*(0.5*pderxy_(1,vi)));
        nabla_s_w_o_n(0,1)=consistencytimefac*(n_(0)*(0.5*pderxy_(1,vi)));

        nabla_s_w_o_n(1,0)=consistencytimefac*(n_(1)*(0.5*pderxy_(0,vi)));
        nabla_s_w_o_n(1,1)=consistencytimefac*(n_(0)*(0.5*pderxy_(0,vi))+n_(1)*(    pderxy_(1,vi)));


        for (int ui=0; ui<piel; ++ui)
        {
          elemat(vi*3    ,ui*3    ) -= pfunct_(ui)*nabla_s_w_o_n(0,0);
          elemat(vi*3    ,ui*3 + 1) -= pfunct_(ui)*nabla_s_w_o_n(0,1);

          elemat(vi*3 + 1,ui*3    ) -= pfunct_(ui)*nabla_s_w_o_n(1,0);
          elemat(vi*3 + 1,ui*3 + 1) -= pfunct_(ui)*nabla_s_w_o_n(1,1);
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
        elevec(vi*3    ) += fac*timefacrhs*2.0*visc*wd_gamma*(
          n_(0)*bvres(0)*(    pderxy_(0,vi))
          +
          n_(1)*bvres(0)*(0.5*pderxy_(1,vi))
          +
          n_(0)*bvres(1)*(0.5*pderxy_(1,vi)));

        elevec(vi*3 + 1) += fac*timefacrhs*2.0*visc*wd_gamma*(
          n_(1)*bvres(0)*(0.5*pderxy_(0,vi))
          +
          n_(0)*bvres(1)*(0.5*pderxy_(0,vi))
          +
          n_(1)*bvres(1)*(    pderxy_(1,vi)));
      }

      //--------------------------------------------------
      // adjoint consistency term, convective part (only on inflow)

      if(flux<0)
      {
        if(complete_linearisation)
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
              elemat(vi*3    ,ui*3    ) -= fac*timefac*pfunct_(vi)*n_(0)*bvres(0);
              elemat(vi*3    ,ui*3 + 1) -= fac*timefac*pfunct_(vi)*n_(1)*bvres(0);

              elemat(vi*3 + 1,ui*3    ) -= fac*timefac*pfunct_(vi)*n_(0)*bvres(1);
              elemat(vi*3 + 1,ui*3 + 1) -= fac*timefac*pfunct_(vi)*n_(1)*bvres(1);
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

          const double fluxtimefac =fac*timefac*flux;

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*3    ,ui*3    ) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
              elemat(vi*3 + 1,ui*3 + 1) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
            }
          }
        } // end if full_linearisation

        const double fluxfac =fac*timefacrhs*flux;

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
          elevec(vi*3    ) += fluxfac*pfunct_(vi)*bvres(0);
          elevec(vi*3 + 1) += fluxfac*pfunct_(vi)*bvres(1);
        }
      } // end if flux<0, i.e. boundary is an inflow boundary

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

      const double penaltytimefac=timefac*tau_B*fac;

      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          const double temp=penaltytimefac*pfunct_(ui)*pfunct_(vi);

          elemat(vi*3    ,ui*3    ) += temp;
          elemat(vi*3 + 1,ui*3 + 1) += temp;
        }
      }

      const double penaltyfac=tau_B*fac*timefacrhs;

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
        elevec(vi*3    ) -= penaltyfac*pfunct_(vi)*bvres(0);
        elevec(vi*3 + 1) -= penaltyfac*pfunct_(vi)*bvres(1);
      }
    } // !onlynormal
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    //          weak boundary conditions only in normal direction
    //---------------------------------------------------------------------
    //---------------------------------------------------------------------
    else
    {
      //--------------------------------------------------
      // partially integrated viscous term

      const double timefacnu=fac*2.0*visc*timefac;

      /*
      // factor: 2*nu
      //
      //    /                             \
      //   |                   s           |
      // + |  v * n , n * nabla  Dacc * n  |
      //   |                               |
      //    \                             / boundaryele
      //
      */
      for (int ui=0; ui<piel; ++ui)
      {
        const double aux=timefacnu*(pderxy_(0,ui)*n_(0)+pderxy_(1,ui)*n_(1));

        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*3    ,ui*3    ) -= pfunct_(vi)*n_(0)*n_(0)*aux;
          elemat(vi*3    ,ui*3 + 1) -= pfunct_(vi)*n_(0)*n_(1)*aux;

          elemat(vi*3 + 1,ui*3    ) -= pfunct_(vi)*n_(1)*n_(0)*aux;
          elemat(vi*3 + 1,ui*3 + 1) -= pfunct_(vi)*n_(1)*n_(1)*aux;
        }
      }

      /*
      // factor: 2*nu
      //
      //    /                             \
      //   |                   s  n+af     |
      // + |  v * n , n * nabla  u    * n  |
      //   |                               |
      //    \                             / boundaryele
      //
      */
      double n_o_nabla_u_o_n =
        vderxyaf_(0,0)*n_(0)*n_(0)
        +
        vderxyaf_(1,1)*n_(1)*n_(1)
        +
        (vderxyaf_(0,1)+vderxyaf_(1,0))*n_(0)*n_(1);

      for (int vi=0; vi<piel; ++vi)
      {
        elevec(vi*3    ) += fac*timefacrhs*2.0*visc*wd_gamma*pfunct_(vi)*n_(0)*n_o_nabla_u_o_n;
        elevec(vi*3 + 1) += fac*timefacrhs*2.0*visc*wd_gamma*pfunct_(vi)*n_(1)*n_o_nabla_u_o_n;
      }

      //--------------------------------------------------
      // (adjoint) consistency term, viscous part

      const double consistencytimefac=fac*2.0*visc*wd_gamma*timefac;

      /*
      // factor: 2*nu*gamma_wd
      //
      //    /                               \
      //   |            s                    |
      // + |   n * nabla  w * n ,  Dacc * n  |
      //   |                                 |
      //    \                               / boundaryele
      //
      */

      for (int vi=0; vi<piel; ++vi)
      {
        for (int ui=0; ui<piel; ++ui)
        {
          const double aux=pderxy_(0,vi)*n_(0)+pderxy_(1,vi)*n_(1);

          elemat(vi*3    ,ui*3    ) -= aux*n_(0)*n_(0)*pfunct_(ui)*consistencytimefac;
          elemat(vi*3    ,ui*3 + 1) -= aux*n_(0)*n_(1)*pfunct_(ui)*consistencytimefac;

          elemat(vi*3 + 1,ui*3    ) -= aux*n_(1)*n_(0)*pfunct_(ui)*consistencytimefac;
          elemat(vi*3 + 1,ui*3 + 1) -= aux*n_(1)*n_(1)*pfunct_(ui)*consistencytimefac;
        }
      }
      /*
      // factor: 2*nu*gamma_wd
      //
      //    /                                        \
      //   |            s          / n+af     \       |
      // + |   n * nabla  w * n , | u    - u   | * n  |
      //   |                       \        b /       |
      //    \                                        / boundaryele
      //
      */

      double bvres_o_n=bvres(0)*n_(0)+bvres(1)*n_(1);

      for (int vi=0; vi<piel; ++vi)
      {
        double aux=(pderxy_(0,vi)*n_(0)+pderxy_(1,vi)*n_(1));

        elevec(vi*3    ) += fac*timefacrhs*2.0*visc*wd_gamma*aux*n_(0)*bvres_o_n;
        elevec(vi*3 + 1) += fac*timefacrhs*2.0*visc*wd_gamma*aux*n_(1)*bvres_o_n;
      }

      //--------------------------------------------------
      // penalty term

      const double penaltytimefac=timefac*tau_B*fac;

      /*
      // factor: nu*Cb/h*afgdt
      //
      //    /                  \
      //   |                    |
      // + |  w o n , Dacc o n  |
      //   |                    |
      //    \                  / boundaryele
      //
      */
      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          elemat(vi*3    ,ui*3    ) += penaltytimefac*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
          elemat(vi*3    ,ui*3 + 1) += penaltytimefac*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);

          elemat(vi*3 + 1,ui*3    ) += penaltytimefac*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
          elemat(vi*3 + 1,ui*3 + 1) += penaltytimefac*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
        }
      }

      const double penaltyfac=tau_B*fac*timefacrhs;

      /*
      // factor: nu*Cb/h
      //
      //    /                           \
      //   |            / n+af   \       |
      // + |   w o n , | u    - u | o n  |
      //   |            \  b     /       |
      //    \                           / boundaryele
      //
      */

      for (int vi=0; vi<piel; ++vi)
      {
        elevec(vi*3    ) -= penaltyfac*pfunct_(vi)*n_(0)*bvres_o_n;
        elevec(vi*3 + 1) -= penaltyfac*pfunct_(vi)*n_(1)*bvres_o_n;
      }

      //--------------------------------------------------
      // adjoint consistency term, convective part (only on inflow)

      if(flux<0)
      {
        if(complete_linearisation)
        {
          /*
          // These linearisations have only to be included if
          // u*n is negative --- otherwise they're nonesense
          //
          // factor: afgdt
          //
          //    /                                              \
          //   |    /        \    /     \     / n+af     \      |
          // - |   | Dacc * n |  | w o n | , | u    - u   | o n |
          //   |    \        /    \     /     \        b /      |
          //    \                                              / boundaryele, inflow
          //
          */

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*3    ,ui*3    ) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
              elemat(vi*3    ,ui*3 + 1) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);

              elemat(vi*3 + 1,ui*3    ) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
              elemat(vi*3 + 1,ui*3 + 1) -= fac*timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
            }
          }

          const double fluxtimefac =fac*timefac*flux;

          /*
          // factor: afgdt
          //
          //    /                                   \
          //   |    / n+af   \   /     \             |
          // - |   | u    * n | | w o n | , Dacc o n |
          //   |    \        /   \     /             |
          //    \  |          |                     / boundaryele, inflow
          //       +----------+
          //           <0
          */
          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*3    ,ui*3    ) -= fluxtimefac*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
              elemat(vi*3    ,ui*3 + 1) -= fluxtimefac*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);

              elemat(vi*3 + 1,ui*3    ) -= fluxtimefac*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
              elemat(vi*3 + 1,ui*3 + 1) -= fluxtimefac*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
            }
          }

        } // end complete_linearisation

        const double fluxfac =fac*timefacrhs*flux;

        /*
        // factor: 1
        //
        //    /                                            \
        //   |    / n+af   \   /     \    / n+af     \      |
        // - |   | u    * n | | w o n |, | u    - u   | o n |
        //   |    \        /   \     /    \        b /      |
        //    \  |          |                              / boundaryele, inflow
        //       +----------+
        //           <0
        */
        for (int vi=0; vi<piel; ++vi)
        {
          elevec(vi*3    ) += fluxfac*pfunct_(vi)*n_(0)*bvres_o_n;
          elevec(vi*3 + 1) += fluxfac*pfunct_(vi)*n_(1)*bvres_o_n;
        }
      } // end if flux<0, i.e. boundary is an inflow boundary

    } // onlynormal

  } // end gaussloop

  return 0;
}

