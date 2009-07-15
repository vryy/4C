#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2_weak_dbc.H"

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
DRT::ELEMENTS::Fluid2LineWeakDBCInterface* DRT::ELEMENTS::Fluid2LineWeakDBCInterface::Impl(
  DRT::ELEMENTS::Fluid2Line* f2line
  )
{
  switch (f2line->Shape())
  {
  case DRT::Element::line2:
  {
    static Fluid2LineWeakDBC<DRT::Element::line2,DRT::Element::quad4>* fline2;

    if(f2line->parent_->Shape()==DRT::Element::quad4)
    {
      if (fline2==NULL)
      {
        fline2 = new Fluid2LineWeakDBC<DRT::Element::line2,DRT::Element::quad4>();
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
    static Fluid2LineWeakDBC<DRT::Element::nurbs3,DRT::Element::nurbs9>* flinen3;

    if(f2line->parent_->Shape()==DRT::Element::nurbs9)
    {
      if (flinen3==NULL)
      {
        flinen3 = new Fluid2LineWeakDBC<DRT::Element::nurbs3,DRT::Element::nurbs9>();
      }
    }
    else
    {
      dserror("expected combination nurbs3/nurbs9 for line/parent pair");
    }

    return flinen3;
  }
  default:
  {
    dserror("shape %d (%d nodes) not supported by weak DBC", f2line->Shape(), f2line->NumNode());
  }
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
DRT::ELEMENTS::Fluid2LineWeakDBC<distype,pdistype>::Fluid2LineWeakDBC()
{
  return;
}

//-----------------------------------------------------------------
//             evaluate implementation for weak dbcs
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype>
int DRT::ELEMENTS::Fluid2LineWeakDBC<distype,pdistype>::EvaluateWeakDBC(
  Fluid2Line*                lineele       ,
  ParameterList&             params        ,
  DRT::Discretization&       discretization,
  vector<int>&               lm            ,
  Epetra_SerialDenseMatrix&  elemat_epetra ,
  Epetra_SerialDenseVector&  elevec_epetra )
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
  // Ciarlet PG. The Finite element method for elliptic
  // problems. Amsterdam: North-Holland; 1978.

  // Bazilevs Michler etal use 4.0 for quadratic nurbs as well
  // (in combination with a dynamic computation of tau, so lets
  // try it as well)
  if(lineele->parent_->Shape()!=Fluid2::quad4 && lineele->parent_->Shape()!=Fluid2::nurbs9)
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

  lineele->parent_->LocationVector(discretization,*plm,*plmowner);

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
  RefCountPtr<const Epetra_Vector> velaf
    =
    discretization.GetState("u and p (n+alpha_F,trial)");
  if (velaf==null)
  {
    dserror("Cannot get state vector 'velaf'");
  }

  vector<double> mypvelaf((*plm).size());
  DRT::UTILS::ExtractMyValues(*velaf,mypvelaf,*plm);

  // velocities (intermediate time step, n+1)
  RefCountPtr<const Epetra_Vector> velnp
    =
    discretization.GetState("u and p (n+1      ,trial)");
  if (velnp==null)
  {
    dserror("Cannot get state vector 'velnp'");
  }

  vector<double> mypvelnp((*plm).size());
  DRT::UTILS::ExtractMyValues(*velnp,mypvelnp,*plm);

  vector<double> myedispnp ((lm  ).size());
  vector<double> mypedispnp((*plm).size());
  if (lineele->parent_->IsAle())
  {
    // mesh displacements, new time step, n+1
    RefCountPtr<const Epetra_Vector> dispnp
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

  if (lineele->parent_->IsAle())
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
    pxyze_(0,i)=lineele->parent_->Nodes()[i]->X()[0];
    pxyze_(1,i)=lineele->parent_->Nodes()[i]->X()[1];
  }

  if (lineele->parent_->IsAle())
  {
    for (int i=0;i<piel;++i)
    {
      pxyze_(0,i) += pedispnp_(0,i);
      pxyze_(1,i) += pedispnp_(1,i);
    }
  }

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = lineele->parent_->Material();

  if( mat->MaterialType()    != INPAR::MAT::m_carreauyasuda
      && mat->MaterialType() != INPAR::MAT::m_modpowerlaw
      && mat->MaterialType() != INPAR::MAT::m_fluid)
          dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(mat->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    visc = actmat->Viscosity();
  }
  else
  {
    dserror("up to now I expect a constant viscosity to inforce weak DBCs\n");
  }

  //--------------------------------------------------
  //          GET BOUNDARY ELEMENT DATA
  //--------------------------------------------------

  // local surface id
  int lineid =lineele->lline_;

  // extract node coords
  for(int i=0;i<iel;++i)
  {
    xyze_(0,i)=lineele->Nodes()[i]->X()[0];
    xyze_(1,i)=lineele->Nodes()[i]->X()[1];
  }

  if (lineele->parent_->IsAle())
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
  if(lineele->Shape()==Fluid2::nurbs3)
  {
    // --------------------------------------------------
    // get knotvector

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    RefCountPtr<DRT::NURBS::Knotvector> knots=(*nurbsdis).GetKnotVector();

    bool zero_sized_parent
      =knots->GetBoundaryEleAndParentKnots(mypknots              ,
                                           myknots               ,
                                           normalfac             ,
                                           lineele->parent_->Id(),
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
        dynamic_cast<DRT::NURBS::ControlPoint* > (lineele->parent_->Nodes()[i]);

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
	    // evaluate function at current gauss point
	    functionfac(dim) = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dim,coordgp.A(),0.0,NULL);
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
              lineele->parent_->Id(),
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
        elemat(vi*3  ,ui*3+2) += fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*3+1,ui*3+2) += fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(1);
      }
    }

    for (int vi=0; vi<piel; ++vi)
    {
      elevec(vi*3    ) -= fac*pfunct_(vi)*n_(0)*prenp_;
      elevec(vi*3 + 1) -= fac*pfunct_(vi)*n_(1)*prenp_;
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
        elemat(vi*3+2,ui*3  ) -= fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*3+2,ui*3+1) -= fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(1);
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
      elevec(vi*3+2) += fac*pfunct_(vi)*
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
      const double timefacnu=fac*2.0*visc*afgdt;

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
        nabla_u_o_n[0]=fac*2.0*visc*(vderxyaf_(0,0)*n_(0)+0.5*(vderxyaf_(0,1)+vderxyaf_(1,0))*n_(1));
        nabla_u_o_n[1]=fac*2.0*visc*(0.5*(vderxyaf_(1,0)+vderxyaf_(0,1))*n_(0)+vderxyaf_(1,1) *n_(1));

        for (int vi=0; vi<piel; ++vi)
        {
          elevec(vi*3    ) += pfunct_(vi)*nabla_u_o_n[0];
          elevec(vi*3 + 1) += pfunct_(vi)*nabla_u_o_n[1];
        }
      }

      //--------------------------------------------------
      // (adjoint) consistency term, viscous part

      const double consistencytimefac=fac*2.0*visc*wd_gamma*afgdt;


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
        elevec(vi*3    ) += fac*2.0*visc*wd_gamma*(
          n_(0)*bvres(0)*(    pderxy_(0,vi))
          +
          n_(1)*bvres(0)*(0.5*pderxy_(1,vi))
          +
          n_(0)*bvres(1)*(0.5*pderxy_(1,vi)));

        elevec(vi*3 + 1) += fac*2.0*visc*wd_gamma*(
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
          const double timefac=fac*afgdt;

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*3    ,ui*3    ) -= timefac*pfunct_(vi)*n_(0)*bvres(0);
              elemat(vi*3    ,ui*3 + 1) -= timefac*pfunct_(vi)*n_(1)*bvres(0);

              elemat(vi*3 + 1,ui*3    ) -= timefac*pfunct_(vi)*n_(0)*bvres(1);
              elemat(vi*3 + 1,ui*3 + 1) -= timefac*pfunct_(vi)*n_(1)*bvres(1);
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

          const double fluxtimefac =fac*afgdt*flux;

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*3    ,ui*3    ) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
              elemat(vi*3 + 1,ui*3 + 1) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
            }
          }
        } // end if full_linearisation

        const double fluxfac =fac*flux;

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

      const double penaltytimefac=afgdt*tau_B*fac;

      for (int ui=0; ui<piel; ++ui)
      {
        for (int vi=0; vi<piel; ++vi)
        {
          const double temp=penaltytimefac*pfunct_(ui)*pfunct_(vi);

          elemat(vi*3    ,ui*3    ) += temp;
          elemat(vi*3 + 1,ui*3 + 1) += temp;
        }
      }

      const double penaltyfac=tau_B*fac;

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

      const double timefacnu=fac*2.0*visc*afgdt;

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
        elevec(vi*3    ) += fac*2.0*visc*wd_gamma*pfunct_(vi)*n_(0)*n_o_nabla_u_o_n;
        elevec(vi*3 + 1) += fac*2.0*visc*wd_gamma*pfunct_(vi)*n_(1)*n_o_nabla_u_o_n;
      }

      //--------------------------------------------------
      // (adjoint) consistency term, viscous part

      const double consistencytimefac=fac*2.0*visc*wd_gamma*afgdt;

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

        elevec(vi*3    ) += fac*2.0*visc*wd_gamma*aux*n_(0)*bvres_o_n;
        elevec(vi*3 + 1) += fac*2.0*visc*wd_gamma*aux*n_(1)*bvres_o_n;
      }

      //--------------------------------------------------
      // penalty term

      const double penaltytimefac=afgdt*tau_B*fac;

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

      const double penaltyfac=tau_B*fac;

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
          const double timefac=fac*afgdt;

          for (int ui=0; ui<piel; ++ui)
          {
            for (int vi=0; vi<piel; ++vi)
            {
              elemat(vi*3    ,ui*3    ) -= timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(0)*pfunct_(ui);
              elemat(vi*3    ,ui*3 + 1) -= timefac*bvres_o_n*pfunct_(vi)*n_(0)*n_(1)*pfunct_(ui);

              elemat(vi*3 + 1,ui*3    ) -= timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(0)*pfunct_(ui);
              elemat(vi*3 + 1,ui*3 + 1) -= timefac*bvres_o_n*pfunct_(vi)*n_(1)*n_(1)*pfunct_(ui);
            }
          }

          const double fluxtimefac =fac*afgdt*flux;

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

        const double fluxfac =fac*flux;

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

#endif
#endif
