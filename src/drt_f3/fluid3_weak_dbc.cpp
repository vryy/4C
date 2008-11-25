#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_weak_dbc.H"



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
DRT::ELEMENTS::Fluid3SurfaceWeakDBCInterface* DRT::ELEMENTS::Fluid3SurfaceWeakDBCInterface::Impl(
  DRT::ELEMENTS::Fluid3Surface* f3surf
  )
{
  switch (f3surf->Shape())
  {
  case DRT::Element::quad4:
  {
    static Fluid3SurfaceWeakDBC<DRT::Element::quad4,DRT::Element::hex8>* fsurfq4;

    if(f3surf->parent_->Shape()==DRT::Element::hex8)
    {
      if (fsurfq4==NULL)
        fsurfq4 = new Fluid3SurfaceWeakDBC<DRT::Element::quad4,DRT::Element::hex8>();
    }
    else
    {
      dserror("expected combination quad4/hex8 for surface/parent pair");
    }

    return fsurfq4;
  }
  default:
    dserror("shape %d (%d nodes) not supported by weak DBC", f3surf->Shape(), f3surf->NumNode());
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
DRT::ELEMENTS::Fluid3SurfaceWeakDBC<distype,pdistype>::Fluid3SurfaceWeakDBC()
{
  return;
}

//-----------------------------------------------------------------
//             evaluate implementation for weak dbcs
//-----------------------------------------------------------------
template <DRT::Element::DiscretizationType distype,
          DRT::Element::DiscretizationType pdistype>
int DRT::ELEMENTS::Fluid3SurfaceWeakDBC<distype,pdistype>::EvaluateWeakDBC(
  Fluid3Surface*             surfele       ,
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
  if(surfele->parent_->Shape()!=Fluid3::hex8)
  {
    dserror("Cb up to now only implemented for trilinear elements");
  }
  const double Cb    = 4.0;

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
      
  surfele->parent_->LocationVector(discretization,*plm,*plmowner);

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

    pevelnp_(0,i) = mypvelaf[  fi];
    pevelnp_(1,i) = mypvelaf[1+fi];
    pevelnp_(2,i) = mypvelaf[2+fi];

    peprenp_(  i) = mypvelaf[3+fi];
  }

  // extract node coords
  for(int i=0;i<piel;++i)
  {
    pxyze_(0,i)=surfele->parent_->Nodes()[i]->X()[0];
    pxyze_(1,i)=surfele->parent_->Nodes()[i]->X()[1];
    pxyze_(2,i)=surfele->parent_->Nodes()[i]->X()[2];
  }

  //--------------------------------------------------
  // get material of volume element this surface belongs to
  RCP<MAT::Material> mat = surfele->parent_->Material();

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

  // local surface id
  int surfaceid =surfele->lsurface_;

  // extract node coords
  for(int i=0;i<iel;++i)
  {
    xyze_(0,i)=surfele->Nodes()[i]->X()[0];
    xyze_(1,i)=surfele->Nodes()[i]->X()[1];
    xyze_(2,i)=surfele->Nodes()[i]->X()[2];
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
  default:
    dserror("invalid discretization type for fluid3surface weak DBC evaluation");
  }

  // gaussian points on surface
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //--------------------------------------------------
  // the gausspoints above have to be mapped to the 
  // parent element to be able to evaluate one sided
  // derivatives on the boundary
  Epetra_SerialDenseMatrix pqxg(intpoints.nquad,3);

  switch(surfaceid)
  {
  case 0:
  {
    // t=-1
    /*
                parent               surface

                 r|                    s|                    
                  |                     |                    
             1         2           3         2               
              +-------+             +-------+                
	      |   |   |      s      |   |   |      r  
              |   +-- |  -----      |   +-- |  -----  
              |       |             |       |                
              +-------+             +-------+                
             0         3           0         1               
    */

    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      pqxg(iquad,0)= intpoints.qxg[1][iquad];
      pqxg(iquad,1)= intpoints.qxg[0][iquad];
      pqxg(iquad,2)=-1.0;
    }
    break;
  }
  case 1:
  {
    // s=-1
    /*
                parent               surface

                 t|                    s|                    
                  |                     |                    
             4         5           3         2               
              +-------+             +-------+                
	      |   |   |      r      |   |   |      r  
              |   +-- |  -----      |   +-- |  -----  
              |       |             |       |                
              +-------+             +-------+                
             0         1           0         1               
    */
    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      pqxg(iquad,0)= intpoints.qxg[0][iquad];
      pqxg(iquad,1)=-1.0;
      pqxg(iquad,2)= intpoints.qxg[1][iquad];
    }
    break;
  }
  case 2:
  {
    // r= 1
    /*
                parent               surface

                 t|                    s|                    
                  |                     |                    
             5         6           3         2               
              +-------+             +-------+                
	      |   |   |      s      |   |   |      r  
              |   +-- |  -----      |   +-- |  -----  
              |       |             |       |                
              +-------+             +-------+                
             1         2           0         1               
    */
    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      pqxg(iquad,0)= 1.0;
      pqxg(iquad,1)= intpoints.qxg[0][iquad];
      pqxg(iquad,2)= intpoints.qxg[1][iquad];
    }
    break;
  }
  case 3:
  {
    // s= 1
    /*
                parent               surface

                 t|                    s|                    
                  |                     |                    
             6         7           3         2               
              +-------+             +-------+                
	r     |   |   |             |   |   |      r  
        ----  |   +-- |             |   +-- |  -----  
              |       |             |       |                
              +-------+             +-------+                
             2         3           0         1               
    */
    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      pqxg(iquad,0)=-intpoints.qxg[0][iquad];
      pqxg(iquad,1)= 1.0;
      pqxg(iquad,2)= intpoints.qxg[1][iquad];
    }
    break;
  }
  case 4:
  {
    // r=-1
    /*
                parent               surface

                 s|                    s|                    
                  |                     |                    
             3         7           3         2               
              +-------+             +-------+                
	      |   |   |      t      |   |   |      r  
              |   +-- |  -----      |   +-- |  -----  
              |       |             |       |                
              +-------+             +-------+                
             0         4           0         1               
    */
    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      pqxg(iquad,0)=-1.0;
      pqxg(iquad,1)= intpoints.qxg[1][iquad];
      pqxg(iquad,2)= intpoints.qxg[0][iquad];
    }
    break;
  }
  case 5:
  {
    // t=1
    /*
                parent               surface

                 s|                    s|                    
                  |                     |                    
             7         6           3         2               
              +-------+             +-------+                
	      |   |   |      r      |   |   |      r  
              |   +-- |  -----      |   +-- |  -----  
              |       |             |       |                
              +-------+             +-------+                
             4         5           0         1               
    */
    for (int iquad=0;iquad<intpoints.nquad;++iquad)
    { 
      pqxg(iquad,0)= intpoints.qxg[0][iquad];
      pqxg(iquad,1)= intpoints.qxg[1][iquad];
      pqxg(iquad,2)= 1.0;
    }
    break;
  }
  default:
    dserror("invalid number of surfaces, unable to determine intpoint in parent");
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

    // ------------------------------------------------
    // shape function derivs of boundary element at gausspoint
    DRT::UTILS::shape_function_2D       (funct_,xi,eta,distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv_,xi,eta,distype);

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
	    // evaluate function at current gauss point
	    functionfac(dim) = DRT::UTILS::FunctionManager::Instance().Funct(functnum-1).Evaluate(dim,coordgp.A());
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
    DRT::UTILS::shape_function_3D       (pfunct_,r,s,t,pdistype);
    DRT::UTILS::shape_function_3D_deriv1(pderiv_,r,s,t,pdistype);


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
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", surfele->parent_->Id(), pdet);
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

    //--------------------------------------------------
    // partially integrated pressure term
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
        elemat(vi*4  ,ui*4+3) += fac*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*4+1,ui*4+3) += fac*pfunct_(vi)*pfunct_(ui)*n_(1);
        elemat(vi*4+2,ui*4+3) += fac*pfunct_(vi)*pfunct_(ui)*n_(2);
      }
    }

    for (int vi=0; vi<piel; ++vi) 
    {
      elevec(vi*4    ) -= fac*pfunct_(vi)*n_(0)*prenp_;
      elevec(vi*4 + 1) -= fac*pfunct_(vi)*n_(1)*prenp_;
      elevec(vi*4 + 2) -= fac*pfunct_(vi)*n_(2)*prenp_;
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
    const double timefacnu=fac*2.0*visc*afgdt;

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
      nabla_u_o_n[0]=fac*2.0*visc*
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
        elemat(vi*4+3,ui*4  ) -= fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(0);
        elemat(vi*4+3,ui*4+1) -= fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(1);
        elemat(vi*4+3,ui*4+2) -= fac*gdt*pfunct_(vi)*pfunct_(ui)*n_(2);
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
      elevec(vi*4+3) += fac*pfunct_(vi)*
        ((velintnp_(0)-(*val)[0]*functionfac(0)*curvefac)*n_(0)
         +
         (velintnp_(1)-(*val)[1]*functionfac(1)*curvefac)*n_(1)
         +
         (velintnp_(2)-(*val)[2]*functionfac(2)*curvefac)*n_(2));
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
    const double consistencytimefac=fac*2.0*visc*wd_gamma*afgdt;

    for (int ui=0; ui<piel; ++ui) 
    {
      for (int vi=0; vi<piel; ++vi) 
      {

	elemat(vi*4    ,ui*4    ) -= consistencytimefac*pfunct_(ui)*
          (n_(0)*(    pderxy_(0,vi))+n_(1)*(0.5*pderxy_(1,vi))+n_(2)*(0.5*pderxy_(2,vi)));
        elemat(vi*4    ,ui*4 + 1) -= consistencytimefac*pfunct_(ui)*
          (n_(0)*(0.5*pderxy_(1,vi)));
        elemat(vi*4    ,ui*4 + 2) -= consistencytimefac*pfunct_(ui)*
          (n_(0)*(0.5*pderxy_(2,vi)));
        
        elemat(vi*4 + 1,ui*4    ) -= consistencytimefac*pfunct_(ui)*
          (n_(1)*(0.5*pderxy_(0,vi)));
        elemat(vi*4 + 1,ui*4 + 1) -= consistencytimefac*pfunct_(ui)*
          (n_(0)*(0.5*pderxy_(0,vi))+n_(1)*(    pderxy_(1,vi))+n_(2)*(0.5*pderxy_(2,vi)));
        elemat(vi*4 + 1,ui*4 + 2) -= consistencytimefac*pfunct_(ui)*
          (n_(1)*(0.5*pderxy_(2,vi)));

	elemat(vi*4 + 2,ui*4    ) -= consistencytimefac*pfunct_(ui)*
          (n_(2)*(0.5*pderxy_(0,vi)));
        elemat(vi*4 + 2,ui*4 + 1) -= consistencytimefac*pfunct_(ui)*
          (n_(2)*(0.5*pderxy_(1,vi)));
        elemat(vi*4 + 2,ui*4 + 2) -= consistencytimefac*pfunct_(ui)*
          (n_(0)*(0.5*pderxy_(0,vi))+n_(1)*(0.5*pderxy_(1,vi))+n_(2)*(    pderxy_(2,vi)));
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
      elevec(vi*4    ) += fac*2.0*visc*wd_gamma*(
        n_(0)*(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac)*(    pderxy_(0,vi))
        +
        n_(1)*(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac)*(0.5*pderxy_(1,vi))
        +
        n_(2)*(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac)*(0.5*pderxy_(2,vi))
        +
        n_(0)*(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac)*(0.5*pderxy_(1,vi))
        +
        n_(0)*(velintaf_(2)-(*val)[2]*functionfac(2)*curvefac)*(0.5*pderxy_(2,vi)));

      elevec(vi*4 + 1) += fac*2.0*visc*wd_gamma*(
        n_(1)*(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac)*(0.5*pderxy_(0,vi))
        +
        n_(0)*(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac)*(0.5*pderxy_(0,vi))
        +
        n_(1)*(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac)*(    pderxy_(1,vi))
        +
        n_(2)*(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac)*(0.5*pderxy_(2,vi))
        +
        n_(1)*(velintaf_(2)-(*val)[2]*functionfac(2)*curvefac)*(0.5*pderxy_(2,vi)));

      elevec(vi*4 + 2) += fac*2.0*visc*wd_gamma*(
        n_(2)*(velintaf_(0)-(*val)[0]*functionfac(0)*curvefac)*(0.5*pderxy_(0,vi))
        +
        n_(2)*(velintaf_(1)-(*val)[1]*functionfac(1)*curvefac)*(0.5*pderxy_(1,vi))
        +
        n_(0)*(velintaf_(2)-(*val)[2]*functionfac(2)*curvefac)*(0.5*pderxy_(0,vi))
        +
        n_(1)*(velintaf_(2)-(*val)[2]*functionfac(2)*curvefac)*(0.5*pderxy_(1,vi))
        +
        n_(2)*(velintaf_(2)-(*val)[2]*functionfac(2)*curvefac)*(    pderxy_(2,vi)));
    }

    //--------------------------------------------------
    // adjoint consistency term, convective part (only on inflow)
    
    double flux=velintaf_(0)*n_(0)+velintaf_(1)*n_(1)+velintaf_(2)*n_(2);

    if(flux<0)
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
          elemat(vi*4    ,ui*4    ) -= timefac*pfunct_(vi)*n_(0)*
            (velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);
          elemat(vi*4    ,ui*4 + 1) -= timefac*pfunct_(vi)*n_(1)*
            (velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);
          elemat(vi*4    ,ui*4 + 2) -= timefac*pfunct_(vi)*n_(2)*
            (velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);

          elemat(vi*4 + 1,ui*4    ) -= timefac*pfunct_(vi)*n_(0)*
            (velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);
          elemat(vi*4 + 1,ui*4 + 1) -= timefac*pfunct_(vi)*n_(1)*
            (velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);
          elemat(vi*4 + 1,ui*4 + 2) -= timefac*pfunct_(vi)*n_(2)*
            (velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);
          
          elemat(vi*4 + 2,ui*4    ) -= timefac*pfunct_(vi)*n_(0)*
            (velintaf_(2)-(*val)[2]*functionfac(2)*curvefac);
          elemat(vi*4 + 2,ui*4 + 1) -= timefac*pfunct_(vi)*n_(1)*
            (velintaf_(2)-(*val)[2]*functionfac(2)*curvefac);
          elemat(vi*4 + 2,ui*4 + 2) -= timefac*pfunct_(vi)*n_(2)*
            (velintaf_(2)-(*val)[2]*functionfac(2)*curvefac);
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
          elemat(vi*4    ,ui*4    ) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
          elemat(vi*4 + 1,ui*4 + 1) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
          elemat(vi*4 + 2,ui*4 + 2) -= fluxtimefac*pfunct_(ui)*pfunct_(vi);
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
      const double fluxfac =fac*flux;
      
      for (int vi=0; vi<piel; ++vi) 
      {
        elevec(vi*4    ) += fluxfac*pfunct_(vi)*
          (velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);
        elevec(vi*4 + 1) += fluxfac*pfunct_(vi)*
          (velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);
        elevec(vi*4 + 2) += fluxfac*pfunct_(vi)*
          (velintaf_(2)-(*val)[2]*functionfac(2)*curvefac);
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

      const double penaltytimefac=afgdt*Cb*visc/h*fac;
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

      const double penaltyfac=Cb*visc/h*fac;
      for (int vi=0; vi<piel; ++vi) 
      {
        elevec(vi*4    ) -= penaltyfac*pfunct_(vi)*
          (velintaf_(0)-(*val)[0]*functionfac(0)*curvefac);
        elevec(vi*4 + 1) -= penaltyfac*pfunct_(vi)*
          (velintaf_(1)-(*val)[1]*functionfac(1)*curvefac);
        elevec(vi*4 + 2) -= penaltyfac*pfunct_(vi)*
          (velintaf_(2)-(*val)[2]*functionfac(2)*curvefac);
      }
    } // end if flux<0, i.e.

  } // end gaussloop

  return 0;
}

#endif
#endif
