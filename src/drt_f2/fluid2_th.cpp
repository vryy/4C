/*----------------------------------------------------------------------*/
/*!
\file fluid2_th.cpp

\brief Internal implementation of Fluid2 element (Taylor Hood Q2Q1)

<pre>
Maintainer: Tobias Wiesner
wiesner@lnm.mw.tum.de
Created on: Jun 3, 2009
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2_th.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_inpar/inpar_fluid.H"

#include "../drt_lib/drt_discret.H"

#include <Epetra_SerialDenseSolver.h>
#include <Epetra_Vector.h>
#include <Epetra_LAPACK.h>


/*----------------------------------------------------------------------*
 * Interface class for Taylor-Hood element
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2THInterface* DRT::ELEMENTS::Fluid2THInterface::Impl(DRT::ELEMENTS::Fluid2* f2)
{
  switch (f2->Shape())
  {
    case DRT::Element::quad4:
    {
      static Fluid2TH<DRT::Element::quad4>* fq4;
      if (fq4==NULL)
        fq4 = new Fluid2TH<DRT::Element::quad4>;
      return fq4;
    }
    case DRT::Element::quad9:
    {
      static Fluid2TH<DRT::Element::quad9>* fq9;
      if (fq9==NULL)
        fq9 = new Fluid2TH<DRT::Element::quad9>;
      return fq9;
    }
    default:
      dserror("shape %d (%d nodes) not supported", f2->Shape(), f2->NumNode());
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/// constructor
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Fluid2TH<distype>::Fluid2TH()
: vart_(),
xyze_(),
edeadng_(),
funct_(),
densfunct_(),
deriv_(),
deriv2_(),
xjm_(),
xji_(),
vderxy_(),
mderxy_(),
fsvderxy_(),
vderxy2_(),
derxy_(),
densderxy_(),
derxy2_(),
bodyforce_(),
velint_(),
velints_(),
fsvelint_(),
gradp_(),
tau_(),
viscs2_(),
conv_c_(),
mdiv_(),
vdiv_(),
rhsmom_(),
conv_old_(),
visc_old_(),
res_old_(),
conv_resM_(),
xder2_()
{
}

/*!
* basic element call for projection methods
* for calculation of mass matrix, gradop and lumped mass matrix
* input: dispnp for ALE case
*
* @param ele Element
* @param params parameter list
* @param discretization
* @param lm vector with ids
* @param elemat1_epetra gradop (o)
* @param elemat2_epetra mass matrix (o)
* @param elevec1_epetra lumped mass matrix as vector (o)
* @return 0
*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid2TH<distype>::CalcGradPAndMassMatrix(Fluid2*	ele, ParameterList&				params, DRT::Discretization&		discretization, vector<int>&				lm,  Epetra_SerialDenseMatrix&	elemat1_epetra, Epetra_SerialDenseMatrix&	elemat2_epetra, Epetra_SerialDenseVector& 	elevec1_epetra)
{
  LINALG::Matrix<idofs,idofs> elemat(elemat1_epetra.A(),true);
  LINALG::Matrix<idofs, 1> elelmass(elevec1_epetra.A(),true);

  int* ndofs = NULL;
  switch(distype)
  {
    // for quad elements: array with startindices per node
    // for quad4 e.g. only first four entries (3 DOFs per node)
    // for quad9 additional 5 entries with 2 DOFs.
    case DRT::Element::quad9:
    case DRT::Element::quad8:
    case DRT::Element::quad4:
    {
      static int ndofarrayquad[9] = {0,3,6,9,12,14,16,18,20};
      ndofs = &ndofarrayquad[0];
    }
    break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      static int ndofarraytri[6] = {0,3,6,9,11,13};
      ndofs = &ndofarraytri[0];
    }
    break;
    default:
      dserror("distype not supported");
  }

  // store node coords
  DRT::Node** nodes = ele->Nodes();
  for(int inode=0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // in ALE case check displacement
  RCP<const Epetra_Vector> dispnp;
  vector<double> mydispnp;
  LINALG::Matrix<2, iel> edispnp;

  if (ele->is_ale_)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

    // assign grid displacement to element arrays
    for (int i=0;i<iel;++i)
    {
      edispnp(0,i) = mydispnp[ndofs[i]];
      edispnp(1,i) = mydispnp[ndofs[i]+1];
    }
    xyze_ += edispnp;
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);

  // integration loop over all gauß points
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    // coords of current gauß point
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    // for tri elements necessary
    pfunct_.PutScalar(0.0);
    pderiv_.PutScalar(0.0);

    // shape functions and derivatives
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
        DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
        DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);
        DRT::UTILS::shape_function_2D(pfunct_,e1,e2,DRT::Element::quad4);	// shape function for pressure (quad9)
        break;
      default: dserror("tri elements not yet implemented");
    }


    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);
    if(det<0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(),det);

    const double fac = intpoints.qwgt[iquad]*det;

    // global first derivative
    derxy_.Multiply(xji_,deriv_);

    // compute gradp
    for(int i=0;i<iel;i++)
    {
      for(int p=0; p<4; p++) // only pressure dofs //TODO no tri elments supported.
      {
        elemat(ndofs[i],3*p+2)   -= fac*pfunct_(p)*derxy_(0,i);
        elemat(ndofs[i]+1,3*p+2) -= fac*pfunct_(p)*derxy_(1,i);
      }
    }

    // mass matrix
    for(int i=0; i<iel; i++)
    {
      for(int j=0; j<iel; j++)
      {
        double temp = fac*funct_(i)*funct_(j);
        elemat(ndofs[i],ndofs[j]) += temp;
        elemat(ndofs[i]+1,ndofs[j]+1) += temp;
      }
    }

  }	// <- end loop over integration points

  // Lumped mass matrix
  for(int i=0; i<iel; i++)
  {
    int row = ndofs[i];
    elelmass(row) = 0;
    elelmass(row+1) = 0;
    for(int j=0; j<iel; j++)
    {
      int col = ndofs[j];
      elelmass(row) += elemat(row,col);
      elelmass(row+1) += elemat(row+1,col+1);
    }
  }

  return 0;
}

/*!
* basic element call for projection methods
* for calculation of impulse equation (implicit)
*
* @param ele Element
* @param params parameter list
* @param discretization
* @param lm vector with ids
* @param elemat1_epetra stiffmat (o)
* @param elevec1_epetra rhs (o)
* @return 0
*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid2TH<distype>::CalcImpulseEqnImplicit(Fluid2* ele,
    ParameterList& params,
    DRT::Discretization& discretization,
    vector<int> lm,
    Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseVector& elevec1_epetra,
    RefCountPtr< MAT::Material >   material/*,
		_MATERIAL *   material*/)
{
  LINALG::Matrix<idofs,idofs> elestiff(elemat1_epetra.A(),true);
  LINALG::Matrix<idofs, 1> eleforce(elevec1_epetra.A(),true);

  int* ndofs = NULL;
  switch(distype)
  {
    case DRT::Element::quad9:
    case DRT::Element::quad8:
    case DRT::Element::quad4:
    {
      static int ndofarrayquad[9] = {0,3,6,9,12,14,16,18,20};
      ndofs = &ndofarrayquad[0];
    }
    break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      static int ndofarraytri[6] = {0,3,6,9,11,13};
      ndofs = &ndofarraytri[0];
    }
    break;
    default:
      dserror("distype not supported");
  }

  // get current state vectors from discretization
  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  RCP<const Epetra_Vector> vedenp = discretization.GetState("scanp");
  RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  RCP<const Epetra_Vector> dispnp;
  RCP<const Epetra_Vector> gridv;

  vector<double> myvelnp(lm.size());
  vector<double> myvedenp(lm.size());
  vector<double> myhist(lm.size());
  vector<double> mydispnp;
  vector<double> mygridv;
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

  // ALE case
  if(ele->is_ale_)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vector 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    gridv = discretization.GetState("gridv");
    if (gridv==null) dserror("Cannot get state vector 'gridv'");
    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
  }

  // get element vectors
  LINALG::Matrix<iel,1>   epren;	// pressure time n (for quad elements 4 and for tri elements 3 entries; rest is set to zero)
  LINALG::Matrix<2,iel>	evelnp;	// velocity n+1
  LINALG::Matrix<iel,1>	edensnp;// density n+1
  LINALG::Matrix<2,iel>	emhist;
  LINALG::Matrix<iel,1>	echist;
  LINALG::Matrix<2,iel>	edispnp;
  LINALG::Matrix<2,iel>	egridv;

  // zero out element vectors
  epren.Clear();
  evelnp.Clear();
  edensnp.PutScalar(1.0);		// density
  emhist.Clear();
  echist.Clear();
  edispnp.Clear();
  egridv.Clear();

  // fill element vectors
  for(int i=0; i<iel; i++)
  {
    evelnp(0,i) = myvelnp[ndofs[i]];
    evelnp(1,i) = myvelnp[ndofs[i]+1];
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
      {
        if(i<4)
        {
          epren(i) 	  = myvelnp[ndofs[i]+2];
          edensnp(i)  = myvedenp[ndofs[i]+2];
        }
      }
      break;
      case DRT::Element::tri6:
      case DRT::Element::tri3:
      {
        if(i<3)
        {
          epren(i)	= myvelnp[ndofs[i]+2];
          edensnp(i)  = myvedenp[ndofs[i]+2];
        }
      }
      break;
      default:
        dserror("distype not supported");
    }

    emhist(0,i) = myhist[ndofs[i]];
    emhist(1,i) = myhist[ndofs[i]+1];

    // ALE part
    if(ele->is_ale_)
    {
      edispnp(0,i) = mydispnp[ndofs[i]];
      edispnp(1,i) = mydispnp[ndofs[i]+1];
      egridv(0,i) = mygridv[ndofs[i]];
      egridv(1,i) = mygridv[ndofs[i]+1];
    }
  }

  // no fine-scale viscosity
  LINALG::Matrix<2,iel> fsevelnp;
  for(int i=0;i<iel;i++)
  {
    fsevelnp(0,i) = 0.0;
    fsevelnp(1,i) = 0.0;
  }
  Fluid2::FineSubgridVisc fssgv = Fluid2::no_fssgv;

  // get control parameters for time integration
  const double time = params.get<double>("total time",-1.0);

  bool newton = false;
  // bool loma = false;
  string newtonstr = params.get<string>("Linearisation");
  // string lomastr = params.get<string>("low-Mach-number solver");
  if(newtonstr=="Newton") newton = true;
  // if(lomastr=="Yes") loma = true;
  INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  // stabilization
  ParameterList& stablist = params.sublist("STABILIZATION");

  Fluid2::TauType whichtau = Fluid2::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");
    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid2::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid2::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid2::codina;
    }
  }

  bool higher_order_ele = ele->isHigherOrderElement(distype);
  if(stablist.get<string>("STABTYPE") == "inconsistent") higher_order_ele = false;

  const double dt = params.get<double>("dt");

  // One-step-Theta: timefac = theta*dt
  // BDF2:           timefac = 2/3 * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // --------------------------------------------------
  // set parameters for classical turbulence models
  // --------------------------------------------------
  //ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");

  // initialise the Smagorinsky constant Cs to zero
  double Cs            = 0.0;
  double visceff       = 0.0;

  // the default action is no model
  Fluid2::TurbModelAction turb_mod_action = Fluid2::no_model;

  // save coordinates of nodes
  DRT::Node** nodes = ele->Nodes();
  for(int inode = 0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // ALE case: add displacement
  if(ele->is_ale_) xyze_+=edispnp;

  // dead load in element nodes
  BodyForce(ele,time);

  // check here, if we really have a fluid !!
  if( material->MaterialType() != INPAR::MAT::m_fluid
      && material->MaterialType() != INPAR::MAT::m_carreauyasuda
      && material->MaterialType() != INPAR::MAT::m_modpowerlaw) dserror("Material law is not a fluid");

  // viscosity
  double visc = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
    visc = actmat->Viscosity();
  }

  Caltau(ele,evelnp,fsevelnp,edensnp,whichtau,material,visc,timefac,dt,turb_mod_action,Cs,visceff,fssgv);

  // in case of viscous stabilization decide whether to use GLS or USFEM
  /*double vstabfac= 0.0;
	if (vstab == Fluid2::viscous_stab_usfem || vstab == Fluid2::viscous_stab_usfem_only_rhs)
	{
		vstabfac =  1.0;
	}
	else if(vstab == Fluid2::viscous_stab_gls || vstab == Fluid2::viscous_stab_gls_only_rhs)
	{
		vstabfac = -1.0;
	}*/

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
#ifdef DEBUG
  if(ele->gaussrule_!=DRT::UTILS::intrule_quad_9point)
    cout << "WARNING: no 9point gaussrule for quadrature! element matrices may be nontrivial disturbed!" << endl;
#endif

  // integration loop for gauss points
  for(int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    pfunct_.PutScalar(0.0);
    pderiv_.PutScalar(0.0);

    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
        DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
        DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);
        DRT::UTILS::shape_function_2D(pfunct_,e1,e2,DRT::Element::quad4);	// shape functions for pressure (quad9!) TODO no tri elements supported
        DRT::UTILS::shape_function_2D_deriv1(pderiv_,e1,e2,DRT::Element::quad4);
        break;
      default: dserror("not yet implemented");
    }

    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);
    if(det<0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(),det);

    const double fac = intpoints.qwgt[iquad]*det;

    derxy_.Multiply(xji_,deriv_);
    pderxy_.Multiply(xji_,pderiv_);

    vderxy_.MultiplyNT(evelnp,derxy_);

    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    // (inverse-)density-weighted shape functions and global derivatives
    for (int inode=0; inode<iel; inode++)
    {
      densfunct_(inode) = edensnp(inode)*funct_(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    if (higher_order_ele)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,e1,e2,distype);
      // initialize everything
      LINALG::Matrix<3,3> bm;

      // calculate elements of jacobian_bar matrix
      bm(0,0) =                     xjm_(0,0)*xjm_(0,0);
      bm(0,1) =                     xjm_(0,1)*xjm_(0,1);
      bm(0,2) =                 2.0*xjm_(0,0)*xjm_(0,1);

      bm(1,0) =                     xjm_(1,0)*xjm_(1,0);
      bm(1,1) =                     xjm_(1,1)*xjm_(1,1);
      bm(1,2) =                 2.0*xjm_(1,1)*xjm_(1,0);

      bm(2,0) =                     xjm_(0,0)*xjm_(1,0);
      bm(2,1) =                     xjm_(0,1)*xjm_(1,1);
      bm(2,2) = xjm_(0,0)*xjm_(1,1)+xjm_(0,1)*xjm_(1,0);

      xder2_.MultiplyNT(deriv2_, xyze_);
      derxy2_.Multiply(-1.0, xder2_, derxy_);
      derxy2_ += deriv2_;

      LINALG::FixedSizeSerialDenseSolver<3,3,iel> solver_bm;
      solver_bm.SetMatrix(bm);                // A = bm
      // X == B is not a problem, derxy2_ will contain the solution.
      solver_bm.SetVectors(derxy2_,derxy2_);  // X = B = derxy2_
      // error code
      int ierr = 0;

      // Factor. Calling this directly is not necessary, only to find
      // out where an error came from.
      ierr = solver_bm.Factor();
      if (ierr!=0)
        dserror("Unable to perform LU factorisation during computation of derxy2");

      solver_bm.Solve();                      // Solve A*X = B
      if (ierr!=0)
        dserror("Unable to perform backward substitution after factorisation of jacobian");

      // calculate 2nd velocity derivatives at integration point
      vderxy2_.MultiplyNT(evelnp, derxy2_);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get momentum (n+g,i) at integration point
    velint_.Multiply(evelnp, densfunct_);

    // history data (n,i) am Integrationspunkt
    histmom_.Multiply(emhist, densfunct_);
    histcon_ = funct_.Dot(echist);

    // get convective velocity at integration point
    // We handle the ale case very implicitely here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_ = velint_;
    if (ele->is_ale_) convvelint_.Multiply(-1.0, egridv, densfunct_, 1.0);

    // get pressure gradient at integration point
    gradp_.Multiply(pderxy_, epren);

    const double timefacfac = timefac * fac;

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_.Multiply(edeadng_, densfunct_);
    // get density at integration point
    double dens = funct_.Dot(edensnp);

    /*------------------------ evaluate rhs vectors at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    rhsmom_(0) = histmom_(0) + bodyforce_(0)*timefac;
    rhsmom_(1) = histmom_(1) + bodyforce_(1)*timefac;
    rhscon_ = histcon_ - dens;

    /*----------------- get numerical representation of single operators ---*/
    /* Convective term  u_old * grad u_old: */
    conv_old_.Multiply(vderxy_, convvelint_);

    /* Viscous term  div epsilon(u_old) */
    if (higher_order_ele)
    {
      visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,2));
      visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,2));
    }
    else
    {
      visc_old_.Clear();
    }

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
		          with  N .. form function matrix                                   */
    conv_c_.MultiplyTN(derxy_, convvelint_);

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
		              1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
		              - |                              |         N_y .. y-line of N
		             2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
		                \                              /                            */

      for (int i=0; i<iel; ++i) viscs2_(0,i) = 0.5 * (2.0 * derxy2_(0,i) + derxy2_(1,i));
      for (int i=0; i<iel; ++i) viscs2_(1,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(2,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(3,i) = 0.5 * (derxy2_(0,i) + 2.0 * derxy2_(1,i));
    }
    else
    {
      viscs2_.Clear();
    }

    /* momentum and velocity divergence: */
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (physicaltype == INPAR::FLUID::loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);
    //if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    {
      ////////////////// GALERKIN ANTEIL ///////////////////////

      /* Konvektionsterm (u*grad(u),v) + Anteil Massenmatrix*/
      for(int ui=0; ui<iel; ++ui)
      {
        const int tui = ndofs[ui];
        const double v = fac*densfunct_(ui)+timefacfac*conv_c_(ui);	//13.11 keine Konvektion auf linker Seite
        for(int vi=0; vi<iel; ++vi)
        {
          const int tvi = ndofs[vi];
          double v2 = v*funct_(vi) ;
          elestiff(tvi,     tui    ) += v2;
          elestiff(tvi + 1, tui + 1) += v2;
        }
      }


      /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
      const double viscefftimefacfac = visceff*timefacfac;
      for (int ui=0; ui<iel; ++ui)
      {
        const int tui  = ndofs[ui];
        const int tuip = tui+1;
        for (int vi=0; vi<iel; ++vi)
        {
          const int tvi  = ndofs[vi];
          const int tvip = tvi+1;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);

          elestiff(tvi,  tui ) += viscefftimefacfac*(2.0*derxy_0ui_0vi
              +
              derxy_1ui_1vi) ;
          elestiff(tvi,  tuip) += viscefftimefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          elestiff(tvip, tui ) += viscefftimefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          elestiff(tvip, tuip) += viscefftimefacfac*(derxy_0ui_0vi
              +
              2.0*derxy_1ui_1vi) ;

        }
      }

      // reaktive Konvektionsterme (ok)
      if (newton)
      {
        for (int vi=0; vi<iel; ++vi)
        {
          const int tvi  = ndofs[vi];
          const int tvip = tvi+1;
          const double v = timefacfac*funct_(vi);
          for (int ui=0; ui<iel; ++ui)
          {
            const int tui  = ndofs[ui];
            const int tuip = tui+1;
            const double v2 = v*densfunct_(ui);

            /*  convection, reactive part

			               /                           \
			               |  /          \   n+1       |
			               | | Du o nabla | u     , v  |
			               |  \          /   (i)       |
			               \                           /
             */
            elestiff(tvi,  tui ) += v2*vderxy_(0, 0) ;
            elestiff(tvi,  tuip) += v2*vderxy_(0, 1) ;
            elestiff(tvip, tui ) += v2*vderxy_(1, 0) ;
            elestiff(tvip, tuip) += v2*vderxy_(1, 1) ;
          }
        }
      }

      /* inertia */
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* inertia */
        const double v = -fac*funct_(vi);
        eleforce(tvi    ) += v*velint_(0) ;
        eleforce(tvi + 1) += v*velint_(1) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* convection */
        double v = -timefacfac*funct_(vi);
        eleforce(tvi    ) += v*(convvelint_(0)*vderxy_(0, 0)
            +
            convvelint_(1)*vderxy_(0, 1)) ;
        eleforce(tvi + 1) += v*(convvelint_(0)*vderxy_(1, 0)
            +
            convvelint_(1)*vderxy_(1, 1)) ;
      }

      // viscosity
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* viscosity */
        eleforce(tvi    ) -= visceff*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
            +
            derxy_(1, vi)*vderxy_(0, 1)
            +
            derxy_(1, vi)*vderxy_(1, 0)) ;
        eleforce(tvi + 1) -= visceff*timefacfac*(derxy_(0, vi)*vderxy_(0, 1)
            +
            derxy_(0, vi)*vderxy_(1, 0)
            +
            2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      // source term
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eleforce(tvi    ) += v*rhsmom_(0) ;
        eleforce(tvi + 1) += v*rhsmom_(1) ;
      }
    } // <- Elementanteile (Galerkin + Stabiliserungsterme)
  } // <- Ende schleife über alle Integrationspunkte
  return 0;
}

/*!
* basic element call for projection methods
* for calculation of impulse equation (semi-implicit)
*
* @param ele Element
* @param params parameter list
* @param discretization
* @param lm vector with ids
* @param elemat1_epetra stiffmat (o)
* @param elevec1_epetra rhs (o)
* @return 0
*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid2TH<distype>::CalcImpulseEqnSemiImplicit(Fluid2* ele,
    ParameterList& params,
    DRT::Discretization& discretization,
    vector<int> lm,
    Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseVector& elevec1_epetra,
    RefCountPtr< MAT::Material >   material
)
{
  LINALG::Matrix<idofs,idofs> elestiff(elemat1_epetra.A(),true);
  LINALG::Matrix<idofs, 1> eleforce(elevec1_epetra.A(),true);

  int* ndofs = NULL;
  switch(distype)
  {
    case DRT::Element::quad9:
    case DRT::Element::quad8:
    case DRT::Element::quad4:
    {
      static int ndofarrayquad[9] = {0,3,6,9,12,14,16,18,20};
      ndofs = &ndofarrayquad[0];
    }
    break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      static int ndofarraytri[6] = {0,3,6,9,11,13};
      ndofs = &ndofarraytri[0];
    }
    break;
    default:
      dserror("distype not supported");
  }

  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  RCP<const Epetra_Vector> veln = discretization.GetState("veln");
  //RCP<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  RCP<const Epetra_Vector> vedenp = discretization.GetState("scanp");
  RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  RCP<const Epetra_Vector> dispnp;
  RCP<const Epetra_Vector> gridv;

  vector<double> myvelnp(lm.size());
  vector<double> myveln(lm.size());
  vector<double> myvedenp(lm.size());
  vector<double> myhist(lm.size());
  vector<double> mydispnp;
  vector<double> mygridv;
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  DRT::UTILS::ExtractMyValues(*veln,myveln,lm);
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

  // ALE case
  if(ele->is_ale_)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    gridv = discretization.GetState("gridv");
    if (gridv==null) dserror("Cannot get state vectors 'gridv'");
    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
  }

  LINALG::Matrix<iel,1>   epren;
  LINALG::Matrix<2,iel>	evelnp;
  LINALG::Matrix<2,iel>	eveln;
  LINALG::Matrix<iel,1>	edensnp;
  LINALG::Matrix<2,iel>	emhist;
  LINALG::Matrix<iel,1>	echist;
  LINALG::Matrix<2,iel>	edispnp;
  LINALG::Matrix<2,iel>	egridv;

  epren.Clear();
  evelnp.Clear();
  eveln.Clear();
  edensnp.PutScalar(1.0);
  emhist.Clear();
  echist.Clear();
  edispnp.Clear();
  egridv.Clear();

  // Elementvektoren füllen
  for(int i=0; i<iel; i++)
  {
    evelnp(0,i) = myvelnp[ndofs[i]];
    evelnp(1,i) = myvelnp[ndofs[i]+1];
    eveln(0,i) = myveln[ndofs[i]];
    eveln(1,i) = myveln[ndofs[i]+1];
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
      {
        if(i<4)
        {
          epren(i) 	  = myvelnp[ndofs[i]+2];
          edensnp(i)  = myvedenp[ndofs[i]+2];
        }
      }
      break;
      case DRT::Element::tri6:
      case DRT::Element::tri3:
      {
        if(i<3)
        {
          epren(i)	= myvelnp[ndofs[i]+2];
          edensnp(i)  = myvedenp[ndofs[i]+2];
        }
      }
      break;
      default:
        dserror("distype not supported");
    }

    emhist(0,i) = myhist[ndofs[i]];
    emhist(1,i) = myhist[ndofs[i]+1];

    // ALE part
    if(ele->is_ale_)
    {
      edispnp(0,i) = mydispnp[ndofs[i]];
      edispnp(1,i) = mydispnp[ndofs[i]+1];
      egridv(0,i) = mygridv[ndofs[i]];
      egridv(1,i) = mygridv[ndofs[i]+1];
    }
  }

  // no fine-scale viscosity
  LINALG::Matrix<2,iel> fsevelnp;
  for(int i=0;i<iel;i++)
  {
    fsevelnp(0,i) = 0.0;
    fsevelnp(1,i) = 0.0;
  }
  Fluid2::FineSubgridVisc fssgv = Fluid2::no_fssgv;

  const double time = params.get<double>("total time",-1.0);

  bool newton = false;
  // bool loma = false;
  string newtonstr = params.get<string>("Linearisation");
  //string lomastr = params.get<string>("low-Mach-number solver");
  if(newtonstr=="Newton") newton = true;
  // if(lomastr=="Yes") loma = true;
  INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  ParameterList& stablist = params.sublist("STABILIZATION");

  Fluid2::TauType whichtau = Fluid2::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");
    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid2::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid2::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid2::codina;
    }
  }

  bool higher_order_ele = ele->isHigherOrderElement(distype);
  if(stablist.get<string>("STABTYPE") == "inconsistent") higher_order_ele = false;

  const double dt = params.get<double>("dt");

  // One-step-Theta: timefac = theta*dt
  // BDF2:           timefac = 2/3 * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // --------------------------------------------------
  // set parameters for classical turbulence models
  // --------------------------------------------------
  //ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");

  // initialise the Smagorinsky constant Cs to zero
  double Cs            = 0.0;
  double visceff       = 0.0;

  // the default action is no model
  Fluid2::TurbModelAction turb_mod_action = Fluid2::no_model;

  // save node coordinates
  DRT::Node** nodes = ele->Nodes();
  for(int inode = 0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // ALE case: add displacements
  if(ele->is_ale_) xyze_+=edispnp;

  // dead load in Elementknoten
  BodyForce(ele,time);

  // check here, if we really have a fluid !!
  if( material->MaterialType() != INPAR::MAT::m_fluid
      && material->MaterialType() != INPAR::MAT::m_carreauyasuda
      && material->MaterialType() != INPAR::MAT::m_modpowerlaw) dserror("Material law is not a fluid");

  // get viscosity
  double visc = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
    visc = actmat->Viscosity();
  }

  Caltau(ele,evelnp,fsevelnp,edensnp,whichtau,material,visc,timefac,dt,turb_mod_action,Cs,visceff,fssgv);

  // in case of viscous stabilization decide whether to use GLS or USFEM
  /*double vstabfac= 0.0;
	if (vstab == Fluid2::viscous_stab_usfem || vstab == Fluid2::viscous_stab_usfem_only_rhs)
	{
		vstabfac =  1.0;
	}
	else if(vstab == Fluid2::viscous_stab_gls || vstab == Fluid2::viscous_stab_gls_only_rhs)
	{
		vstabfac = -1.0;
	}*/

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
#ifdef DEBUG
  if(ele->gaussrule_!=DRT::UTILS::intrule_quad_9point)
    cout << "WARNING: no 9point gaussrule for quadrature! element matrices may be nontrivial disturbed!" << endl;
#endif

  // integration loop over all gauss points
  for(int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    pfunct_.PutScalar(0.0);
    pderiv_.PutScalar(0.0);

    // shape functions and their derivatives
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
        DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
        DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);
        DRT::UTILS::shape_function_2D(pfunct_,e1,e2,DRT::Element::quad4);
        DRT::UTILS::shape_function_2D_deriv1(pderiv_,e1,e2,DRT::Element::quad4);
        break;
      default: dserror("not yet implemented");
    }

    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);
    if(det<0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(),det);

    const double fac = intpoints.qwgt[iquad]*det;

    derxy_.Multiply(xji_,deriv_);
    pderxy_.Multiply(xji_,pderiv_);

    vderxy_.MultiplyNT(evelnp,derxy_);

    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    // (inverse-)density-weighted shape functions and global derivatives
    for (int inode=0; inode<iel; inode++)
    {
      densfunct_(inode) = edensnp(inode)*funct_(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    if (higher_order_ele)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,e1,e2,distype);
      // initialize everything
      LINALG::Matrix<3,3> bm;

      // calculate elements of jacobian_bar matrix
      bm(0,0) =                     xjm_(0,0)*xjm_(0,0);
      bm(0,1) =                     xjm_(0,1)*xjm_(0,1);
      bm(0,2) =                 2.0*xjm_(0,0)*xjm_(0,1);

      bm(1,0) =                     xjm_(1,0)*xjm_(1,0);
      bm(1,1) =                     xjm_(1,1)*xjm_(1,1);
      bm(1,2) =                 2.0*xjm_(1,1)*xjm_(1,0);

      bm(2,0) =                     xjm_(0,0)*xjm_(1,0);
      bm(2,1) =                     xjm_(0,1)*xjm_(1,1);
      bm(2,2) = xjm_(0,0)*xjm_(1,1)+xjm_(0,1)*xjm_(1,0);

      xder2_.MultiplyNT(deriv2_, xyze_);
      derxy2_.Multiply(-1.0, xder2_, derxy_);
      derxy2_ += deriv2_;

      LINALG::FixedSizeSerialDenseSolver<3,3,iel> solver_bm;
      solver_bm.SetMatrix(bm);                // A = bm
      // X == B is not a problem, derxy2_ will contain the solution.
      solver_bm.SetVectors(derxy2_,derxy2_);  // X = B = derxy2_
      // error code
      int ierr = 0;

      // Factor. Calling this directly is not necessary, only to find
      // out where an error came from.
      ierr = solver_bm.Factor();
      if (ierr!=0)
        dserror("Unable to perform LU factorisation during computation of derxy2");

      solver_bm.Solve();                      // Solve A*X = B
      if (ierr!=0)
        dserror("Unable to perform backward substitution after factorisation of jacobian");

      // calculate 2nd velocity derivatives at integration point
      vderxy2_.MultiplyNT(evelnp, derxy2_);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get momentum (n+g,i) at integration point
    velint_.Multiply(evelnp, densfunct_);

    // get momentum (n,i) at integration point
    velints_.Multiply(eveln, densfunct_);

    // history data (n,i) am Integrationspunkt
    histmom_.Multiply(emhist, densfunct_);
    histcon_ = funct_.Dot(echist);

    // get convective velocity at integration point
    // We handle the ale case very implicitely here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_ = velints_;	// here's the difference!
    if (ele->is_ale_) convvelint_.Multiply(-1.0, egridv, densfunct_, 1.0);

    // get pressure gradient at integration point
    gradp_.Multiply(pderxy_, epren);

    //--------------------------------------------------------------
    // perform integration for entire matrix and rhs
    //--------------------------------------------------------------

    const double timefacfac = timefac * fac;

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_.Multiply(edeadng_, densfunct_);
    // get density at integration point
    double dens = funct_.Dot(edensnp);

    /*------------------------ evaluate rhs vectors at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    rhsmom_(0) = histmom_(0) + bodyforce_(0)*timefac;
    rhsmom_(1) = histmom_(1) + bodyforce_(1)*timefac;
    rhscon_ = histcon_ - dens;

    /*----------------- get numerical representation of single operators ---*/
    /* Convective term  u_old * grad u_old: */
    conv_old_.Multiply(vderxy_, convvelint_);

    /* Viscous term  div epsilon(u_old) */
    if (higher_order_ele)
    {
      visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,2));
      visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,2));
    }
    else
    {
      visc_old_.Clear();
    }

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
				  with  N .. form function matrix                                   */
    conv_c_.MultiplyTN(derxy_, convvelint_);

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
					  1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
					  - |                              |         N_y .. y-line of N
					 2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
						\                              /                            */

      for (int i=0; i<iel; ++i) viscs2_(0,i) = 0.5 * (2.0 * derxy2_(0,i) + derxy2_(1,i));
      for (int i=0; i<iel; ++i) viscs2_(1,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(2,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(3,i) = 0.5 * (derxy2_(0,i) + 2.0 * derxy2_(1,i));
    }
    else
    {
      viscs2_.Clear();
    }

    /* momentum and velocity divergence: */
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (physicaltype == INPAR::FLUID::loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);
    //if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);


    {
      ////////////////// GALERKIN ANTEIL ///////////////////////

      /* Konvektionsterm (u*grad(u),v) + Anteil Massenmatrix*/
      for(int ui=0; ui<iel; ++ui)
      {
        const int tui = ndofs[ui];
        const double v = fac*densfunct_(ui)+timefacfac*conv_c_(ui);
        for(int vi=0; vi<iel; ++vi)
        {
          const int tvi = ndofs[vi];
          double v2 = v*funct_(vi) ;
          elestiff(tvi,     tui    ) += v2;
          elestiff(tvi + 1, tui + 1) += v2;
        }
      }


      /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
      const double viscefftimefacfac = visceff*timefacfac;
      for (int ui=0; ui<iel; ++ui)
      {
        const int tui  = ndofs[ui];
        const int tuip = tui+1;
        for (int vi=0; vi<iel; ++vi)
        {
          const int tvi  = ndofs[vi];
          const int tvip = tvi+1;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);

          elestiff(tvi,  tui ) += viscefftimefacfac*(2.0*derxy_0ui_0vi
              +
              derxy_1ui_1vi) ;
          elestiff(tvi,  tuip) += viscefftimefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          elestiff(tvip, tui ) += viscefftimefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          elestiff(tvip, tuip) += viscefftimefacfac*(derxy_0ui_0vi
              +
              2.0*derxy_1ui_1vi) ;

        }
      }

      // reaktive Konvektionsterme
      if (newton)
      {
        for (int vi=0; vi<iel; ++vi)
        {
          const int tvi  = ndofs[vi];
          const int tvip = tvi+1;
          const double v = timefacfac*funct_(vi);
          for (int ui=0; ui<iel; ++ui)
          {
            const int tui  = ndofs[ui];
            const int tuip = tui+1;
            const double v2 = v*densfunct_(ui);

            /*  convection, reactive part

						   /                           \
						   |  /          \   n+1       |
						   | | Du o nabla | u     , v  |
						   |  \          /   (i)       |
						   \                           /
             */
            elestiff(tvi,  tui ) += v2*vderxy_(0, 0) ;
            elestiff(tvi,  tuip) += v2*vderxy_(0, 1) ;
            elestiff(tvip, tui ) += v2*vderxy_(1, 0) ;
            elestiff(tvip, tuip) += v2*vderxy_(1, 1) ;
          }
        }
      }

      /* inertia */
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* inertia */
        const double v = -fac*funct_(vi);
        eleforce(tvi    ) += v*velint_(0) ;
        eleforce(tvi + 1) += v*velint_(1) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* convection */
        double v = -timefacfac*funct_(vi);
        eleforce(tvi    ) += v*(convvelint_(0)*vderxy_(0, 0)
            +
            convvelint_(1)*vderxy_(0, 1)) ;
        eleforce(tvi + 1) += v*(convvelint_(0)*vderxy_(1, 0)
            +
            convvelint_(1)*vderxy_(1, 1)) ;
      }

      // viscosity
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* viscosity */
        eleforce(tvi    ) -= visceff*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
            +
            derxy_(1, vi)*vderxy_(0, 1)
            +
            derxy_(1, vi)*vderxy_(1, 0)) ;
        eleforce(tvi + 1) -= visceff*timefacfac*(derxy_(0, vi)*vderxy_(0, 1)
            +
            derxy_(0, vi)*vderxy_(1, 0)
            +
            2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      // source term
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eleforce(tvi    ) += v*rhsmom_(0) ;
        eleforce(tvi + 1) += v*rhsmom_(1) ;
      }
    } // <- end of galerkin part
  } // <- end loop over all gausss points
  return 0;
}

/*!
* basic element call for implicit standard fluid solver
* inf-sup stable -> no stabilization terms
*
* @param ele Element
* @param params parameter list
* @param discretization
* @param lm vector with ids
* @param elemat1_epetra stiffmat (o)
* @param elevec1_epetra rhs (o)
* @return 0
*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid2TH<distype>::CalcSysmatAndResidual(Fluid2* ele,
    ParameterList& params,
    DRT::Discretization& discretization,
    vector<int> lm,
    Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseVector& elevec1_epetra,
    RefCountPtr< MAT::Material >   material/*,
	    							_MATERIAL *   material*/)
	    							{
  LINALG::Matrix<idofs,idofs> elemat(elemat1_epetra.A(),true);
  LINALG::Matrix<idofs, 1> eleres(elevec1_epetra.A(),true);


  int* ndofs = NULL;
  switch(distype)
  {
    case DRT::Element::quad9:
    case DRT::Element::quad8:
    case DRT::Element::quad4:
    {
      static int ndofarrayquad[9] = {0,3,6,9,12,14,16,18,20};
      ndofs = &ndofarrayquad[0];
    }
    break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      static int ndofarraytri[6] = {0,3,6,9,11,13};
      ndofs = &ndofarraytri[0];
    }
    break;
    default:
      dserror("distype not supported");
  }

  //RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  RCP<const Epetra_Vector> velnp = discretization.GetState("velaf");
  //RCP<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  RCP<const Epetra_Vector> vescnp = discretization.GetState("scaaf");
  RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  RCP<const Epetra_Vector> dispnp;
  RCP<const Epetra_Vector> gridv;

  vector<double> myvelnp(lm.size());
  vector<double> myvescnp(lm.size());
  vector<double> myhist(lm.size());
  vector<double> mydispnp;
  vector<double> mygridv;
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  DRT::UTILS::ExtractMyValues(*vescnp,myvescnp,lm);
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

  // ALE case
  if(ele->is_ale_)
  {
    //cout << "ALE must not be true!!" << endl;
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    gridv = discretization.GetState("gridv");
    if (gridv==null) dserror("Cannot get state vectors 'gridv'");
    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
  }

  LINALG::Matrix<iel,1>   eprenp;
  LINALG::Matrix<2,iel>	evelnp;
  LINALG::Matrix<iel,1>	edensnp;
  LINALG::Matrix<2,iel>	emhist;
  LINALG::Matrix<iel,1>	echist;
  LINALG::Matrix<2,iel>	edispnp;
  LINALG::Matrix<2,iel>	egridv;

  eprenp.Clear();
  evelnp.Clear();
  edensnp.PutScalar(1.0);
  emhist.Clear();
  echist.Clear();
  edispnp.Clear();
  egridv.Clear();

  for(int i=0; i<iel; i++)
  {
    evelnp(0,i) = myvelnp[ndofs[i]];
    evelnp(1,i) = myvelnp[ndofs[i]+1];
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
      {
        if(i<4)
        {
          eprenp(i) 	  = myvelnp[ndofs[i]+2];
          edensnp(i)  = myvescnp[ndofs[i]+2];

          echist(i) = myhist[ndofs[i]+2];
        }
      }
      break;
      case DRT::Element::tri6:
      case DRT::Element::tri3:
      {
        if(i<3)
        {
          eprenp(i)	= myvelnp[ndofs[i]+2];
          edensnp(i)  = myvescnp[ndofs[i]+2];
        }
      }
      break;
      default:
        dserror("distype not supported");
    }

    emhist(0,i) = myhist[ndofs[i]];
    emhist(1,i) = myhist[ndofs[i]+1];

    // ALE part
    if(ele->is_ale_)
    {
      edispnp(0,i) = mydispnp[ndofs[i]];
      edispnp(1,i) = mydispnp[ndofs[i]+1];
      egridv(0,i) = mygridv[ndofs[i]];
      egridv(1,i) = mygridv[ndofs[i]+1];
    }
  }

  // no fine-scale viscosity
  LINALG::Matrix<2,iel> fsevelnp;
  for(int i=0;i<iel;i++)
  {
    fsevelnp(0,i) = 0.0;
    fsevelnp(1,i) = 0.0;
  }
  Fluid2::FineSubgridVisc fssgv = Fluid2::no_fssgv;

  const double time = params.get<double>("total time",-1.0);

  bool newton = false;
  // bool loma = false;
  string newtonstr = params.get<string>("Linearisation");
  if(newtonstr=="Newton") newton = true;
  //string lomastr = params.get<string>("low-Mach-number solver");  // isn't supported anymore
  //if(lomastr=="Yes") loma = true;
  // loma = false;
  INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  ParameterList& stablist = params.sublist("STABILIZATION");

  Fluid2::TauType whichtau = Fluid2::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");
    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid2::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid2::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid2::codina;
    }
  }

  bool higher_order_ele = ele->isHigherOrderElement(distype);
  if(stablist.get<string>("STABTYPE") == "inconsistent") higher_order_ele = false;

  // time step
  const double dt = params.get<double>("dt");

  // One-step-Theta: timefac = theta*dt
  // BDF2:           timefac = 2/3 * dt
  const double theta = params.get<double>("theta",-1.0);
  const double timefac = dt * theta;
  if (timefac < 0.0) dserror("Negative time-integration parameter or time-step length supplied");

  // initialise the Smagorinsky constant Cs to zero
  double Cs            = 0.0;
  double visceff       = 0.0;

  // the default action is no model
  Fluid2::TurbModelAction turb_mod_action = Fluid2::no_model;

  // save node coordinates
  DRT::Node** nodes = ele->Nodes();
  for(int inode = 0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // ALE case: add displacements
  if(ele->is_ale_) xyze_+=edispnp;

  // dead load in element nodes
  BodyForce(ele,time);

  // check here, if we really have a fluid !!
  if( material->MaterialType() != INPAR::MAT::m_fluid
      && material->MaterialType() != INPAR::MAT::m_carreauyasuda
      && material->MaterialType() != INPAR::MAT::m_modpowerlaw) dserror("Material law is not a fluid");

  // viscosity
  double visc = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
    visc = actmat->Viscosity();
  }
  Caltau(ele,evelnp,fsevelnp,edensnp,whichtau,material,visc,timefac,dt,turb_mod_action,Cs,visceff,fssgv);

  // in case of viscous stabilization decide whether to use GLS or USFEM
  /*double vstabfac= 0.0;
	if (vstab == Fluid2::viscous_stab_usfem || vstab == Fluid2::viscous_stab_usfem_only_rhs)
	{
		vstabfac =  1.0;
	}
	else if(vstab == Fluid2::viscous_stab_gls || vstab == Fluid2::viscous_stab_gls_only_rhs)
	{
		vstabfac = -1.0;
	}*/

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
#ifdef DEBUG
  if(ele->gaussrule_!=DRT::UTILS::intrule_quad_9point)
    cout << "WARNING: no 9point gaussrule for quadrature! element matrices may be nontrivial disturbed!" << endl;
#endif

  // integration loop over all gausss points
  for(int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    pfunct_.PutScalar(0.0);
    pderiv_.PutScalar(0.0);

    // shape functions and their derivatives
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
        DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
        DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);
        DRT::UTILS::shape_function_2D(pfunct_,e1,e2,DRT::Element::quad4);
        DRT::UTILS::shape_function_2D_deriv1(pderiv_,e1,e2,DRT::Element::quad4);
        break;
      default: dserror("not yet implemented");
    }

    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);
    if(det<0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(),det);

    const double fac = intpoints.qwgt[iquad]*det;

    derxy_.Multiply(xji_,deriv_);
    pderxy_.Multiply(xji_,pderiv_);

    vderxy_.MultiplyNT(evelnp,derxy_);
    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    // (inverse-)density-weighted shape functions and global derivatives
    for (int inode=0; inode<iel; inode++)
    {
      densfunct_(inode) = edensnp(inode)*funct_(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    if (higher_order_ele)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,e1,e2,distype);
      // initialize everything
      LINALG::Matrix<3,3> bm;

      // calculate elements of jacobian_bar matrix
      bm(0,0) =                     xjm_(0,0)*xjm_(0,0);
      bm(0,1) =                     xjm_(0,1)*xjm_(0,1);
      bm(0,2) =                 2.0*xjm_(0,0)*xjm_(0,1);

      bm(1,0) =                     xjm_(1,0)*xjm_(1,0);
      bm(1,1) =                     xjm_(1,1)*xjm_(1,1);
      bm(1,2) =                 2.0*xjm_(1,1)*xjm_(1,0);

      bm(2,0) =                     xjm_(0,0)*xjm_(1,0);
      bm(2,1) =                     xjm_(0,1)*xjm_(1,1);
      bm(2,2) = xjm_(0,0)*xjm_(1,1)+xjm_(0,1)*xjm_(1,0);

      xder2_.MultiplyNT(deriv2_, xyze_);
      derxy2_.Multiply(-1.0, xder2_, derxy_);
      derxy2_ += deriv2_;

      LINALG::FixedSizeSerialDenseSolver<3,3,iel> solver_bm;
      solver_bm.SetMatrix(bm);                // A = bm
      // X == B is not a problem, derxy2_ will contain the solution.
      solver_bm.SetVectors(derxy2_,derxy2_);  // X = B = derxy2_
      // error code
      int ierr = 0;

      // Factor. Calling this directly is not necessary, only to find
      // out where an error came from.
      ierr = solver_bm.Factor();
      if (ierr!=0)
        dserror("Unable to perform LU factorisation during computation of derxy2");

      solver_bm.Solve();                      // Solve A*X = B
      if (ierr!=0)
        dserror("Unable to perform backward substitution after factorisation of jacobian");

      // calculate 2nd velocity derivatives at integration point
      vderxy2_.MultiplyNT(evelnp, derxy2_);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get momentum (n+g,i) at integration point
    velint_.Multiply(evelnp, densfunct_);

    // Moementenableitung zum Zeitpunkt n+1 am Integrationspunkt
    mderxy_.MultiplyNT(evelnp,densderxy_);

    // history data (n,i) am Integrationspunkt
    histmom_.Multiply(emhist, densfunct_);
    histcon_ = funct_.Dot(echist);

    // get convective velocity at integration point
    // We handle the ale case very implicitely here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_ = velint_;
    if (ele->is_ale_) convvelint_.Multiply(-1.0, egridv, densfunct_, 1.0);

    // get pressure gradient at integration point
    gradp_.Multiply(pderxy_, eprenp);	// ACHTUNG pderxy_!

    // get pressure at integration point
    double press = pfunct_.Dot(eprenp);

    //--------------------------------------------------------------
    // perform integration for entire matrix and rhs
    //--------------------------------------------------------------

    const double timefacfac = timefac * fac;

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_.Multiply(edeadng_, densfunct_);
    // get density at integration point
    double dens = funct_.Dot(edensnp);

    /*------------------------ evaluate rhs vectors at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    rhsmom_(0) = histmom_(0) + bodyforce_(0)*timefac;
    rhsmom_(1) = histmom_(1) + bodyforce_(1)*timefac;
    rhscon_ = histcon_ - dens;

    /*----------------- get numerical representation of single operators ---*/
    /* Convective term  u_old * grad u_old: */
    conv_old_.Multiply(vderxy_, convvelint_);

    /* Viscous term  div epsilon(u_old) */
    if (higher_order_ele)
    {
      visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,2));
      visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,2));
    }
    else
    {
      visc_old_.Clear();
    }

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
		          with  N .. form function matrix                                   */
    conv_c_.MultiplyTN(derxy_, convvelint_);

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
		              1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
		              - |                              |         N_y .. y-line of N
		             2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
		                \                              /                            */

      for (int i=0; i<iel; ++i) viscs2_(0,i) = 0.5 * (2.0 * derxy2_(0,i) + derxy2_(1,i));
      for (int i=0; i<iel; ++i) viscs2_(1,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(2,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(3,i) = 0.5 * (derxy2_(0,i) + 2.0 * derxy2_(1,i));
    }
    else
    {
      viscs2_.Clear();
    }

    // momentum and velocity divergence:
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (physicaltype == INPAR::FLUID::loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);
    // if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    /////////////////////////////////
    // build fluid matrix
    {
      ////////////////// GALERKIN ANTEIL ///////////////////////

      /* Konvektionsterm (u*grad(u),v) + Anteil Massenmatrix*/
      for(int ui=0; ui<iel; ++ui)
      {
        const int tui = ndofs[ui];
        const double v = fac*densfunct_(ui)+timefacfac*conv_c_(ui);
        for(int vi=0; vi<iel; ++vi)
        {
          const int tvi = ndofs[vi];
          double v2 = v*funct_(vi) ;
          elemat(tvi,     tui    ) += v2;
          elemat(tvi + 1, tui + 1) += v2;
        }
      }


      /* Viskositätsterm (2*nu*epsilon(u),epsilon(v)) */
      const double viscefftimefacfac = visceff*timefacfac;
      for (int ui=0; ui<iel; ++ui)
      {
        const int tui  = ndofs[ui];
        const int tuip = tui+1;
        for (int vi=0; vi<iel; ++vi)
        {
          const int tvi  = ndofs[vi];
          const int tvip = tvi+1;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);

          elemat(tvi,  tui ) += viscefftimefacfac*(2.0*derxy_0ui_0vi
              +
              derxy_1ui_1vi) ;
          elemat(tvi,  tuip) += viscefftimefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          elemat(tvip, tui ) += viscefftimefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          elemat(tvip, tuip) += viscefftimefacfac*(derxy_0ui_0vi
              +
              2.0*derxy_1ui_1vi) ;

        }
      }

      // Druckterm (gradp) entspricht Delta t theta (Dp, nabla o v)
      for(int i=0;i<iel;i++)
      {
        for(int p=0; p<4; p++) // only pressure dofs //TODO no tri elements supported!
        {
          elemat(ndofs[i],3*p+2)   -= timefacfac*pfunct_(p)*derxy_(0,i);	// 09.01 VZ verändert
          elemat(ndofs[i]+1,3*p+2) -= timefacfac*pfunct_(p)*derxy_(1,i);
        }
      }

      // Divergenzfreiheit
      for(int p=0; p<4; p++)	// Schleife über alle Druckzeilen
      {
        for(int i=0; i<iel; i++)	// Schleife über alle Geschwindigkeitsspalten
        {
          elemat(3*p+2,ndofs[i])   += timefacfac*pfunct_(p)*derxy_(0,i);
          elemat(3*p+2,ndofs[i]+1) += timefacfac*pfunct_(p)*derxy_(1,i);
        }
      }

      // reaktive Konvektionsterme
      if (newton)
      {
        for (int vi=0; vi<iel; ++vi)
        {
          const int tvi  = ndofs[vi];
          const int tvip = tvi+1;
          const double v = timefacfac*funct_(vi);
          for (int ui=0; ui<iel; ++ui)
          {
            const int tui  = ndofs[ui];
            const int tuip = tui+1;
            const double v2 = v*densfunct_(ui);

            /*  convection, reactive part

			               /                           \
			               |  /          \   n+1       |
			               | | Du o nabla | u     , v  |
			               |  \          /   (i)       |
			               \                           /
             */
            elemat(tvi,  tui ) += v2*vderxy_(0, 0) ;
            elemat(tvi,  tuip) += v2*vderxy_(0, 1) ;
            elemat(tvip, tui ) += v2*vderxy_(1, 0) ;
            elemat(tvip, tuip) += v2*vderxy_(1, 1) ;
          }
        }
      }
    }	// <- end build fluid matrix

    /////////////////////////////////
    // build rhs
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* inertia */
        const double v = -fac*funct_(vi);
        eleres(tvi    ) += v*velint_(0) ;
        eleres(tvi + 1) += v*velint_(1) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* convection */
        double v = -timefacfac*funct_(vi);
        eleres(tvi    ) += v*(convvelint_(0)*vderxy_(0, 0)
            +
            convvelint_(1)*vderxy_(0, 1)) ;
        eleres(tvi + 1) += v*(convvelint_(0)*vderxy_(1, 0)
            +
            convvelint_(1)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<iel; ++vi)	// this term only here and not in CalcIntUsfem
      {
        const int tvi = ndofs[vi];
        /* pressure */
        double v = press*timefacfac;
        eleres(tvi    ) += v*derxy_(0, vi) ;
        eleres(tvi + 1) += v*derxy_(1, vi) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* viscosity */
        eleres(tvi    ) -= visceff*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
            +
            derxy_(1, vi)*vderxy_(0, 1)
            +
            derxy_(1, vi)*vderxy_(1, 0)) ;
        eleres(tvi + 1) -= visceff*timefacfac*(derxy_(0, vi)*vderxy_(0, 1)
            +
            derxy_(0, vi)*vderxy_(1, 0)
            +
            2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eleres(tvi    ) += v*rhsmom_(0) ;
        eleres(tvi + 1) += v*rhsmom_(1) ;
      }

      {
        //const double timefacfac_mdiv = timefacfac * mdiv_;
        for (int vi=0; vi<iel; ++vi)
        {
          // continuity equation
          if(vi < 4)		// TODO no tri elements supported
          {
            eleres(ndofs[vi]+2) -= timefacfac*pfunct_(vi)*vdiv_;
          }
        }
      }
    } // end of rhs
  } // <- end loop over all gauss points
  return 0;
	    							}

/*!
* basic element call for
* calculation of rhs only
*
* @param ele Element
* @param params parameter list
* @param discretization
* @param lm vector with ids
* @param elevec1_epetra rhs (o)
* @return 0
*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Fluid2TH<distype>::CalcResidual(Fluid2* ele,
    ParameterList& params,
    DRT::Discretization& discretization,
    vector<int> lm,
    Epetra_SerialDenseVector& elevec_epetra,
    RefCountPtr< MAT::Material >   material/*,
						_MATERIAL *   material*/)
						{

  LINALG::Matrix<idofs, 1> eleres(elevec_epetra.A(),true);

  int* ndofs = NULL;
  switch(distype)
  {
    case DRT::Element::quad9:
    case DRT::Element::quad8:
    case DRT::Element::quad4:
    {
      static int ndofarrayquad[9] = {0,3,6,9,12,14,16,18,20};
      ndofs = &ndofarrayquad[0];
    }
    break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
    {
      static int ndofarraytri[6] = {0,3,6,9,11,13};
      ndofs = &ndofarraytri[0];
    }
    break;
    default:
      dserror("distype not supported");
  }

  RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
  //RCP<const Epetra_Vector> vedenp = discretization.GetState("vedenp");
  RCP<const Epetra_Vector> vedenp = discretization.GetState("scanp");
  RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  RCP<const Epetra_Vector> dispnp;
  RCP<const Epetra_Vector> gridv;

  vector<double> myvelnp(lm.size());
  vector<double> myvedenp(lm.size());
  vector<double> myhist(lm.size());
  vector<double> mydispnp;
  vector<double> mygridv;
  DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
  DRT::UTILS::ExtractMyValues(*vedenp,myvedenp,lm);
  DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

  // ALE case
  if(ele->is_ale_)
  {
    dispnp = discretization.GetState("dispnp");
    if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
    mydispnp.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    gridv = discretization.GetState("gridv");
    if (gridv==null) dserror("Cannot get state vectors 'gridv'");
    mygridv.resize(lm.size());
    DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
  }

  LINALG::Matrix<iel,1>   epren;
  LINALG::Matrix<2,iel>	evelnp;
  LINALG::Matrix<iel,1>	edensnp;
  LINALG::Matrix<2,iel>	emhist;
  LINALG::Matrix<iel,1>	echist;
  LINALG::Matrix<2,iel>	edispnp;
  LINALG::Matrix<2,iel>	egridv;

  epren.Clear();
  evelnp.Clear();
  edensnp.PutScalar(1.0);
  emhist.Clear();
  echist.Clear();
  edispnp.Clear();
  egridv.Clear();

  for(int i=0; i<iel; i++)
  {
    evelnp(0,i) = myvelnp[ndofs[i]];
    evelnp(1,i) = myvelnp[ndofs[i]+1];
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
      {
        if(i<4)
        {
          epren(i) 	  = myvelnp[ndofs[i]+2];
          edensnp(i)  = myvedenp[ndofs[i]+2];
        }
      }
      break;
      case DRT::Element::tri6:
      case DRT::Element::tri3:
      {
        if(i<3)
        {
          epren(i)	= myvelnp[ndofs[i]+2];
          edensnp(i)  = myvedenp[ndofs[i]+2];
        }
      }
      break;
      default:
        dserror("distype not supported");
    }

    emhist(0,i) = myhist[ndofs[i]];
    emhist(1,i) = myhist[ndofs[i]+1];

    // ALE case
    if(ele->is_ale_)
    {
      edispnp(0,i) = mydispnp[ndofs[i]];
      edispnp(1,i) = mydispnp[ndofs[i]+1];
      egridv(0,i) = mygridv[ndofs[i]];
      egridv(1,i) = mygridv[ndofs[i]+1];
    }
  }

  // no fine-scale viscosity
  LINALG::Matrix<2,iel> fsevelnp;
  for(int i=0;i<iel;i++)
  {
    fsevelnp(0,i) = 0.0;
    fsevelnp(1,i) = 0.0;
  }
  Fluid2::FineSubgridVisc fssgv = Fluid2::no_fssgv;

  const double time = params.get<double>("total time",-1.0);

  bool newton = false;
  // bool loma = false;
  string newtonstr = params.get<string>("Linearisation");
  // string lomastr = params.get<string>("low-Mach-number solver");
  if(newtonstr=="Newton") newton = true;
  // if(lomastr=="Yes") loma = true;
  INPAR::FLUID::PhysicalType physicaltype = params.get<INPAR::FLUID::PhysicalType>("Physical Type");

  ParameterList& stablist = params.sublist("STABILIZATION");

  //Fluid2::StabilisationAction cross    = ele->ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
  //Fluid2::StabilisationAction reynolds = ele->ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

  Fluid2::TauType whichtau = Fluid2::tau_not_defined;
  {
    const string taudef = stablist.get<string>("DEFINITION_TAU");
    if(taudef == "Barrenechea_Franca_Valentin_Wall")
    {
      whichtau = Fluid2::franca_barrenechea_valentin_wall;
    }
    else if(taudef == "Bazilevs")
    {
      whichtau = Fluid2::bazilevs;
    }
    else if(taudef == "Codina")
    {
      whichtau = Fluid2::codina;
    }
  }

  bool higher_order_ele = ele->isHigherOrderElement(distype);
  if(stablist.get<string>("STABTYPE") == "inconsistent") higher_order_ele = false;

  const double dt = params.get<double>("dt");

  // One-step-Theta: timefac = theta*dt
  // BDF2:           timefac = 2/3 * dt
  const double timefac = params.get<double>("thsl",-1.0);
  if (timefac < 0.0) dserror("No thsl supplied");

  // --------------------------------------------------
  // set parameters for classical turbulence models
  // --------------------------------------------------
  //ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");

  // initialise the Smagorinsky constant Cs to zero
  double Cs            = 0.0;
  double visceff       = 0.0;

  // the default action is no model
  Fluid2::TurbModelAction turb_mod_action = Fluid2::no_model;

  // save node coordinates
  DRT::Node** nodes = ele->Nodes();
  for(int inode = 0; inode<iel; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // ALE case
  if(ele->is_ale_) xyze_+=edispnp;

  // dead load in element nodes
  BodyForce(ele,time);

  // check here, if we really have a fluid !!
  if( material->MaterialType() != INPAR::MAT::m_fluid
      && material->MaterialType() != INPAR::MAT::m_carreauyasuda
      && material->MaterialType() != INPAR::MAT::m_modpowerlaw) dserror("Material law is not a fluid");

  // viscosity
  double visc = 0.0;
  if(material->MaterialType() == INPAR::MAT::m_fluid)
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());
    visc = actmat->Viscosity();
  }

  Caltau(ele,evelnp,fsevelnp,edensnp,whichtau,material,visc,timefac,dt,turb_mod_action,Cs,visceff,fssgv);


  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
#ifdef DEBUG
  if(ele->gaussrule_!=DRT::UTILS::intrule_quad_9point)
    cout << "WARNING: no 9point gaussrule for quadrature! element matrices may be nontrivial disturbed!" << endl;
#endif

  // integration loop over all gauss points
  for(int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];

    pfunct_.PutScalar(0.0);
    pderiv_.PutScalar(0.0);

    // shape functions and their derivative
    switch(distype)
    {
      case DRT::Element::quad9:
      case DRT::Element::quad8:
      case DRT::Element::quad4:
        DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
        DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);
        DRT::UTILS::shape_function_2D(pfunct_,e1,e2,DRT::Element::quad4);
        DRT::UTILS::shape_function_2D_deriv1(pderiv_,e1,e2,DRT::Element::quad4);
        break;
      default: dserror("not yet implemented");
    }

    xjm_.MultiplyNT(deriv_,xyze_);
    const double det = xji_.Invert(xjm_);
    if(det<0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f",ele->Id(),det);

    const double fac = intpoints.qwgt[iquad]*det;

    derxy_.Multiply(xji_,deriv_);
    pderxy_.Multiply(xji_,pderiv_);


    vderxy_.MultiplyNT(evelnp,derxy_);
    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    // (inverse-)density-weighted shape functions and global derivatives
    for (int inode=0; inode<iel; inode++)
    {
      densfunct_(inode) = edensnp(inode)*funct_(inode);

      densderxy_(0,inode) = edensnp(inode)*derxy_(0,inode);
      densderxy_(1,inode) = edensnp(inode)*derxy_(1,inode);
    }

    if (higher_order_ele)
    {
      // get values of shape functions and derivatives in the gausspoint
      DRT::UTILS::shape_function_2D_deriv2(deriv2_,e1,e2,distype);
      // initialize everything
      LINALG::Matrix<3,3> bm;

      // calculate elements of jacobian_bar matrix
      bm(0,0) =                     xjm_(0,0)*xjm_(0,0);
      bm(0,1) =                     xjm_(0,1)*xjm_(0,1);
      bm(0,2) =                 2.0*xjm_(0,0)*xjm_(0,1);

      bm(1,0) =                     xjm_(1,0)*xjm_(1,0);
      bm(1,1) =                     xjm_(1,1)*xjm_(1,1);
      bm(1,2) =                 2.0*xjm_(1,1)*xjm_(1,0);

      bm(2,0) =                     xjm_(0,0)*xjm_(1,0);
      bm(2,1) =                     xjm_(0,1)*xjm_(1,1);
      bm(2,2) = xjm_(0,0)*xjm_(1,1)+xjm_(0,1)*xjm_(1,0);

      xder2_.MultiplyNT(deriv2_, xyze_);
      derxy2_.Multiply(-1.0, xder2_, derxy_);
      derxy2_ += deriv2_;

      LINALG::FixedSizeSerialDenseSolver<3,3,iel> solver_bm;
      solver_bm.SetMatrix(bm);                // A = bm
      // X == B is not a problem, derxy2_ will contain the solution.
      solver_bm.SetVectors(derxy2_,derxy2_);  // X = B = derxy2_
      // error code
      int ierr = 0;

      // Factor. Calling this directly is not necessary, only to find
      // out where an error came from.
      ierr = solver_bm.Factor();
      if (ierr!=0)
        dserror("Unable to perform LU factorisation during computation of derxy2");

      solver_bm.Solve();                      // Solve A*X = B
      if (ierr!=0)
        dserror("Unable to perform backward substitution after factorisation of jacobian");

      // calculate 2nd velocity derivatives at integration point
      vderxy2_.MultiplyNT(evelnp, derxy2_);
    }
    else
    {
      derxy2_  = 0.;
      vderxy2_ = 0.;
    }

    // get momentum (n+g,i) at integration point
    velint_.Multiply(evelnp, densfunct_);

    // history data (n,i) am Integrationspunkt
    histmom_.Multiply(emhist, densfunct_);
    histcon_ = funct_.Dot(echist);

    // get convective velocity at integration point
    // We handle the ale case very implicitely here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_ = velint_;
    if (ele->is_ale_) convvelint_.Multiply(-1.0, egridv, densfunct_, 1.0);

    // get pressure gradient at integration point
    gradp_.Multiply(pderxy_, epren);

    // get pressure at integration point
    double press = pfunct_.Dot(epren);

    //--------------------------------------------------------------
    // perform integration for entire matrix and rhs
    //--------------------------------------------------------------

    const double timefacfac = timefac * fac;

    // get (density-weighted) bodyforce in gausspoint
    bodyforce_.Multiply(edeadng_, densfunct_);
    // get density at integration point
    double dens = funct_.Dot(edensnp);

    /*------------------------ evaluate rhs vectors at integration point ---*/
    // no switch here at the moment w.r.t. is_ale
    rhsmom_(0) = histmom_(0) + bodyforce_(0)*timefac;
    rhsmom_(1) = histmom_(1) + bodyforce_(1)*timefac;
    rhscon_ = histcon_ - dens;

    /*----------------- get numerical representation of single operators ---*/
    /* Convective term  u_old * grad u_old: */
    conv_old_.Multiply(vderxy_, convvelint_);

    /* Viscous term  div epsilon(u_old) */
    if (higher_order_ele)
    {
      visc_old_(0) = vderxy2_(0,0) + 0.5 * (vderxy2_(0,1) + vderxy2_(1,2));
      visc_old_(1) = vderxy2_(1,1) + 0.5 * (vderxy2_(1,0) + vderxy2_(0,2));
    }
    else
    {
      visc_old_.Clear();
    }

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
		          with  N .. form function matrix                                   */
    conv_c_.MultiplyTN(derxy_, convvelint_);

    if (higher_order_ele)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
		              1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
		              - |                              |         N_y .. y-line of N
    	             2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
		                \                              /                            */

      for (int i=0; i<iel; ++i) viscs2_(0,i) = 0.5 * (2.0 * derxy2_(0,i) + derxy2_(1,i));
      for (int i=0; i<iel; ++i) viscs2_(1,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(2,i) = 0.5 *  derxy2_(2,i);
      for (int i=0; i<iel; ++i) viscs2_(3,i) = 0.5 * (derxy2_(0,i) + 2.0 * derxy2_(1,i));
    }
    else
    {
      viscs2_.Clear();
    }

    /* momentum and velocity divergence: */
    mdiv_ = mderxy_(0, 0) + mderxy_(1, 1);
    if (physicaltype == INPAR::FLUID::loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);
    // if (loma) vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    /////////////////////////////////
    // get rhs
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* inertia */
        const double v = -fac*funct_(vi);
        eleres(tvi    ) += v*velint_(0) ;
        eleres(tvi + 1) += v*velint_(1) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* convection */
        double v = -timefacfac*funct_(vi);
        eleres(tvi    ) += v*(convvelint_(0)*vderxy_(0, 0)
            +
            convvelint_(1)*vderxy_(0, 1)) ;
        eleres(tvi + 1) += v*(convvelint_(0)*vderxy_(1, 0)
            +
            convvelint_(1)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<iel; ++vi)	// this term only here and not in CalcIntUsfem
      {
        const int tvi = ndofs[vi];
        /* pressure */
        double v = press*timefacfac;
        eleres(tvi    ) += v*derxy_(0, vi) ;
        eleres(tvi + 1) += v*derxy_(1, vi) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        /* viscosity */
        eleres(tvi    ) -= visceff*timefacfac*(2.0*derxy_(0, vi)*vderxy_(0, 0)
            +
            derxy_(1, vi)*vderxy_(0, 1)
            +
            derxy_(1, vi)*vderxy_(1, 0)) ;
        eleres(tvi + 1) -= visceff*timefacfac*(derxy_(0, vi)*vderxy_(0, 1)
            +
            derxy_(0, vi)*vderxy_(1, 0)
            +
            2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      for (int vi=0; vi<iel; ++vi)
      {
        const int tvi = ndofs[vi];
        // source term of the right hand side
        const double v = fac*funct_(vi);
        eleres(tvi    ) += v*rhsmom_(0) ;
        eleres(tvi + 1) += v*rhsmom_(1) ;
      }

      {		// this term only here and not in CalcIntUsfem
        //const double timefacfac_mdiv = timefacfac * mdiv_;
        for (int vi=0; vi<iel; ++vi)
        {
          // continuity equation
          if(vi < 4)		// TODO tri elements not suported
          {
            eleres(ndofs[vi]+2) -= timefacfac*pfunct_(vi)*vdiv_;
          }
        }
      }
    } // end Galerkin part
  } // <- end loop over all gauss points
  return 0;
						}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2TH<distype>::Caltau(
    Fluid2*                                           ele,
    const LINALG::Matrix<2,iel>&  evelnp,
    const LINALG::Matrix<2,iel>&  fsevelnp,
    const LINALG::Matrix<iel,1>&  edensnp,
    const enum Fluid2::TauType                        whichtau,
    Teuchos::RCP<const MAT::Material>                 material,
    double&                                           visc,
    const double                                      timefac,
    const double                                      dt,
    const enum Fluid2::TurbModelAction                turb_mod_action,
    double&                                           Cs,
    double&                                           visceff,
    const enum Fluid2::FineSubgridVisc      		  fssgv
)
{

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule2D integrationrule_stabili=DRT::UTILS::intrule2D_undefined;
  switch (distype)
  {
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      integrationrule_stabili = DRT::UTILS::intrule_quad_1point;
      break;
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      integrationrule_stabili = DRT::UTILS::intrule_tri_1point;
      break;
    default:
      dserror("invalid discretization type for fluid2");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double wquad = intpoints.qwgt[0];

  DRT::UTILS::shape_function_2D(funct_,e1,e2,distype);
  DRT::UTILS::shape_function_2D_deriv1(deriv_,e1,e2,distype);

  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
    case DRT::Element::tri3:
    case DRT::Element::quad4:
      mk = 0.333333333333333333333;
      break;
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::tri6:
      mk = 0.083333333333333333333;
      break;
    default:
      dserror("type unknown!\n");
  }

  // get velocities at element center
  velint_.Multiply(evelnp, funct_);

  // get density at element center
  const double dens = funct_.Dot(edensnp);

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_, xyze_);
  // inverse of jacobian
  const double det = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det < 0.0) dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);

  // compute element area
  const double area = wquad*det;

  // get characteristic element length: square root of element area
  const double hk = sqrt(area);

  //             compute global first derivates
  //
  // this is necessary only for the calculation of the
  // streamlength (required by the quasistatic formulation)
  //
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

   */

  // compute global derivates
  derxy_.Multiply(xji_, deriv_);

  // get velocity (np,i) derivatives at integration point
  vderxy_.MultiplyNT(evelnp, derxy_);

  // get velocity norm
  const double vel_norm = velint_.Norm2();

  // normed velocity at element centre (currently not used)
  //if (vel_norm>=1e-6)
  //{
  //  velino_ = velint_/vel_norm;
  //}
  //else
  //{
  //  velino_ = 0.;
  //  velino_(0) = 1;
  //}

  // get streamlength (currently not used)
  //const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
  //const double strle = 2.0/val;

  /*------------------------------------------------------------------*/
  /*                                                                  */
  /*                 GET EFFECTIVE VISCOSITY IN GAUSSPOINT            */
  /*                                                                  */
  /* This part is used to specify an effective viscosity. This eff.   */
  /* viscosity may be caused by a Smagorinsky model                   */
  /*                                                                  */
  /*          visc    = visc + visc                                   */
  /*              eff              turbulent                          */
  /*                                                                  */
  /* here, the latter turbulent viscosity is not a material thing,    */
  /* but a flow feature!                                              */
  /*                                                                  */
  /* Another cause for the necessity of an effective viscosity might  */
  /* be the use of a shear thinning Non-Newtonian fluid               */
  /*                                                                  */
  /*                            /         \                           */
  /*            visc    = visc | shearrate |                          */
  /*                eff         \         /                           */
  /*                                                                  */
  /*                                                                  */
  /* Mind that at the moment all stabilization (tau and viscous test  */
  /* functions if applied) are based on the material viscosity not    */
  /* the effective viscosity. We do this since we do not evaluate the */
  /* stabilisation parameter in the gausspoints but just once in the  */
  /* middle of the element.                                           */
  /*------------------------------------------------------------------*/

  // compute nonlinear viscosity according to the Carreau-Yasuda model
  //if( material->mattyp != m_fluid ) CalVisc( material, visc);
  if ( material->MaterialType() != INPAR::MAT::m_fluid )
    CalVisc(material, visc);


  if (turb_mod_action == Fluid2::smagorinsky)
  {
    //
    // SMAGORINSKY MODEL
    // -----------------
    //                                   +-                                 -+ 1
    //                               2   |          / h \           / h \    | -
    //    visc          = dens * lmix  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent           |      |          \   / ij        \   / ij |
    //                            |      +-                                 -+
    //                            |
    //                            |      |                                   |
    //                            |      +-----------------------------------+
    //                            |           'resolved' rate of strain
    //                         mixing length
    //

    double rateofstrain = 0;
    {
      double epsilon;
      for(int rr=0;rr<2;rr++)
      {
        for(int mm=0;mm<rr;mm++)
        {
          epsilon = vderxy_(rr,mm) + vderxy_(mm,rr);
          rateofstrain += epsilon*epsilon;
        }
        rateofstrain += 2.0 * vderxy_(rr,rr) * vderxy_(rr,rr);
      }
      rateofstrain = sqrt(rateofstrain);
    }
    //
    // Choices of the Smagorinsky constant Cs:
    //
    //             Cs = 0.17   (Lilly --- Determined from filter
    //                          analysis of Kolmogorov spectrum of
    //                          isotropic turbulence)
    //
    //             0.1 < Cs < 0.24 (depending on the flow)
    //
    // mixing length set proportional to grid witdh
    //
    //                     lmix = Cs * hk

    double lmix = Cs * hk;

    //
    //          visc    = visc + visc
    //              eff              turbulent

    visceff = visc + dens * lmix * lmix * rateofstrain;
  }
  else
  {
    visceff = visc;
  }

  // calculate tau

  if (whichtau == Fluid2::franca_barrenechea_valentin_wall)
  {
    /*----------------------------------------------------- compute tau_Mu ---*/
    /* stability parameter definition according to

	    Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
	    element method for a generalized Stokes problem. Numerische
	    Mathematik, Vol. 92, pp. 652-677, 2002.
	    http://www.lncc.br/~valentin/publication.htm

	    and:

	    Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
	    Finite Element Method for the Advective-Reactive-Diffusive
	    Equation. Computer Methods in Applied Mechanics and Enginnering,
	    Vol. 190, pp. 1785-1800, 2000.
	    http://www.lncc.br/~valentin/publication.htm                   */


    /* viscous : reactive forces */
    const double re1 = 4.0 * timefac * visceff / (mk * dens * DSQR(hk));
    //cout << "mk: " << mk << " dens: " << dens << " hk: " << hk << endl;

    /* convective : viscous forces */
    const double re2 = mk * dens * vel_norm * hk / (2.0 * visceff);
    //cout << " re2 " << re2 << " visceff " << visceff << endl;

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = DSQR(hk)/(DSQR(hk)*dens*xi1+(4.0*timefac*visceff/mk)*xi2);
    tau_(1) = tau_(0);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    //tau[2] = sqrt(DSQR(visc)+DSQR(0.5*vel_norm*hk));

    // Wall Diss. 99
    /*
	                      xi2 ^
	                          |
	                        1 |   +-----------
	                          |  /
	                          | /
	                          |/
	                          +--------------> Re2
	                              1
     */
    const double xi_tau_c = DMIN(re2,1.0);
    tau_(2) = vel_norm * hk * 0.5 * xi_tau_c / ( timefac * dens );
  }
  else if(whichtau == Fluid2::bazilevs)
  {
    /* INSTATIONARY FLOW PROBLEM, ONE-STEP-THETA, BDF2

	    tau_M: Bazilevs et al.
	                                                               1.0
	                 +-                                       -+ - ---
	                 |                                         |   2.0
	                 | 4.0    n+1       n+1          2         |
	          tau  = | --- + u     * G u     + C * nu  * G : G |
	             M   |   2           -          I        -   - |
	                 | dt            -                   -   - |
	                 +-                                       -+

	   tau_C: Bazilevs et al., derived from the fine scale complement Shur
	          operator of the pressure equation


	                                  1.0
	                    tau  = -----------------
	                       C            /     \
	                            tau  * | g * g |
	                               M    \-   -/
     */

    /*            +-           -+   +-           -+   +-           -+
	                  |             |   |             |   |             |
	                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
	            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
	             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
	                  |    i     j  |   |    i     j  |   |    i     j  |
	                  +-           -+   +-           -+   +-           -+
     */
    double g;

    /*            +----
	                   \
	          G : G =   +   G   * G
	          -   -    /     ij    ij
	          -   -   +----
	                   i,j
     */
    double normG = 0;

    /*                      +----
	           n+1       n+1     \     n+1          n+1
	          u     * G u     =   +   u    * G   * u
	                  -          /     i     -ij    j
	                  -         +----        -
	                             i,j
     */
    double Gnormu = 0;

    const double dens_sqr = dens*dens;
    for (int nn=0;nn<2;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*velint_(nn);
      for (int rr=0;rr<2;++rr)
      {
        g = xji_(nn,0)*xji_(rr,0);
        for (int mm=1;mm<3;++mm)
        {
          g += xji_(nn,mm)*xji_(rr,mm);
        }
        normG += g*g;
        Gnormu+=dens_sqr_velint_nn*g*velint_(rr);
      }
    }

    // definition of constant
    // (Akkerman et al. (2008) used 36.0 for quadratics, but Stefan
    //  brought 144.0 from Austin...)
    const double CI = 12.0/mk;

    /*                                                         1.0
	                 +-                                       -+ - ---
	                 |                                         |   2.0
	                 | 4.0    n+1       n+1          2         |
	          tau  = | --- + u     * G u     + C * nu  * G : G |
	             M   |   2           -          I        -   - |
	                 | dt            -                   -   - |
	                 +-                                       -+
     */
    tau_(0) = 1.0/(timefac*sqrt((4.0*dens_sqr)/(dt*dt)+Gnormu+CI*visceff*visceff*normG));
    tau_(1) = tau_(0);

    /*           +-     -+   +-     -+   +-     -+
	                 |       |   |       |   |       |
	                 |  dr   |   |  ds   |   |  dt   |
	            g  = |  ---  | + |  ---  | + |  ---  |
	             i   |  dx   |   |  dx   |   |  dx   |
	                 |    i  |   |    i  |   |    i  |
	                 +-     -+   +-     -+   +-     -+
     */
    const double g0 = xji_(0,0) + xji_(0,1);
    const double g1 = xji_(1,0) + xji_(1,1);

    /*           +----
	                  \
	         g * g =   +   g * g
	         -   -    /     i   i
	                 +----
	                   i
     */
    const double normgsq = g0*g0+g1*g1;

    /*
	                                1.0
	                  tau  = -----------------
	                     C            /     \
	                          tau  * | g * g |
	                             M    \-   -/
     */
    tau_(2) = 1./(tau_(0)*normgsq*timefac*timefac*dens_sqr);

  }
  else if(whichtau == Fluid2::codina)
  {
    /*----------------------------------------------------- compute tau_Mu ---*/
    /* stability parameter definition according to

	    Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
	    element method for a generalized Stokes problem. Numerische
	    Mathematik, Vol. 92, pp. 652-677, 2002.
	    http://www.lncc.br/~valentin/publication.htm

	    and:

	    Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
	    Finite Element Method for the Advective-Reactive-Diffusive
	    Equation. Computer Methods in Applied Mechanics and Enginnering,
	    Vol. 190, pp. 1785-1800, 2000.
	    http://www.lncc.br/~valentin/publication.htm                   */


    /* viscous : reactive forces */
    const double re1 = 4.0 * timefac * visceff / (mk * dens * DSQR(hk));

    /* convective : viscous forces */
    const double re2 = mk * dens * vel_norm * hk / (2.0 * visceff);

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = DSQR(hk)/(DSQR(hk)*dens*xi1+(4.0*timefac*visceff/mk)*xi2);
    tau_(1) = tau_(0);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    tau_(2) = sqrt(DSQR(visceff)+DSQR(0.5*dens*vel_norm*hk)) / ( timefac*dens*dens );

  }
  else
  {
    dserror("unknown definition of tau\n");
  }

  /*------------------------------------------- compute subgrid viscosity ---*/
  if (fssgv != Fluid2::no_fssgv)
  {
    if (fssgv == Fluid2::smagorinsky_all or
        fssgv == Fluid2::smagorinsky_small)
    {
      //
      // SMAGORINSKY MODEL
      // -----------------
      //                                      +-                                 -+ 1
      //                                  2   |          / h \           / h \    | -
      //    visc          = dens * (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
      //        turbulent                     |          \   / ij        \   / ij |
      //                                      +-                                 -+
      //                                      |                                   |
      //                                      +-----------------------------------+
      //                                            'resolved' rate of strain
      //

      double rateofstrain = 0.0;
      {
        // get fine-scale or all-scale velocity (np,i) derivatives at element center
        if (fssgv == Fluid2::smagorinsky_small)
          fsvderxy_.MultiplyNT(fsevelnp, derxy_);
        else fsvderxy_.MultiplyNT(evelnp, derxy_);

        double epsilon;
        for(int rr=0;rr<2;rr++)
        {
          for(int mm=0;mm<rr;mm++)
          {
            epsilon = fsvderxy_(rr,mm) + fsvderxy_(mm,rr);
            rateofstrain += epsilon*epsilon;
          }
          rateofstrain += 2.0 * fsvderxy_(rr,rr) * fsvderxy_(rr,rr);
        }
        rateofstrain = sqrt(rateofstrain);
      }
      //
      // Choices of the fine-scale Smagorinsky constant Cs:
      //
      //             Cs = 0.17   (Lilly --- Determined from filter
      //                          analysis of Kolmogorov spectrum of
      //                          isotropic turbulence)
      //
      //             0.1 < Cs < 0.24 (depending on the flow)

      vart_ = dens * Cs * Cs * hk * hk * rateofstrain;
    }
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2TH<distype>::CalVisc(
    Teuchos::RCP<const MAT::Material> material,
    double&                 visc)
    {
  // compute shear rate
  double rateofshear = 0.0;

  double epsilon;
  for(int rr=0;rr<2;rr++) {
    for(int mm=0;mm<rr;mm++) {
      epsilon = vderxy_(rr,mm) + vderxy_(mm,rr);
      rateofshear += epsilon*epsilon;
    }
    rateofshear += 2.0 * vderxy_(rr,rr) * vderxy_(rr,rr);
  }

  rateofshear = sqrt(rateofshear);

  if(material->MaterialType() == INPAR::MAT::m_carreauyasuda)
  {
    const MAT::CarreauYasuda* actmat = static_cast<const MAT::CarreauYasuda*>(material.get());

    double nu_0   = actmat->Nu0();    // parameter for zero-shear viscosity
    double nu_inf = actmat->NuInf();  // parameter for infinite-shear viscosity
    double lambda = actmat->Lambda(); // parameter for characteristic time
    double a    = actmat->AParam(); // constant parameter
    double b      = actmat->BParam(); // constant parameter

    // compute viscosity according to the Carreau-Yasuda model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    const double tmp = pow(lambda*rateofshear,b);
    visc = nu_inf + ((nu_0 - nu_inf)/pow((1 + tmp),a));
  }
  else if(material->MaterialType() == INPAR::MAT::m_modpowerlaw)
  {
    const MAT::ModPowerLaw* actmat = static_cast<const MAT::ModPowerLaw*>(material.get());

    // get material parameters
    double m     = actmat->MCons();     // consistency constant
    double delta = actmat->Delta();      // safety factor
    double a     = actmat->AExp();      // exponent

    // compute viscosity according to a modified power law model for shear-thinning fluids
    // see Dhruv Arora, Computational Hemodynamics: Hemolysis and Viscoelasticity,PhD, 2005
    visc = m * pow((delta + rateofshear), (-1)*a);
  }
  else
    dserror("material type is not yet implemented");
    }



/*----------------------------------------------------------------------*
	 |  get the body force in the nodes of the element (private) gammi 04/07|
	 |  the Neumann condition associated with the nodes is stored in the    |
	 |  array edeadng only if all nodes have a VolumeNeumann condition      |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid2TH<distype>::BodyForce(Fluid2*      ele,    const double time)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique surface Neumann condition
  DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);

  if (myneumcond.size()>1)
    dserror("more than one SurfaceNeumann cond on one node");

  if (myneumcond.size()==1)
  {
    // find out whether we will use a time curve
    const vector<int>* curve  = myneumcond[0]->Get<vector<int> >("curve");
    int curvenum = -1;

    if (curve) curvenum = (*curve)[0];

    // initialisation
    double curvefac    = 0.0;

    if (curvenum >= 0) // yes, we have a timecurve
    {
      // time factor for the intermediate step
      if(time >= 0.0)
      {
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      }
      else
      {
        // A negative time value indicates an error.
        dserror("Negative time value in body force calculation: time = %f",time);
      }
    }
    else // we do not have a timecurve --- timefactors are constant equal 1
    {
      curvefac = 1.0;
    }

    // get values and switches from the condition
    const vector<int>*    onoff = myneumcond[0]->Get<vector<int> >   ("onoff");
    const vector<double>* val   = myneumcond[0]->Get<vector<double> >("val"  );
    const vector<int>*    functions = myneumcond[0]->Get<vector<int> >("funct");

    // factor given by spatial function
    double functionfac = 1.0;
    int functnum = -1;

    // set this condition to the edeadng array
    for (int jnode=0; jnode<iel; jnode++)
    {
      const double* x = (ele->Nodes()[jnode])->X();
      for(int isd=0;isd<2;isd++)
      {
        // get factor given by spatial function
        if (functions)
          functnum = (*functions)[isd];
        else
          functnum = -1;

        if (functnum>0)
        {
          // evaluate function at the position of the current node
          functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,x,time,NULL);
        }
        else
          functionfac = 1.0;

        // compute and store the (normalized) bodyforce value
        edeadng_(isd,jnode) = (*onoff)[isd]*(*val)[isd]*curvefac*functionfac;
      }
    }
  }
  else
  {
    // we have no dead load
    edeadng_.Clear();
  }
}


#endif
#endif /* D_FLUID2 */
