/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_impl.cpp

\brief Internal implementation of scalar transport boundary elements

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/

#if defined(D_FLUID3)
#ifdef CCADISCRET

// general Butler-Volmer is activated if the define-flag ButlerVolmer_Shifted is off
// the shifted version of Butler-Volmer is activated if the define-flag ButlerVolmer_Shifted is on
// define-flag PERCENT: how much is the curve shifted

//#define BUTLERVOLMER_SHIFTED
//#define PERCENT 0.01

#include <cstdlib>
#include "scatra_ele_boundary_impl.H"
#include "scatra_ele_impl.H"
#include "scatra_element.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_geometry/position_array.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
// material headers
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_spec.H"
#include "../drt_mat/arrhenius_temp.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::ScaTraBoundaryImplInterface* DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(
    const DRT::Element* ele,
    const enum INPAR::SCATRA::ScaTraType scatratype)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = ele->NumDofPerNode(*(ele->Nodes()[0]));
  int numscal = numdofpernode;
  if (scatratype == INPAR::SCATRA::scatratype_elch_enc)
    numscal -= 1;

  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return ScaTraBoundaryImpl<DRT::Element::quad4>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::quad8:
  {
    return ScaTraBoundaryImpl<DRT::Element::quad8>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::quad9:
  {
    return ScaTraBoundaryImpl<DRT::Element::quad9>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::tri3:
  {
    return ScaTraBoundaryImpl<DRT::Element::tri3>::Instance(numdofpernode,numscal);
  }
  /*  case DRT::Element::tri6:
  {
    return ScaTraBoundaryImpl<DRT::Element::tri6>::Instance(numdofpernode,numscal);
  }*/
  case DRT::Element::line2:
  {
    return ScaTraBoundaryImpl<DRT::Element::line2>::Instance(numdofpernode,numscal);
  }/*
  case DRT::Element::line3:
  {
    return ScaTraBoundaryImpl<DRT::Element::line3>::Instance(numdofpernode,numscal);
  }*/
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs2>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs3>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs4>::Instance(numdofpernode,numscal);
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return ScaTraBoundaryImpl<DRT::Element::nurbs9>::Instance(numdofpernode,numscal);
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryImpl<distype> * DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::Instance(
    const int numdofpernode,
    const int numscal
    )
{
  static ScaTraBoundaryImpl<distype> * instance;
  if ( instance==NULL )
    instance = new ScaTraBoundaryImpl<distype>(numdofpernode,numscal);
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ScaTraBoundaryImpl
(int numdofpernode,
    int numscal)
    : numdofpernode_(numdofpernode),
    numscal_(numscal),
    isale_(false),
    xyze_(true),  // initialize to zero
    edispnp_(true),
    diffus_(numscal_,0),
    valence_(numscal_,0),
    shcacp_(0.0),
    xsi_(true),
    funct_(true),
    deriv_(true),
    derxy_(true),
    normal_(true),
    velint_(true),
    metrictensor_(true)
    {
        return;
    }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::Evaluate(
    DRT::ELEMENTS::TransportBoundary* ele,
    ParameterList&                    params,
    DRT::Discretization&              discretization,
    vector<int>&                      lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseMatrix&         elemat2_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    Epetra_SerialDenseVector&         elevec2_epetra,
    Epetra_SerialDenseVector&         elevec3_epetra
)
{
  // First, do the things that are needed for all actions:

  // get the parent element including its material
  DRT::ELEMENTS::Transport* parentele = ele->ParentElement();
  RCP<MAT::Material> mat = parentele->Material();

  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = params.get<INPAR::SCATRA::ScaTraType>("scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // get node coordinates (we have a nsd_+1 dimensional domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // get additional state vector for ALE case: grid displacement
  isale_ = params.get<bool>("isale");
  if (isale_)
  {
    const RCP<Epetra_MultiVector> dispnp = params.get< RCP<Epetra_MultiVector> >("dispnp",null);
    if (dispnp==null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_+1);
    // add nodal displacements
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,myknots_,weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  // Now, check for the action parameter
  const string action = params.get<string>("action","none");
  if (action == "calc_normal_vectors")
  {
    // access the global vector
    const RCP<Epetra_MultiVector> normals = params.get< RCP<Epetra_MultiVector> >("normal vectors",null);
    if (normals == Teuchos::null) dserror("Could not access vector 'normal vectors'");

    // determine constant normal to this element
    GetConstNormal(normal_,xyze_);

    // loop over the element nodes
    for (int j=0;j<nen_;j++)
    {
      const int nodegid = (ele->Nodes()[j])->Id();
      if (normals->Map().MyGID(nodegid) )
      { // OK, the node belongs to this processor

        // scaling to a unit vector is performed on the global level after
        // assembly of nodal contributions since we have no reliable information
        // about the number of boundary elements adjacent to a node
        for (int dim=0; dim<(nsd_+1); dim++)
        {
          normals->SumIntoGlobalValue(nodegid,dim,normal_(dim));
        }
      }
      //else: the node belongs to another processor; the ghosted
      //      element will contribute the right value on that proc
    }
  }
  else if (action =="calc_elch_electrode_kinetics")
  {
    // get actual values of transported scalars
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // get current condition
    Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'ElectrodeKinetics'");

    // access parameters of the condition
    const std::string* kinetics = cond->Get<std::string>("kinetic model");
    const int    speciesid = cond->GetInt("species");
    if (speciesid < 1) dserror("species number is not >= 1");
    double       pot0 = cond->GetDouble("pot");
    const int    curvenum = cond->GetInt("curve");
    const double alphaa = cond->GetDouble("alpha_a");
    const double alphac = cond->GetDouble("alpha_c");
    double       i0 = cond->GetDouble("i0");
    if (i0 >EPS14) dserror("i0 is positive, ergo not pointing INTO the domain: %f",i0);
    const double gamma = cond->GetDouble("gamma");
    const double refcon = cond->GetDouble("refcon");
    if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
    const double dlcapacitance = cond->GetDouble("dlcap");
    if (dlcapacitance < -EPS12) dserror("double-layer capacitance is negative: %f",dlcapacitance);
    const double frt = params.get<double>("frt"); // = F/RT

    // get control parameter from parameter list
    bool iselch(true);
    if (scatratype != INPAR::SCATRA::scatratype_elch_enc)
      iselch = false;
    const bool   is_stationary = params.get<bool>("using stationary formulation");
    const double time = params.get<double>("total time");
    double       timefac = 1.0;
    double       alphaF  = 1.0;
    // find out whether we shell use a time curve and get the factor
    // this feature can be also used for stationary "pseudo time loops"
    if (curvenum>=0)
    {
      const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      // adjust potential at metal side accordingly
      pot0 *= curvefac;
    }

    const bool calc_status = params.get<bool>("calc_status",false);
    if (!calc_status)
    {
      if (not is_stationary)
      {
        // One-step-Theta:    timefac = theta*dt
        // BDF2:              timefac = 2/3 * dt
        // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
        timefac = params.get<double>("time factor");
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
        if (timefac < 0.0) dserror("time factor is negative.");
        // we multiply i0 with timefac. So we do not have to give down this paramater
        // and perform automatically the multiplication of matrix and rhs with timefac
        i0 *= timefac;
      }

# if 0
      // print all parameters read from the current condition
      cout<<"kinetic model  = "<<*kinetics<<endl;
      cout<<"react. species = "<<speciesid<<endl;
      cout<<"pot0(mod.)     = "<<pot0<<endl;
      cout<<"pot0n          = "<<pot0n<<endl;
      cout<<"curvenum       = "<<curvenum<<endl;
      cout<<"alpha_a        = "<<alphaa<<endl;
      cout<<"alpha_c        = "<<alphac<<endl;
      cout<<"i0(mod.)       = "<<i0<<endl;
      cout<<"gamma          = "<<gamma<<endl;
      cout<<"refcon         = "<<refcon<<endl;
      cout<<"F/RT           = "<<frt<<endl<<endl;
      cout<<"time factor    = "<<timefac<<endl;
      cout<<"alpha_F        = "<<alphaF<<endl;
#endif

      EvaluateElectrodeKinetics(
          ele,
          elemat1_epetra,
          elevec1_epetra,
          ephinp,
          mat,
          speciesid,
          *kinetics,
          pot0,
          alphaa,
          alphac,
          i0,
          frt,
          gamma,
          refcon,
          iselch
      );
    }
    else
    {
      // NOTE: add integral value only for elements which are NOT ghosted!
      if(ele->Owner() == discretization.Comm().MyPID())
      {
        // get actual values of transported scalars
        RefCountPtr<const Epetra_Vector> phidtnp = discretization.GetState("timederivative");
        if (phidtnp==null) dserror("Cannot get state vector 'ephidtnp'");
        // extract local values from the global vector
        vector<double> ephidtnp(lm.size());
        DRT::UTILS::ExtractMyValues(*phidtnp,ephidtnp,lm);

        double pot0hist(0.0);
        if (not is_stationary)
        {
          // One-step-Theta:    timefac = theta*dt
          // BDF2:              timefac = 2/3 * dt
          // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
          timefac = params.get<double>("time factor");
          alphaF = params.get<double>("alpha_F");
          timefac *= alphaF;
          if (timefac < 0.0) dserror("time factor is negative.");
          // access history of electrode potential
          if (dlcapacitance > EPS12)
          {
            pot0hist = cond->GetDouble("pothist");
          }
        }

        ElectrodeStatus(
            ele,
            params,
            ephinp,
            ephidtnp,
            *kinetics,
            speciesid,
            pot0,
            pot0hist,
            alphaa,
            alphac,
            i0,
            frt,
            gamma,
            refcon,
            iselch,
            timefac,
            dlcapacitance);}
    }
    // Warning: If ButlerVolmer_Shifted defines-flag is on
    #ifdef BUTLERVOLMER_SHIFTED
      if (time < 5*timefac)
        if (ele->Id()== 0)
          if (discretization.Comm().MyPID() == 0)
            cout<<" Warning: ButlerVolmer is shifted to the right about " << PERCENT*100 <<"% !! " << endl;
    #endif
  }
  else if (action =="calc_therm_press")
  {
    // we dont know the parent element's lm vector; so we have to build it here
    const int nenparent = parentele->NumNode();
    vector<int> lmparent(nenparent);
    vector<int> lmparentowner;
    vector<int> lmparentstride;
    parentele->LocationVector(discretization, lmparent, lmparentowner,lmparentstride);

    // get velocity values at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    // we deal with a (nsd_+1)-dimensional flow field
    Epetra_SerialDenseVector evel((nsd_+1)*nenparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

    // get values of scalar
    RefCountPtr<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vectors for the parent(!) element
    vector<double> myphinp(lmparent.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lmparent);

    // define vector for normal fluxes
    vector<double> mydiffflux(lm.size());
    vector<double> mydivu(lm.size());

    // determine normal to this element
    GetConstNormal(normal_,xyze_);

    // extract temperature flux vector for each node of the parent element
    LINALG::SerialDenseMatrix eflux(3,nenparent);
    DRT::Element* peleptr = (DRT::Element*) parentele;
    int k=numscal_-1;     // temperature is always last degree of freedom!!
    ostringstream temp;
    temp << k;
    string name = "flux_phi_"+temp.str();
    // try to get the pointer to the entry (and check if type is RCP<Epetra_MultiVector>)
    RCP<Epetra_MultiVector>* f = params.getPtr< RCP<Epetra_MultiVector> >(name);
    if (f!= NULL) // field has been set and is not of type Teuchos::null
    {
      DRT::UTILS::ExtractMyNodeBasedValues(peleptr,eflux,*f,3);
    }
    else
      dserror("MultiVector %s has not been found!",name.c_str());

    for (int i=0; i<nen_; ++i)
    {
      for(int j = 0; j<nenparent;++j)
      {
        // calculate normal diffusive flux and velocity div. at present node
        mydiffflux[i] = 0.0;
        mydivu[i]     = 0.0;
        for (int l=0; l<nsd_+1; l++)
        {
          mydiffflux[i] += eflux(l,j)*normal_(l);
          mydivu[i]     += evel[i*(nsd_+1)+l]*normal_(l);
        }
      }
    }

    // calculate integral of normal diffusive flux and velocity divergence
    // NOTE: add integral value only for elements which are NOT ghosted!
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      DifffluxAndDivuIntegral(ele,params,mydiffflux,mydivu);
    }
  }
  else if (action =="integrate_shape_functions")
  {
    // NOTE: add area value only for elements which are NOT ghosted!
    const bool addarea = (ele->Owner() == discretization.Comm().MyPID());
    IntegrateShapeFunctions(ele,params,elevec1_epetra,addarea);
  }
  else if (action =="calc_Neumann_inflow")
  {
    // get control parameters
    is_stationary_  = params.get<bool>("using stationary formulation");
    is_genalpha_    = params.get<bool>("using generalized-alpha time integration");
    is_incremental_ = params.get<bool>("incremental solver");

    // get time factor and alpha_F if required
    // one-step-Theta:    timefac = theta*dt
    // BDF2:              timefac = 2/3 * dt
    // generalized-alpha: timefac = alphaF * (gamma*/alpha_M) * dt
    double timefac = 1.0;
    double alphaF  = 1.0;
    if (not is_stationary_)
    {
      timefac = params.get<double>("time factor");
      if (is_genalpha_)
      {
        alphaF = params.get<double>("alpha_F");
        timefac *= alphaF;
      }
      if (timefac < 0.0) dserror("time factor is negative.");
    }

    // set thermodynamic pressure and its time derivative as well as
    // flag for turbulence model if required
    thermpress_ = 0.0;
    if (scatratype==INPAR::SCATRA::scatratype_loma)
      thermpress_ = params.get<double>("thermodynamic pressure");

    // we dont know the parent element's lm vector; so we have to build it here
    const int nenparent = parentele->NumNode();
    vector<int> lmparent(nenparent);
    vector<int> lmparentowner;
    vector<int> lmparentstride;
    parentele->LocationVector(discretization, lmparent, lmparentowner,lmparentstride);

    // get velocity values at nodes
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    // we deal with a (nsd_+1)-dimensional flow field
    Epetra_SerialDenseVector evel((nsd_+1)*nenparent);
    DRT::UTILS::ExtractMyNodeBasedValues(parentele,evel,velocity,nsd_+1);

    // get values of scalar
    RefCountPtr<const Epetra_Vector> phinp  = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vectors for the parent(!) element
    vector<double> myphinp(lmparent.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lmparent);

    // create object for density and solution array
    vector<LINALG::Matrix<nen_,1> > ephinp(numscal_);
    LINALG::Matrix<nsd_+1,nen_>     evelnp;

    // insert into element arrays
    for (int i=0;i<nen_;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
      }

      // insert velocity field into element array
      for (int idim=0 ; idim < nsd_+1 ; idim++)
      {
        evelnp(idim,i) = evel[idim + i*(nsd_+1)];
      }
    }

    NeumannInflow(ele,
                  mat,
                  ephinp,
                  evelnp,
                  elemat1_epetra,
                  elevec1_epetra,
                  timefac,
                  alphaF);

  }
  else if (action =="WeakDirichlet")
  {
    switch (distype)
    {
      if (numscal_>1) dserror("not yet implemented for more than one scalar\n");

      // 2D:
      case DRT::Element::line2:
      {
        if(ele->ParentElement()->Shape()==DRT::Element::quad4)
        {
          WeakDirichlet<DRT::Element::line2,DRT::Element::quad4>(ele,
                                                                 params,
                                                                 discretization,
                                                                 mat,
                                                                 elemat1_epetra,
                                                                 elevec1_epetra);
        }
        else
        {
          dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
        }
        break;
      }
      // 3D:
      case DRT::Element::quad4:
      {
        if(ele->ParentElement()->Shape()==DRT::Element::hex8)
        {
          WeakDirichlet<DRT::Element::quad4,DRT::Element::hex8>(ele,
                                                                params,
                                                                discretization,
                                                                mat,
                                                                elemat1_epetra,
                                                                elevec1_epetra);
        }
        else
        {
          dserror("expected combination quad4/hex8 or line2/quad4 for surface/parent pair");
        }
        break;
      }
      default:
      {
        dserror("not implemented yet\n");
      }
      }
  }
  else
    dserror("Unknown type of action for Scatra boundary impl.: %s",action.c_str());

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateNeumann(
    DRT::Element*             ele,
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  // get node coordinates (we have a nsd_+1 dimensional computational domain!)
  GEO::fillInitialPositionArray<distype,nsd_+1,LINALG::Matrix<nsd_+1,nen_> >(ele,xyze_);

  // get additional state vector for ALE case: grid displacement
  isale_ = params.get<bool>("isale");
  if (isale_)
  {
    const RCP<Epetra_MultiVector> dispnp = params.get< RCP<Epetra_MultiVector> >("dispnp",null);
    if (dispnp==null) dserror("Cannot get state vector 'dispnp'");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,edispnp_,dispnp,nsd_);
    // add nodal displacements to point coordinates
    xyze_ += edispnp_;
  }
  else edispnp_.Clear();

  // Now do the nurbs specific stuff (for isogeometric elements)
  if(DRT::NURBS::IsNurbs(distype))
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis
    = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    bool zero_size(false);
    // get local knot vector entries and check for zero sized elements
    zero_size = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots_,ele->Id());

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
    {
      return(0);
    }
    // you are still here? So get the node weights for nurbs elements as well
    DRT::Node** nodes = ele->Nodes();
    for (int inode=0; inode<nen_; inode++)
    {
      DRT::NURBS::ControlPoint* cp
      = dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
      weights_(inode) = cp->W();
    }
  } // Nurbs specific stuff

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

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
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );
  const vector<int>*    func  = condition.Get<vector<int> >   ("funct");

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());

    // multiply integration factor with the timecurve factor
    fac *= curvefac;

    // factor given by spatial function
    double functfac = 1.0;
    // determine global coordinates of current Gauss point
    double coordgp[nsd_];
    for (int i = 0; i< nsd_; i++)
    {
      coordgp[i] = 0.0;
      for (int j = 0; j < nen_; j++)
      {
        coordgp[i] += xyze_(i,j) * funct_(j);
      }
    }

    int functnum = -1;
    const double* coordgpref = &coordgp[0]; // needed for function evaluation


    for(int dof=0;dof<numdofpernode_;dof++)
    {
      if ((*onoff)[dof]) // is this dof activated?
      {
        // factor given by spatial function
        if (func) functnum = (*func)[dof];
        {
          if (functnum>0)
          {
            // evaluate function at current gauss point
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof,coordgpref,0.0,NULL);
          }
          else
            functfac = 1.0;
        }

        const double val_fac_functfac = (*val)[dof]*fac*functfac;

        for (int node=0;node<nen_;++node)
        {
          elevec1[node*numdofpernode_+dof] += funct_(node)*val_fac_functfac;
        }
      } // if ((*onoff)[dof])
    }
  } //end of loop over integration points

  return 0;
}


/*----------------------------------------------------------------------*
 | calculate Neumann inflow boundary conditions                vg 03/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::NeumannInflow(
    const DRT::Element*                   ele,
    Teuchos::RCP<const MAT::Material>     material,
    const vector<LINALG::Matrix<nen_,1> >& ephinp,
    const LINALG::Matrix<nsd_+1,nen_>&     evelnp,
    Epetra_SerialDenseMatrix&             emat,
    Epetra_SerialDenseVector&             erhs,
    const double                          timefac,
    const double                          alphaF)
{
  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // determine constant normal to this element
  GetConstNormal(normal_,xyze_);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,iquad,ele->Id());

    for(int k=0;k<numdofpernode_;++k)
    {
      // get velocity at integration point
      velint_.Multiply(evelnp,funct_);

      // normal velocity
      const double normvel = velint_.Dot(normal_);

      if (normvel<-0.0001)
      {
        // set density to 1.0
        double dens = 1.0;

        // get density if not constant
        if (material->MaterialType() == INPAR::MAT::m_matlist)
        {
          const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

          const int matid = actmat->MatID(0);
          Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

          if (singlemat->MaterialType() == INPAR::MAT::m_arrhenius_temp)
          {
            const MAT::ArrheniusTemp* actsinglemat = static_cast<const MAT::ArrheniusTemp*>(singlemat.get());

            // compute temperature
            const double temp = funct_.Dot(ephinp[k]);

            // compute density based on temperature and thermodynamic pressure
            dens = actsinglemat->ComputeDensity(temp,thermpress_);
          }
          else dserror("type of material found in material list is not supported");
        }
        else if (material->MaterialType() == INPAR::MAT::m_mixfrac)
        {
          const MAT::MixFrac* actmat = static_cast<const MAT::MixFrac*>(material.get());

          // compute mixture fraction
          const double mixfrac = funct_.Dot(ephinp[k]);

          // compute density based on mixture fraction
          dens = actmat->ComputeDensity(mixfrac);
        }
        else if (material->MaterialType() == INPAR::MAT::m_sutherland)
        {
          const MAT::Sutherland* actmat = static_cast<const MAT::Sutherland*>(material.get());

          // compute temperature
          const double temp = funct_.Dot(ephinp[k]);

          // compute density based on temperature and thermodynamic pressure
          dens = actmat->ComputeDensity(temp,thermpress_);
        }
        else if (material->MaterialType() == INPAR::MAT::m_arrhenius_pv)
        {
          const MAT::ArrheniusPV* actmat = static_cast<const MAT::ArrheniusPV*>(material.get());

          // compute progress variable
          const double provar = funct_.Dot(ephinp[k]);

          // compute density
          dens = actmat->ComputeDensity(provar);
        }
        else if (material->MaterialType() == INPAR::MAT::m_ferech_pv)
        {
          const MAT::FerEchPV* actmat = static_cast<const MAT::FerEchPV*>(material.get());

          // compute progress variable
          const double provar = funct_.Dot(ephinp[k]);

          // compute density
          dens = actmat->ComputeDensity(provar);
        }
        else dserror("Material type is not supported");

        // integration factor for left-hand side
        const double lhsfac = dens*normvel*timefac*fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if (is_incremental_ and is_genalpha_)
          rhsfac = lhsfac/alphaF;
        else if (not is_incremental_ and is_genalpha_)
          rhsfac = lhsfac*(1.0-alphaF)/alphaF;
        else if (is_incremental_ and not is_genalpha_)
          rhsfac = lhsfac;

        // matrix
        for (int vi=0; vi<nen_; ++vi)
        {
          const double vlhs = lhsfac*funct_(vi);

          const int fvi = vi*numdofpernode_+k;

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) -= vlhs*funct_(ui);
          }
        }

        // scalar at integration point
        const double phi = funct_.Dot(ephinp[k]);

        // rhs
        const double vrhs = rhsfac*phi;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          erhs[fvi] += vrhs*funct_(vi);
        }
      }
    }
  }

  return;

} //ScaTraBoundaryImpl<distype>::NeumannInflow


/*----------------------------------------------------------------------*
 | evaluate shape functions and int. factor at int. point     gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvalShapeFuncAndIntFac(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    const int                                    iquad,      ///< id of current Gauss point
    const int                                    eleid       ///< the element id
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
  {xsi_(idim) = gpcoord[idim];}

  if(not DRT::NURBS::IsNurbs(distype))
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  }
  else // nurbs elements are always somewhat special...
  {
    DRT::NURBS::UTILS::nurbs_get_funct_deriv(
        funct_  ,
        deriv_  ,
        xsi_    ,
        myknots_,
        weights_,
        distype );
  }

  // the metric tensor and the area of an infinitesimal surface/line element
  double drs(0.0);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs);

  // return the integration factor
  return intpoints.IP().qwgt[iquad] * drs;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition (private) gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::EvaluateElectrodeKinetics(
    const DRT::Element*        ele,
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const vector<double>&   ephinp,
    Teuchos::RCP<const MAT::Material> material,
    const int            speciesid,
    const std::string     kinetics,
    const double              pot0,
    const double            alphaa,
    const double            alphac,
    const double                i0,
    const double               frt,
    const double             gamma,
    const double            refcon,
    const bool              iselch
)
{
  //for pre-multiplication of i0 with 1/(F z_k)
  double fz = 1.0/96485.3399; // unit of F: C/mol or mC/mmol or muC / mumol

  // index of reactive species (starting from zero)
  const int k = speciesid-1;

  if (iselch) // this is not necessary for secondary current distributions
  {
    // get valence of the single(!) reactant
    if (material->MaterialType() == INPAR::MAT::m_matlist)
    {
      const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

      const int matid = actmat->MatID(k);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);
      if (singlemat->MaterialType() == INPAR::MAT::m_ion)
      {
        const MAT::Ion* actsinglemat = static_cast<const MAT::Ion*>(singlemat.get());
        fz = fz/actsinglemat->Valence();
      }
      else
        dserror("single material type is not 'ion'");
    }
    else
      dserror("material type is not a 'matlist' material");
  }

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  LINALG::Matrix<nen_,1> conreact(true);

  // el. potential values at element nodes
  LINALG::Matrix<nen_,1> pot(true);
  if(iselch)
  {
    for (int inode=0; inode< nen_;++inode)
    {
      conreact(inode) = ephinp[inode*numdofpernode_+k];
      pot(inode) = ephinp[inode*numdofpernode_+numscal_];
    }
  }
  else
  {
    for (int inode=0; inode< nen_;++inode)
    {
      conreact(inode) = 1.0;
      pot(inode) = ephinp[inode*numdofpernode_];
    }
  }

  // concentration of active species at integration point
  double conint(0.0);
  // el. potential at integration point
  double potint(0.0);
  // a 'working variable'
  double fac_fz_i0_funct_vi(0.0);


 /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // elch-specific values at integration point:
    conint = funct_.Dot(conreact);
    potint = funct_.Dot(pot);

    // surface overpotential eta at integration point
    const double eta = (pot0 - potint);

    if (iselch) // tertiary current distribution
    {
      // concentration-dependent Butler-Volmer law(s)
      if ((kinetics=="Butler-Volmer") or (kinetics=="Butler-Volmer-Yang1997"))
      {
        double pow_conint_gamma_k = 0.0;

// standard Butler-Volmer
#ifndef BUTLERVOLMER_SHIFTED
#ifdef DEBUG
        // some safety checks/ user warnings
        if ((alphaa*frt*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_a...) in Butler-Volmer law is near overflow!"
          <<exp(alphaa*frt*eta)<<endl;
        if (((-alphac)*frt*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
          <<exp((-alphac)*frt*eta)<<endl;
#endif
        if ((conint/refcon) < EPS13)
        {
          pow_conint_gamma_k = pow(EPS13,gamma);
#ifdef DEBUG
          cout<<"WARNING: Rel. Conc. in Butler-Volmer formula is zero/negative: "<<(conint/refcon)<<endl;
          cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<endl;
#endif
        }
        else
          pow_conint_gamma_k = pow(conint/refcon,gamma);
        if (kinetics=="Butler-Volmer")
        {
        // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
        const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

        double concterm = 0.0;
        if (conint > EPS13)
          concterm = gamma*pow(conint,(gamma-1.0))/pow(refcon,gamma);
        else
          concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

        for (int vi=0; vi<nen_; ++vi)
        {
          fac_fz_i0_funct_vi = fac*fz*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += fac_fz_i0_funct_vi*concterm*funct_(ui)*expterm;
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac_fz_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_+k] -= fac_fz_i0_funct_vi*pow_conint_gamma_k*expterm;
        }
        } // if (kinetics=="Butler-Volmer")
#else
// Butler-Volmer Shifted (Ehrl)
// c=0, gamma not 1 or 0: Butler-Volmer is not defined
// Concentration dependency pre-factor of the Butler-Volmer law is adapted by a linear part (Taylor series) below a chosen limit
// Therefore, the pre-factor curve is shifted to the right to ensure the functional value zero for c=0;
// it is possible to compute limiting currents in the case of a galvanostatic boundary condition without neg. concentrations
        double limit = refcon * PERCENT;
        double TaylorAbleitung = gamma*pow(limit,gamma-1)/pow(refcon,gamma);
        double Versatz = -pow(limit/refcon,gamma)/TaylorAbleitung+limit;

        if (conint < (limit-Versatz))
        {
          pow_conint_gamma_k = TaylorAbleitung*conint;
#ifdef DEBUG
          cout<<"WARNING: Rel. Conc. in Butler-Volmer formula is zero/negative: "<<(conint/refcon)<<endl;
          cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<endl;
#endif
        }
        else
          pow_conint_gamma_k = pow((conint+Versatz)/refcon,gamma);
          //pow_conint_gamma_k = pow(conint/refcon,gamma);
        if (kinetics=="Butler-Volmer")
        {
        // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
        const double expterm = exp(alphaa*frt*eta)-exp((-alphac)*frt*eta);

        double concterm = 0.0;
        if (conint > limit-Versatz)
          concterm = gamma*pow(conint+Versatz,(gamma-1.0))/pow(refcon,gamma);
        else
          concterm = TaylorAbleitung;

        for (int vi=0; vi<nen_; ++vi)
        {
          fac_fz_i0_funct_vi = fac*fz*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += fac_fz_i0_funct_vi*concterm*funct_(ui)*expterm;
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac_fz_i0_funct_vi*pow_conint_gamma_k*(((-alphaa)*frt*exp(alphaa*frt*eta))+((-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_+k] -= fac_fz_i0_funct_vi*pow_conint_gamma_k*expterm;
        }
        } // if (kinetics=="Butler-Volmer")
#endif
        if (kinetics=="Butler-Volmer-Yang1997")
        {
          // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
          double concterm = 0.0;
          if ((conint/refcon) > EPS13)
            concterm = gamma*pow(conint,(gamma-1.0))/pow(refcon,gamma);
          else
            concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

          for (int vi=0; vi<nen_; ++vi)
          {
            fac_fz_i0_funct_vi = fac*fz*i0*funct_(vi);
            // ---------------------matrix
            for (int ui=0; ui<nen_; ++ui)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += fac_fz_i0_funct_vi*funct_(ui)*(-(concterm*exp((-alphac)*frt*eta)));
              emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac_fz_i0_funct_vi*(((-alphaa)*frt*exp(alphaa*frt*eta))+(pow_conint_gamma_k*(-alphac)*frt*exp((-alphac)*frt*eta)))*funct_(ui);
            }
            // ------------right-hand-side
            erhs[vi*numdofpernode_+k] -= fac_fz_i0_funct_vi*(exp(alphaa*frt*eta)-(pow_conint_gamma_k*exp((-alphac)*frt*eta)));
          }
        } // if (kinetics=="Butler-Volmer-Yang1997")
      }
      else if(kinetics=="Tafel") // Tafel law (= remove anodic term from Butler-Volmer!)
      {
        // concentration-dependent Tafel law
        double pow_conint_gamma_k(0.0);

  #ifdef DEBUG
        // some safety checks/ user warnings
        if (((-alphac)*frt*eta) > 100.0)
            cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
            <<exp((-alphac)*frt*eta)<<endl;
  #endif
          if ((conint/refcon) < EPS13)
          {
            pow_conint_gamma_k = pow(EPS13,gamma);
  #ifdef DEBUG
            cout<<"WARNING: Rel. Conc. in Tafel formula is zero/negative: "<<(conint/refcon)<<endl;
            cout<<"-> Replacement value: pow(EPS,gamma) = "<< pow_conint_gamma_k <<endl;
  #endif
          }
          else
            pow_conint_gamma_k = pow(conint/refcon,gamma);

          // note: gamma==0 deactivates concentration dependency in Butler-Volmer!
          const double expterm = -exp((-alphac)*frt*eta);

          double concterm = 0.0;
          if (conint > EPS13)
            concterm = gamma*pow(conint,(gamma-1.0))/pow(refcon,gamma);
          else
            concterm = gamma*pow(EPS13,(gamma-1.0))/pow(refcon,gamma);

          for (int vi=0; vi<nen_; ++vi)
          {
            fac_fz_i0_funct_vi = fac*fz*i0*funct_(vi);
            // ---------------------matrix
            for (int ui=0; ui<nen_; ++ui)
            {
              emat(vi*numdofpernode_+k,ui*numdofpernode_+k) += fac_fz_i0_funct_vi*concterm*funct_(ui)*expterm;
              emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac_fz_i0_funct_vi*pow_conint_gamma_k*(-alphac)*frt*exp((-alphac)*frt*eta)*funct_(ui);
            }
            // ------------right-hand-side
            erhs[vi*numdofpernode_+k] -= fac_fz_i0_funct_vi*pow_conint_gamma_k*expterm;
          }
      }
      else if(kinetics=="linear") // linear law:  i_n = i_0*(alphaa*frt*(V_M - phi) + 1.0)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          fac_fz_i0_funct_vi = fac*fz*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += fac_fz_i0_funct_vi*(-alphaa)*frt*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_+k] -= fac_fz_i0_funct_vi*(alphaa*frt*eta + 1.0);
        }
      }
      else
        dserror("Kinetic model not implemented: %s",kinetics.c_str());
    }
    else // secondary current distribution
    {
      if (kinetics=="Butler-Volmer")  // concentration-dependent Butler-Volmer law
      {
#ifdef DEBUG
        // some safety checks/ user warnings
        if ((alphaa*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_a...) in Butler-Volmer law is near overflow!"
          <<exp(alphaa*eta)<<endl;
        if (((-alphac)*eta) > 100.0)
          cout<<"WARNING: Exp(alpha_c...) in Butler-Volmer law is near overflow!"
          <<exp((-alphac)*eta)<<endl;
#endif
        // Butler-Volmer kinetics
        const double expterm = exp(alphaa*eta)-exp((-alphac)*eta);
        const double exptermderiv = (((-alphaa)*exp(alphaa*eta))+((-alphac)*exp((-alphac)*eta)));

        for (int vi=0; vi<nen_; ++vi)
        {
          const double fac_i0_funct_vi = fac*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_,ui*numdofpernode_) += fac_i0_funct_vi*exptermderiv*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_] -= fac_i0_funct_vi*expterm;
        }
      }
      else if(kinetics == "Tafel") // Tafel kinetics
      {
#ifdef DEBUG
        // some safety checks/ user warnings
        if ((-alphac)*eta > 100.0)
          cout<<"WARNING: Exp(alpha_c...) in Tafel law is near overflow!"
          <<exp((-alphac)*eta)<<endl;
#endif
        const double expterm = -exp((-alphac)*eta);
        const double exptermderiv = alphac*expterm; // do not forget the (-1) from differentiation of eta!

        for (int vi=0; vi<nen_; ++vi)
        {
          const double fac_i0_funct_vi = fac*i0*funct_(vi);
          // ---------------------matrix
          for (int ui=0; ui<nen_; ++ui)
          {
            emat(vi*numdofpernode_,ui*numdofpernode_) += fac_i0_funct_vi*exptermderiv*funct_(ui);
          }
          // ------------right-hand-side
          erhs[vi*numdofpernode_] -= fac_i0_funct_vi*expterm;
        }
      }
      else
        dserror("Kinetic model not implemented: %s",kinetics.c_str());
    } // if iselch

  } // end of loop over integration points gpid

  return;

} // ScaTraBoundaryImpl<distype>::EvaluateElectrodeKinetics()


/*----------------------------------------------------------------------*
 | calculate electrode kinetics status information             gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::ElectrodeStatus(
    const DRT::Element*        ele,
    ParameterList&          params,
    const vector<double>&   ephinp,
    const vector<double>& ephidtnp,
    const std::string     kinetics,
    const int            speciesid,
    const double              pot0,
    const double          pot0hist,
    const double            alphaa,
    const double            alphac,
    const double                i0,
    const double               frt,
    const double             gamma,
    const double            refcon,
    const bool              iselch,
    const double           timefac,
    const double     dlcapacitance
)
{
  // get variables with their current values
  double currentintegral  = params.get<double>("currentintegral");
  double boundaryint      = params.get<double>("boundaryintegral");
  double overpotentialint = params.get<double>("overpotentialintegral");
  double concentrationint = params.get<double>("concentrationintegral");
  double currderiv        = params.get<double>("currentderiv");
  double currentresidual  = params.get<double>("currentresidual");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  LINALG::Matrix<nen_,1> conreact(true);

  // el. potential values at element nodes
  LINALG::Matrix<nen_,1> pot(true);
  LINALG::Matrix<nen_,1> potdtnp(true);
  if(iselch)
  {
    for (int inode=0; inode< nen_;++inode)
    {
      conreact(inode) = ephinp[inode*numdofpernode_+(speciesid-1)];
      pot(inode) = ephinp[inode*numdofpernode_+numscal_];
      potdtnp(inode) = ephidtnp[inode*numdofpernode_+numscal_];
    }
  }
  else
  {
    for (int inode=0; inode< nen_;++inode)
    {
      conreact(inode) = 1.0;
      pot(inode) = ephinp[inode*numdofpernode_];
      potdtnp(inode) = ephidtnp[inode*numdofpernode_];
    }
  }

  // concentration of active species at integration point
  double conint;
  // el. potential at integration point
  double potint;
  // history term of el. potential at integration point
  double potdtnpint;

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // elch-specific values at integration point:
    conint = funct_.Dot(conreact);
    potint = funct_.Dot(pot);
    potdtnpint = funct_.Dot(potdtnp);

    // surface overpotential eta at integration point
    const double eta = (pot0 - potint);

    // linearization of current w.r.t applied electrode potential "pot0"
    double linea(0.0);

    if ((kinetics=="Butler-Volmer") or (kinetics=="Butler-Volmer-Yang1997"))
    {
      double expterm(0.0);
      if (iselch)
      {
        if (kinetics=="Butler-Volmer")
        {
// general Butler-Volmer
#ifndef BUTLERVOLMER_SHIFTED
          expterm = pow(conint/refcon,gamma) * (exp(alphaa*frt*eta)-exp((-alphac)*frt*eta));
          linea = pow(conint/refcon,gamma) * frt*((alphaa*exp(alphaa*frt*eta)) + (alphac*exp((-alphac)*frt*eta)));
#else
// Butler-Volmer Shifted: details see in function EvaluateElectrodeKinetics()
          double limit = refcon*PERCENT;
          double TaylorAbleitung = gamma*pow(limit,gamma-1)/pow(refcon,gamma);
          double Versatz = -pow(limit/refcon,gamma)/TaylorAbleitung+limit;

          double concterm = 0.0;
          if (conint > limit-Versatz)
            concterm = pow((conint+Versatz)/refcon,gamma);
          else
            concterm = TaylorAbleitung*conint;

          expterm = concterm * (exp(alphaa*frt*eta)-exp((-alphac)*frt*eta));
          linea = concterm * frt*((alphaa*exp(alphaa*frt*eta)) + (alphac*exp((-alphac)*frt*eta)));
#endif
        }
        if (kinetics=="Butler-Volmer-Yang1997")
        {
          if (((conint/refcon)<EPS13) && (gamma < 1.0))
          {// prevents NaN's in the current density evaluation
            expterm = (exp(alphaa*frt*eta)-(pow(EPS13/refcon,gamma)*exp((-alphac)*frt*eta)));
            linea = ((alphaa)*frt*exp(alphaa*frt*eta))+(pow(EPS13/refcon,gamma)*alphac*frt*exp((-alphac)*frt*eta));
          }
          else
          {
            expterm = (exp(alphaa*frt*eta)-(pow(conint/refcon,gamma)*exp((-alphac)*frt*eta)));
            linea = ((alphaa)*frt*exp(alphaa*frt*eta))+(pow(conint/refcon,gamma)*alphac*frt*exp((-alphac)*frt*eta));
          }
        }
      }
      else // secondary current distribution
      {
        expterm = exp(alphaa*eta)-exp((-alphac)*eta);
        linea = (alphaa*exp(alphaa*eta))+(alphac*exp((-alphac)*eta));
      }

      // scan for NaNs due to negative concentrations under exponent gamma
      if (std::isnan(expterm) or std::isnan(linea))
        dserror("NaN detected in electrode status calculation");

      // compute integrals
      overpotentialint += eta * fac;
      currentintegral += (-i0) * expterm * fac; // the negative(!) normal flux density
      boundaryint += fac;
      concentrationint += conint*fac;

      // tangent and rhs (= negative residual) for galvanostatic equation
      currderiv += i0*linea*timefac*fac;
      currentresidual += (-i0) * expterm * timefac *fac;

      if (dlcapacitance > EPS12)
      {
        // add contributions due to double-layer capacitance
        currderiv -= fac*dlcapacitance;
        currentresidual += fac*dlcapacitance*(pot0-pot0hist-(timefac*potdtnpint));
#if 0
        cout<<"pot0-pot0hist = "<<pot0 << " - "<<pot0hist<<endl;
        cout<<"- timefac*potdtnpint = -"<<timefac<<" * "<<potdtnpint<<endl;
#endif
      }

    }
    else if (iselch && (kinetics=="linear")) // linear: i_n = i_0*(alphaa*frt*eta + 1.0)
    {
      // compute integrals
      overpotentialint += eta * fac;
      currentintegral += (-i0) * (alphaa*frt*eta + 1.0) * fac; // the negative(!) normal flux density
      boundaryint += fac;
      concentrationint += conint*fac;
    }
    else if (iselch && (kinetics=="Tafel")) // concentration-dependent Tafel kinetics
    {
      const double expterm = pow(conint/refcon,gamma) * (-exp((-alphac)*frt*eta));
      linea = pow(conint/refcon,gamma) * frt*(alphac*exp((-alphac)*frt*eta));
      // compute integrals
      overpotentialint += eta * fac;
      currentintegral += (-i0) * expterm * fac; // the negative(!) normal flux density
      boundaryint += fac;
      concentrationint += conint*fac;
    }
    else if ((!iselch) && (kinetics=="Tafel"))
    {
      // secondary current distribution with Tafel kinetics
      double expterm = -exp((-alphac)*eta);
      //linea = (-alphac)*exp((-alphac)*eta);

      // compute integrals
      overpotentialint += eta * fac;
      currentintegral += (-i0) * expterm * fac; // the negative(!) normal flux density
      boundaryint += fac;
      //concentrationint += conint*fac_;
    }
    else
      dserror("Kinetic model not implemented: %s",kinetics.c_str());

  } // loop over integration points

  // add contributions to the global values
  params.set<double>("currentintegral",currentintegral);
  params.set<double>("boundaryintegral",boundaryint);
  params.set<double>("overpotentialintegral",overpotentialint);
  params.set<double>("concentrationintegral",concentrationint);
  params.set<double>("currentderiv",currderiv);
  params.set<double>("currentresidual",currentresidual);

} //ScaTraBoundaryImpl<distype>::ElectrodeStatus


/*----------------------------------------------------------------------*
 | get constant normal                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::GetConstNormal(
    LINALG::Matrix<nsd_+1,1>&          normal,
    const LINALG::Matrix<nsd_+1,nen_>&  xyze
    )
{
  // determine normal to this element
  switch(nsd_)
  {
  case 2:
  {
    LINALG::Matrix<3,1> dist1(true), dist2(true);
    for (int i=0; i<3; i++)
    {
      dist1(i) = xyze(i,1)-xyze(i,0);
      dist2(i) = xyze(i,2)-xyze(i,0);
    }

    normal(0) = dist1(1)*dist2(2) - dist1(2)*dist2(1);
    normal(1) = dist1(2)*dist2(0) - dist1(0)*dist2(2);
    normal(2) = dist1(0)*dist2(1) - dist1(1)*dist2(0);
  }
  break;
  case 1:
  {
    normal(0) = xyze(1,1) - xyze(1,0);
    normal(1) = (-1.0)*(xyze(0,1) - xyze(0,0));
  }
  break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_);
  } // switch(nsd)

  // length of normal to this element
  const double length = normal.Norm2();
  // outward-pointing normal of length 1.0
  normal.Scale(1/length);

  return;
} // ScaTraBoundaryImpl<distype>::


/*----------------------------------------------------------------------*
 |  Integrate shapefunctions over surface (private)           gjb 02/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::IntegrateShapeFunctions(
    const DRT::Element*        ele,
    ParameterList&             params,
    Epetra_SerialDenseVector&  elevec1,
    const bool                 addarea
)
{
  // access boundary area variable with its actual value
  double boundaryint = params.get<double>("boundaryint");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of shape functions
    for (int node=0;node<nen_;++node)
    {
      for (int k=0; k< numscal_; k++)
      {
        elevec1[node*numdofpernode_+k] += funct_(node) * fac;
      }
    }

    if (addarea)
    {
      // area calculation
      boundaryint += fac;
    }

  } //loop over integration points

  // add contribution to the global value
  params.set<double>("boundaryint",boundaryint);

  return;

} //ScaTraBoundaryImpl<distype>::IntegrateShapeFunction


/*----------------------------------------------------------------------*
 | calculate integral of normal diffusive flux + velocity div.  vg 09/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::DifffluxAndDivuIntegral(
    const DRT::Element*             ele,
    ParameterList&                  params,
    const vector<double>&           ediffflux,
    const vector<double>&           edivu
)
{
  // get variables for integrals of normal diffusive flux and velocity div.
  double difffluxintegral = params.get<double>("diffusive-flux integral");
  double divuintegral     = params.get<double>("velocity-divergence integral");

  // integrations points and weights
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());

    // compute integral of normal flux
    for (int node=0;node<nen_;++node)
    {
      difffluxintegral += funct_(node) * ediffflux[node] * fac;
      divuintegral     += funct_(node) * edivu[node] * fac;
    }
  } // loop over integration points

  // add contributions to the global values
  params.set<double>("diffusive-flux integral",difffluxintegral);
  params.set<double>("velocity-divergence integral",divuintegral);

  return;

} //ScaTraBoundaryImpl<distype>::DifffluxAndDivuIntegral


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType bdistype,
          DRT::Element::DiscretizationType pdistype>
   void  DRT::ELEMENTS::ScaTraBoundaryImpl<distype>::WeakDirichlet(
     DRT::ELEMENTS::TransportBoundary*  ele,
     ParameterList&                     params,
     DRT::Discretization&               discretization,
     Teuchos::RCP<const MAT::Material>  material,
     Epetra_SerialDenseMatrix&          elemat_epetra,
     Epetra_SerialDenseVector&          elevec_epetra)
{
  //------------------------------------------------------------------------
  // control parameters for time integration
  //------------------------------------------------------------------------
  is_stationary_  = params.get<bool>("using stationary formulation");
  is_incremental_ = params.get<bool>("incremental solver");

  // get time factor and alpha_F if required
  // one-step-Theta:    timefac = theta*dt
  // BDF2:              timefac = 2/3 * dt
  // generalized-alpha: timefac = alphaF * (gamma/alpha_M) * dt
  double timefac = 1.0;
  double alphaF  = 1.0;
  if (not is_stationary_)
  {
    timefac = params.get<double>("time factor");
    if (is_genalpha_)
    {
      alphaF = params.get<double>("alpha_F");
      timefac *= alphaF;
    }
    if (timefac < 0.0) dserror("time factor is negative.");
  }

  //------------------------------------------------------------------------
  // Dirichlet boundary condition
  //------------------------------------------------------------------------
  RCP<DRT::Condition> dbc = params.get<RCP<DRT::Condition> >("condition");

  // check of total time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = (*dbc).Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // get values and spatial functions from condition
  // (assumed to be constant on element boundary)
  const vector<double>* val  = (*dbc).Get<vector<double> >("val"  );
  const vector<int>*    func = (*dbc).Get<vector<int> >   ("funct");

  // assign boundary value multiplied by time-curve factor
  double dirichval=(*val)[0]*curvefac;

  // spatial function number
  const int funcnum = (*func)[0];

  //------------------------------------------------------------------------
  // preliminary definitions for (boundary) and parent element and
  // evaluation of nodal values of velocity and scalar based on parent
  // element nodes
  //------------------------------------------------------------------------
  // get the parent element
  DRT::ELEMENTS::Transport* pele = ele->ParentElement();

  // number of spatial dimensions regarding (boundary) element
  static const int bnsd = DRT::UTILS::DisTypeToDim<bdistype>::dim;

  // number of spatial dimensions regarding parent element
  static const int pnsd = DRT::UTILS::DisTypeToDim<pdistype>::dim;

  // number of (boundary) element nodes
  static const int bnen = DRT::UTILS::DisTypeToNumNodePerEle<bdistype>::numNodePerElement;

  // number of parent element nodes
  static const int pnen = DRT::UTILS::DisTypeToNumNodePerEle<pdistype>::numNodePerElement;

  // parent element lm vector (vectors plm and plmowner allocated outside in
  // EvaluateConditionUsingParentData)
  RCP<vector<int> > plm      = params.get<RCP<vector<int> > >("plm");
  RCP<vector<int> > plmowner = params.get<RCP<vector<int> > >("plmowner");
  RCP<vector<int> > plmstride = params.get<RCP<vector<int> > >("plmstride");
  pele->LocationVector(discretization,*plm,*plmowner,*plmstride);

  // get velocity values at parent element nodes
  const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
  Epetra_SerialDenseVector evel(pnsd*pnen);
  DRT::UTILS::ExtractMyNodeBasedValues(pele,evel,velocity,pnsd);

  // get scalar values at parent element nodes
  RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp==null) dserror("Cannot get state vector 'phinp'");

  // extract local values from global vectors for parent element
  vector<double> myphinp(plm->size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,*plm);

  // matrix and vector definition
  LINALG::Matrix<pnsd,pnen>       evelnp;
  vector<LINALG::Matrix<pnen,1> > ephinp(numscal_);

  // insert into element arrays
  for (int i=0;i<pnen;++i)
  {
    for (int k = 0; k< numscal_; ++k)
    {
      // split for each tranported scalar, insert into element arrays
      ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
    }

    // insert velocity field into element array
    for (int idim=0 ; idim < pnsd; idim++)
    {
      evelnp(idim,i) = evel[idim + i*pnsd];
    }
  }

  //------------------------------------------------------------------------
  // preliminary definitions for integration loop
  //------------------------------------------------------------------------
  // reshape element matrices and vectors and init to zero, construct views
  elemat_epetra.Shape(pnen,pnen);
  elevec_epetra.Size (pnen);
  LINALG::Matrix<pnen,pnen> emat(elemat_epetra.A(),true);
  LINALG::Matrix<pnen,   1> erhs(elevec_epetra.A(),true);

  // (boundary) element local node coordinates
  LINALG::Matrix<pnsd,bnen>  bxyze(true);
  GEO::fillInitialPositionArray<bdistype,pnsd,LINALG::Matrix<pnsd,bnen> >(ele,bxyze);

  // parent element local node coordinates
  LINALG::Matrix<pnsd,pnen>  pxyze(true);
  GEO::fillInitialPositionArray<pdistype,pnsd,LINALG::Matrix<pnsd,pnen> >(pele,pxyze);

  // coordinates of integration points for (boundary) and parent element
  LINALG::Matrix<bnsd,   1>  bxsi(true);
  LINALG::Matrix<pnsd,   1>  pxsi(true);

  // transposed jacobian "dx/ds" and inverse of transposed jacobian "ds/dx"
  // for parent element
  LINALG::Matrix<pnsd,pnsd>  pxjm(true);
  LINALG::Matrix<pnsd,pnsd>  pxji(true);

  // metric tensor for (boundary) element
  LINALG::Matrix<bnsd,bnsd>  bmetrictensor(true);

  // (outward-pointing) unit normal vector to (boundary) element
  LINALG::Matrix<pnsd,   1>  bnormal(true);

  // velocity vector at integration point
  LINALG::Matrix<pnsd,   1>  velint;

  // gradient of scalar value at integration point
  LINALG::Matrix<pnsd,1> gradphi;

  // (boundary) element shape functions, local and global derivatives
  LINALG::Matrix<bnen,   1>  bfunct(true);
  LINALG::Matrix<bnsd,bnen>  bderiv(true);
  LINALG::Matrix<bnsd,bnen>  bderxy(true);

  // parent element shape functions, local and global derivatives
  LINALG::Matrix<pnen,   1>  pfunct(true);
  LINALG::Matrix<pnsd,pnen>  pderiv(true);
  LINALG::Matrix<pnsd,pnen>  pderxy(true);

  //------------------------------------------------------------------------
  // additional matrices and vectors for mixed-hybrid formulation
  //------------------------------------------------------------------------
  // for volume integrals
  LINALG::Matrix<pnsd*pnen,pnsd*pnen> mat_s_q(true);
  LINALG::Matrix<pnsd*pnen,     pnen> mat_s_gradphi(true);

  LINALG::Matrix<pnsd*pnen,        1> vec_s_gradphi(true);

  // for boundary integrals
  LINALG::Matrix<     pnen,pnsd*pnen> mat_w_q_o_n(true);
  LINALG::Matrix<pnsd*pnen,     pnen> mat_s_o_n_phi(true);

  LINALG::Matrix<pnsd*pnen,        1> vec_s_o_n_phi_minus_g(true);

  // inverse matrix
  LINALG::Matrix<pnsd*pnen,pnsd*pnen> inv_s_q(true);

  //------------------------------------------------------------------------
  // check whether Nitsche (default) or mixed-hybrid formulation as well as
  // preliminary definitions and computations for Nitsche stabilization term
  //------------------------------------------------------------------------
  // default is Nitsche formulation
  bool mixhyb = false;

  // stabilization parameter for Nitsche term
  const double nitsche_stab_para = (*dbc).GetDouble("TauBscaling");

  // if stabilization parameter negative: mixed-hybrid formulation
  if (nitsche_stab_para < 0.0) mixhyb = true;

  // pre-factor for adjoint-consistency term:
  // either 1.0 (adjoint-consistent, default) or -1.0 (adjoint-inconsistent)
  double gamma = 1.0;
  const string* consistency = (*dbc).Get<string>("Choice of gamma parameter");
  if      (*consistency=="adjoint-consistent") gamma = 1.0;
  else if (*consistency=="diffusive-optimal")  gamma = -1.0;
  else dserror("unknown definition for gamma parameter: %s",(*consistency).c_str());

  // use one-point Gauss rule to do calculations at element center
  DRT::UTILS::IntPointsAndWeights<bnsd> intpoints_tau(SCATRA::DisTypeToStabGaussRule<bdistype>::rule);

  // element surface area (1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double* gpcoord = (intpoints_tau.IP().qxg)[0];
  for (int idim=0;idim<bnsd;idim++)
  {
    bxsi(idim) = gpcoord[idim];
  }
  DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);
  double drs = 0.0;
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);
  const double area = intpoints_tau.IP().qwgt[0]*drs;

  // get number of dimensions for (boundary) element (convert from int to double)
  const double dim = (double) bnsd;

  // computation of characteristic length of (boundary) element
  // (2D: square root of element area, 1D: element length)
  const double h = pow(area,(1.0/dim));

  //------------------------------------------------------------------------
  // preliminary computations for integration loop
  //------------------------------------------------------------------------
  // integrations points and weights for (boundary) element and parent element
  const DRT::UTILS::IntPointsAndWeights<bnsd> bintpoints(SCATRA::DisTypeToOptGaussRule<bdistype>::rule);

  const DRT::UTILS::IntPointsAndWeights<pnsd> pintpoints(SCATRA::DisTypeToOptGaussRule<pdistype>::rule);

  // transfer integration-point coordinates of (boundary) element to parent element
  Epetra_SerialDenseMatrix pqxg(pintpoints.IP().nquad,pnsd);
  {
    Epetra_SerialDenseMatrix gps(bintpoints.IP().nquad,bnsd);

    for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
    {
      const double* gpcoord = (bintpoints.IP().qxg)[iquad];

      for (int idim=0;idim<bnsd ;idim++)
      {
        gps(iquad,idim) = gpcoord[idim];
      }
    }
    DRT::UTILS::BoundaryGPToParentGP<pnsd>(pqxg,gps,pdistype,bdistype,ele->BeleNumber());
  }

  //------------------------------------------------------------------------
  // integration loop 1: volume integrals (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    for (int iquad=0; iquad<pintpoints.IP().nquad; ++iquad)
    {
      // reference coordinates of integration point from (boundary) element
      const double* gpcoord = (pintpoints.IP().qxg)[iquad];
      for (int idim=0;idim<pnsd;idim++)
      {
        pxsi(idim) = gpcoord[idim];
      }

      // parent element shape functions and local derivatives
      DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
      DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

      // Jacobian matrix and determinant of parent element (including check)
      pxjm.MultiplyNT(pderiv,pxyze);
      const double det = pxji.Invert(pxjm);
      if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

      // compute integration factor
      const double fac = pintpoints.IP().qwgt[iquad]*det;

      // compute global derivatives
      pderxy.Multiply(pxji,pderiv);

      //--------------------------------------------------------------------
      // loop over scalars (not yet implemented for more than one scalar)
      //--------------------------------------------------------------------
      // for(int k=0;k<numdofpernode_;++k)
      int k=0;
      {
        // get viscosity
        if (material->MaterialType() == INPAR::MAT::m_scatra)
        {
          const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

          dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

          // get constant diffusivity
          diffus_[k] = actmat->Diffusivity();
        }
        else dserror("Material type is not supported");

        // gradient of current scalar value
        gradphi.Multiply(pderxy,ephinp[k]);

        // integration factor for left-hand side
        const double lhsfac = timefac*fac;

        // integration factor for right-hand side
        double rhsfac = 0.0;
        if (is_incremental_ and is_genalpha_)
          rhsfac = lhsfac/alphaF;
        else if (not is_incremental_ and is_genalpha_)
          rhsfac = lhsfac*(1.0-alphaF)/alphaF;
        else if (is_incremental_ and not is_genalpha_)
          rhsfac = lhsfac;

        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*
                       /         \
                  1   |   h   h  |
              - ----- |  s , q   |
                kappa |          |
                      \          / Omega
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          //const double vlhs = lhsfac*pfunct(vi);
          const double vlhs = lhsfac*(1.0/diffus_[k])*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_q(fvi*pnsd+i,fui*pnsd+i) -= vlhs*pfunct(ui);
            }
          }
        }

        /*
                       /                  \
                      |  h         /   h\  |
                    + | s  , grad | phi  | |
                      |            \    /  |
                       \                  / Omega
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          //const double vlhs = lhsfac*diffus_[k]*pfunct(vi);
          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_gradphi(fvi*pnsd+i,fui) += vlhs*pderxy(i,ui);
            }
          }

          //const double vrhs = rhsfac*diffus_[k]*pfunct(vi);
          const double vrhs = rhsfac*pfunct(vi);

          for(int i=0;i<pnsd;++i)
          {
            vec_s_gradphi(fvi*pnsd+i) += vrhs*gradphi(i);
          }
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // integration loop 2: boundary integrals
  //------------------------------------------------------------------------
  for (int iquad=0; iquad<bintpoints.IP().nquad; ++iquad)
  {
    // reference coordinates of integration point from (boundary) element
    const double* gpcoord = (bintpoints.IP().qxg)[iquad];
    for (int idim=0;idim<bnsd;idim++)
    {
      bxsi(idim) = gpcoord[idim];
    }

    // (boundary) element shape functions
    DRT::UTILS::shape_function       <bdistype>(bxsi,bfunct);
    DRT::UTILS::shape_function_deriv1<bdistype>(bxsi,bderiv);

    // global coordinates of current integration point from (boundary) element
    LINALG::Matrix<pnsd,1> coordgp(true);
    for (int A=0;A<bnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        coordgp(j)+=bxyze(j,A)*bfunct(A);
      }
    }

    // reference coordinates of integration point from parent element
    for (int idim=0;idim<pnsd;idim++)
    {
      pxsi(idim) = pqxg(iquad,idim);
    }

    // parent element shape functions and local derivatives
    DRT::UTILS::shape_function       <pdistype>(pxsi,pfunct);
    DRT::UTILS::shape_function_deriv1<pdistype>(pxsi,pderiv);

    // Jacobian matrix and determinant of parent element (including check)
    pxjm.MultiplyNT(pderiv,pxyze);
    const double det = pxji.Invert(pxjm);
    if (det < 1E-16) dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", pele->Id(), det);

    // compute measure tensor for surface element, infinitesimal area element drs
    // and (outward-pointing) unit normal vector
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<bdistype>(bxyze,bderiv,bmetrictensor,drs,&bnormal);

    // compute integration factor
    const double fac = bintpoints.IP().qwgt[iquad]*drs;

    // compute global derivatives
    pderxy.Multiply(pxji,pderiv);

#if 1
    //--------------------------------------------------------------------
    // check whether integration-point coordinates evaluated from
    // (boundary) and parent element match
    //--------------------------------------------------------------------
    LINALG::Matrix<pnsd,1> check(true);
    LINALG::Matrix<pnsd,1> diff(true);

    for (int A=0;A<pnen;++A)
    {
      for(int j=0;j<pnsd;++j)
      {
        check(j)+=pxyze(j,A)*pfunct(A);
      }
    }

    diff=check;
    diff-=coordgp;

    const double norm=diff.Norm2();

    if (norm>1e-9)
    {
      for (int j=0;j<pnsd;++j)
      {
        printf("%12.5e %12.5e\n",check(j),coordgp(j));
      }
      dserror("Gausspoint matching error %12.5e\n",norm);
    }
#endif

    //--------------------------------------------------------------------
    // factor for Dirichlet boundary condition given by spatial function
    //--------------------------------------------------------------------
    double functfac = 1.0;
    if (funcnum > 0)
    {
      // evaluate function at current integration point
      functfac = DRT::Problem::Instance()->Funct(funcnum-1).Evaluate(0,coordgp.A(),0.0,NULL);
    }
    else functfac = 1.0;
    dirichval *= functfac;

    //--------------------------------------------------------------------
    // loop over scalars (not yet implemented for more than one scalar)
    //--------------------------------------------------------------------
    // for(int k=0;k<numdofpernode_;++k)
    int k=0;
    {
      // get viscosity
      if (material->MaterialType() == INPAR::MAT::m_scatra)
      {
        const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());

        dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");

        // get constant diffusivity
        diffus_[k] = actmat->Diffusivity();
      }
      else dserror("Material type is not supported");

      // get scalar value at integration point
      const double phi = pfunct.Dot(ephinp[k]);

      // integration factor for left-hand side
      const double lhsfac = timefac*fac;

      // integration factor for right-hand side
      double rhsfac = 0.0;
      if (is_incremental_ and is_genalpha_)
        rhsfac = lhsfac/alphaF;
      else if (not is_incremental_ and is_genalpha_)
        rhsfac = lhsfac*(1.0-alphaF)/alphaF;
      else if (is_incremental_ and not is_genalpha_)
        rhsfac = lhsfac;

      if (mixhyb)
      {
        //--------------------------------------------------------------------
        //  matrix and vector additions due to mixed-hybrid formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /           \
                   |  h   h     |
                 - | w , q  o n |
                   |            |
                   \            / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_w_q_o_n(fvi,fui*pnsd+i) -= vlhs*pfunct(ui)*bnormal(i);
            }
          }
        }

        /*  adjoint consistency term
                    /                 \
                   |  h          h    |
                 - | s  o n , phi - g |
                   |                  |
                   \                  / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              mat_s_o_n_phi(fvi*pnsd+i,fui) -= vlhs*pfunct(ui)*bnormal(i);
            }
          }

          for(int i=0;i<pnsd;++i)
          {
            vec_s_o_n_phi_minus_g(fvi*pnsd+i) -= pfunct(vi)*bnormal(i)*(rhsfac*phi - timefac*fac*dirichval);
          }
        }
      }
      else
      {
        // parameter alpha for Nitsche stabilization term
        const double alpha = nitsche_stab_para*diffus_[k]/h;

        // get velocity at integration point
        velint.Multiply(evelnp,pfunct);

        // normal velocity
        const double normvel = velint.Dot(bnormal);

        // gradient of current scalar value
        gradphi.Multiply(pderxy,ephinp[k]);

        // gradient of current scalar value in normal direction
        const double gradphi_norm = bnormal.Dot(gradphi);

        //--------------------------------------------------------------------
        //  matrix and vector additions due to Nitsche formulation
        //--------------------------------------------------------------------
        /*  consistency term
                    /                           \
                   |  h                  h      |
                 - | w , kappa * grad(phi ) o n |
                   |                            |
                   \                            / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double vlhs = lhsfac*pfunct(vi)*diffus_[k];

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            for(int i=0;i<pnsd;++i)
            {
              emat(fvi,fui) -= vlhs*pderxy(i,ui)*bnormal(i);
            }
          }

          const double vrhs = rhsfac*diffus_[k];

          erhs(fvi) += vrhs*pfunct(vi)*gradphi_norm;
        }

        /*  adjoint consistency term, inflow/outflow part
              / --          --                                        \
             |  |         h  |                      h           h     |
           - |  |(a o n) w  +| gamma * kappa *grad(w ) o n , phi - g  |
             |  |            |                                        |
             \  --          --                                        / Gamma_in/out
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          // compute diffusive part
          double prefac = 0.0;
          for(int i=0;i<pnsd;++i)
          {
            prefac += gamma*diffus_[k]*pderxy(i,vi)*bnormal(i);
          }

          // add convective part in case of inflow boundary
          if (normvel<-0.0001) prefac += normvel*pfunct(vi);

          const double vlhs = lhsfac*prefac;

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) -= vlhs*pfunct(ui);
          }

          erhs(fvi) += prefac*(rhsfac*phi - timefac*fac*dirichval);
        }

        /*  stabilization term
                            /             \
                           |  h     h     |
                 + alpha * | w , phi - g  |
                           |              |
                           \              / Gamma
        */
        for (int vi=0; vi<pnen; ++vi)
        {
          const int fvi = vi*numdofpernode_+k;

          const double prefac = alpha*pfunct(vi);

          for (int ui=0; ui<pnen; ++ui)
          {
            const int fui = ui*numdofpernode_+k;

            emat(fvi,fui) += lhsfac*prefac*pfunct(ui);
          }

          erhs(fvi) -= prefac*(rhsfac*phi - timefac*fac*dirichval);
        }
      }
    }
  }

  //------------------------------------------------------------------------
  // local condensation (only for mixed-hybrid formulation)
  //------------------------------------------------------------------------
  if (mixhyb)
  {
    // matrix inversion of flux-flux block
    inv_s_q = mat_s_q;

    LINALG::FixedSizeSerialDenseSolver<pnsd*pnen,pnsd*pnen> solver;

    solver.SetMatrix(inv_s_q);
    solver.Invert();

    // computation of matrix-matrix and matrix vector products, local assembly
    for (int vi=0; vi<pnen; ++vi)
    {
      for (int ui=0; ui<pnen; ++ui)
      {
        for(int rr=0; rr<pnsd*pnen; ++rr)
        {
          for(int mm=0; mm<pnsd*pnen; ++mm)
          {
            emat(vi,ui)
              -=mat_w_q_o_n(vi,rr)*inv_s_q(rr,mm)*(mat_s_gradphi(mm,ui)+mat_s_o_n_phi(mm,ui));
          }
        }
      }
    }

    for (int vi=0; vi<pnen; ++vi)
    {
      for(int rr=0; rr<pnsd*pnen; ++rr)
      {
        for(int mm=0; mm<pnsd*pnen; ++mm)
        {
          erhs(vi)-= mat_w_q_o_n(vi,rr)*inv_s_q(rr,mm)*(-vec_s_o_n_phi_minus_g(mm)-vec_s_gradphi(mm));
        }
      }
    }
  }

  return;
}

#endif // CCADISCRET
#endif // D_FLUID3
