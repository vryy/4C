/*----------------------------------------------------------------------*/
/*!
\file condif3_impl.cpp

\brief Internal implementation of Condif3 element

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "condif3_impl.H"
#include "condif3_utils.H"
#include "../drt_mat/convecdiffus.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Epetra_SerialDenseSolver.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Condif3ImplInterface* DRT::ELEMENTS::Condif3ImplInterface::Impl(DRT::ELEMENTS::Condif3* c3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = c3->NumDofPerNode(*(c3->Nodes()[0]));
  int numscal = numdofpernode;
  if (DRT::Problem::Instance()->ProblemType() == "elch")
    numscal -= 1;

  switch (c3->Shape())
  {
  case DRT::Element::hex8:
  {
    static Condif3Impl<DRT::Element::hex8>* ch8;
    if (ch8==NULL)
      ch8 = new Condif3Impl<DRT::Element::hex8>(numdofpernode,numscal);
    return ch8;
  }
  case DRT::Element::hex20:
  {
    static Condif3Impl<DRT::Element::hex20>* ch20;
    if (ch20==NULL)
      ch20 = new Condif3Impl<DRT::Element::hex20>(numdofpernode,numscal);
    return ch20;
  }
  case DRT::Element::hex27:
  {
    static Condif3Impl<DRT::Element::hex27>* ch27;
    if (ch27==NULL)
      ch27 = new Condif3Impl<DRT::Element::hex27>(numdofpernode,numscal);
    return ch27;
  }
  case DRT::Element::tet4:
  {
    static Condif3Impl<DRT::Element::tet4>* ct4;
    if (ct4==NULL)
      ct4 = new Condif3Impl<DRT::Element::tet4>(numdofpernode,numscal);
    return ct4;
  }
  case DRT::Element::tet10:
  {
    static Condif3Impl<DRT::Element::tet10>* ct10;
    if (ct10==NULL)
      ct10 = new Condif3Impl<DRT::Element::tet10>(numdofpernode,numscal);
    return ct10;
  }
  case DRT::Element::wedge6:
  {
    static Condif3Impl<DRT::Element::wedge6>* cw6;
    if (cw6==NULL)
      cw6 = new Condif3Impl<DRT::Element::wedge6>(numdofpernode,numscal);
    return cw6;
  }
  case DRT::Element::wedge15:
  {
    static Condif3Impl<DRT::Element::wedge15>* cw15;
    if (cw15==NULL)
      cw15 = new Condif3Impl<DRT::Element::wedge15>(numdofpernode,numscal);
    return cw15;
  }
  case DRT::Element::pyramid5:
  {
    static Condif3Impl<DRT::Element::pyramid5>* cp5;
    if (cp5==NULL)
      cp5 = new Condif3Impl<DRT::Element::pyramid5>(numdofpernode,numscal);
    return cp5;
  }

  default:
    dserror("shape %d (%d nodes) not supported", c3->Shape(), c3->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Condif3Impl<distype>::Condif3Impl(int numdofpernode, int numscal)
  : numdofpernode_(numdofpernode),
    numscal_(numscal),
    xyze_(),
    bodyforce_(numdofpernode_),
    diffus_(numscal_),
    valence_(numscal_),
    shcacp_(0.0),
    funct_(),
    densfunct_(),
    deriv_(),
    deriv2_(),
    xjm_(),
    xij_(),
    derxy_(),
    derxy2_(),
    rhs_(numdofpernode_),
    hist_(numdofpernode_),
    velint_(),
    migvelint_(),
    tau_(numscal_),
    kart_(numscal_),
    xder2_(),
    fac_(0.0),
    conv_(),
    diff_(true),
    mig_(),
    gradpot_(),
    conint_(numscal_)
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Condif3Impl<distype>::Evaluate(
    Condif3*                   ele,
    ParameterList&             params,
    DRT::Discretization&       discretization,
    vector<int>&               lm,
    Epetra_SerialDenseMatrix&  elemat1_epetra,
    Epetra_SerialDenseMatrix&  elemat2_epetra,
    Epetra_SerialDenseVector&  elevec1_epetra,
    Epetra_SerialDenseVector&  elevec2_epetra,
    Epetra_SerialDenseVector&  elevec3_epetra,
    RefCountPtr<MAT::Material> mat,
    MATERIAL*                  actmat)
    {
  const string action = params.get<string>("action","none");
  if (action=="calc_condif_systemmat_and_residual")
  {
    // need current history vector and density vector
    RefCountPtr<const Epetra_Vector> hist = discretization.GetState("hist");
    RefCountPtr<const Epetra_Vector> densnp = discretization.GetState("densnp");
    if (hist==null || densnp==null)
      dserror("Cannot get state vector 'hist' and/or 'densnp'");

    // extract local values from the global vector
    vector<double> myhist(lm.size());
    vector<double> mydensnp(lm.size());
    DRT::UTILS::ExtractMyValues(*hist,myhist,lm);
    DRT::UTILS::ExtractMyValues(*densnp,mydensnp,lm);

    // get control parameter
    const bool is_stationary = params.get<bool>("using stationary formulation");
    const double time = params.get<double>("total time");

    // One-step-Theta: timefac = theta*dt
    // BDF2:           timefac = 2/3 * dt
    double timefac = 0.0;
    if (not is_stationary)
    {
      timefac = params.get<double>("thsl");
      if (timefac < 0.0) dserror("thsl is negative.");
    }

    // get (weighted) velocity at the nodes
    // compare also with DRT::UTILS::ExtractMyValues()
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    //const int iel = ele->NumNode();
    const int nsd=3;
    Epetra_SerialDenseVector evel(nsd*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity);

    // get flag for fine-scale subgrid diffusivity
    string fssgd = params.get<string>("fs subgrid diffusivity","No");

    // set flag for type of scalar whether it is temperature or not
    string scaltypestr=params.get<string>("problem type");
    bool temperature = false;
    if(scaltypestr =="loma") temperature = true;

    // paramters needed for ELCH ;-)
    vector<double> myphinp(lm.size());
    double frt(0.0);
    if(scaltypestr =="elch")
    {
      RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phinp==null) dserror("Cannot get state vector 'phinp'");
      // extract local values from the global vector
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // create objects for element arrays
    vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> > ephinp(numscal_);
    vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> > ehist(numdofpernode_);
    LINALG::FixedSizeSerialDenseMatrix<iel,1> edensnp;
    LINALG::FixedSizeSerialDenseMatrix<3,iel> evelnp;
    LINALG::FixedSizeSerialDenseMatrix<iel,1> epotnp;

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephinp[k](i,0) = myphinp[k+(i*numdofpernode_)];
      }
      for (int k = 0; k< numdofpernode_; ++k)
      {
        // the history vectors contains information of time step t_n
        ehist[k](i,0) = myhist[k+(i*numdofpernode_)];
      }

      // split for each tranported scalar, insert into element arrays
      evelnp(0,i) = evel[   i*3 ];
      evelnp(1,i) = evel[1+(i*3)];
      evelnp(2,i) = evel[2+(i*3)];

      // insert density vector into element array 
      // (only take values belonging to the first transported scalar!)
      edensnp(i,0) = mydensnp[0+(i*numdofpernode_)];

      if(scaltypestr =="elch")
      {
        // get values for el. potential at element nodes
        epotnp(i) = myphinp[i*numdofpernode_+numscal_];
      }
      else
        epotnp(i) = 0.0;
    } // for i

    // calculate element coefficient matrix and rhs
    Sysmat(
        ele,
        ephinp,
        ehist,
        edensnp,
        epotnp,
        elemat1_epetra,
        elevec1_epetra,
        elevec2_epetra,
        actmat,
        time,
        timefac,
        evelnp,
        temperature,
        fssgd,
        is_stationary,
        frt);
  }
  else if (action =="initialize_one_step_theta")
    // calculate time derivative for time value t_0
  {
    const double time = params.get<double>("total time");
    const double timefac = params.get<double>("thsl");

    // need initial field
    RefCountPtr<const Epetra_Vector> phi0 = discretization.GetState("phi0");
    RefCountPtr<const Epetra_Vector> dens0 = discretization.GetState("dens0");
    if (phi0==null || dens0==null)
      dserror("Cannot get state vector 'phi0' and/or 'densnp'");

    // extract local values from the global vector
    vector<double> myphi0(lm.size());
    vector<double> mydens0(lm.size());
    DRT::UTILS::ExtractMyValues(*phi0,myphi0,lm);
    DRT::UTILS::ExtractMyValues(*dens0,mydens0,lm);

    // get initial velocity values at the nodes
    // compare also with DRT::UTILS::ExtractMyValues()
    const RCP<Epetra_MultiVector> velocity = params.get< RCP<Epetra_MultiVector> >("velocity field",null);
    const int nsd=3;
    Epetra_SerialDenseVector evel(nsd*iel);
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evel,velocity);

    // get flag for fine-scale subgrid diffusivity
    string fssgd = params.get<string>("fs subgrid diffusivity","No");

    // set flag for type of scalar whether it is temperature or not
    string scaltypestr=params.get<string>("problem type");
    bool temperature = false;
    if(scaltypestr =="loma") temperature = true;


    // create objects for element arrays
    vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> > ephi0(numscal_);
    LINALG::FixedSizeSerialDenseMatrix<iel,1> edens0;
    LINALG::FixedSizeSerialDenseMatrix<3,iel> evel0;
    LINALG::FixedSizeSerialDenseMatrix<iel,1> epot0;

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< numscal_; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephi0[k](i,0) = myphi0[k+(i*numdofpernode_)];
      }

      // split for each tranported scalar, insert into element arrays
      evel0(0,i) = evel[   i*3 ];
      evel0(1,i) = evel[1+(i*3)];
      evel0(2,i) = evel[2+(i*3)];

      // insert density vector into element array 
      // (only take values belonging to the first transported scalar!)
      edens0(i,0) = mydens0[0+(i*numdofpernode_)];

      if(scaltypestr =="elch")
      {
        // get values for el. potential at element nodes
        epot0(i) = myphi0[i*numdofpernode_+numscal_];
      }
      else
        epot0(i) = 0.0;
    } // for i

    double frt(0.0);
    if(scaltypestr =="elch")
    {
      // get parameter F/RT
      frt = params.get<double>("frt");
    }

    // calculate mass matrix and rhs
    InitializeOST(
        ele,
        ephi0,
        edens0,
        epot0,
        elemat1_epetra,
        elevec1_epetra,
        elevec2_epetra,
        actmat,
        time,
        timefac,
        evel0,
        temperature,
        fssgd,
        frt);
  }
  else if (action =="calc_subgrid_diffusivity_matrix")
  // calculate normalized subgrid-diffusivity matrix
  {
    // get control parameter
    const bool is_stationary = params.get<bool>("using stationary formulation");

    // One-step-Theta: timefac = theta*dt
    // BDF2:           timefac = 2/3 * dt
    double timefac = 0.0;
    if (not is_stationary)
    {
      timefac = params.get<double>("thsl");
      if (timefac < 0.0) dserror("No thsl supplied");
    }

    // calculate mass matrix and rhs
    CalcSubgridDiffMatrix(
        ele,
        elemat1_epetra,
        timefac,
        is_stationary);
  }
  else if (action=="calc_elch_kwok_error")
  {
    // need current solution
    RefCountPtr<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

    // create objects for element arrays
    LINALG::FixedSizeSerialDenseMatrix<iel,2> ephinp;
    LINALG::FixedSizeSerialDenseMatrix<iel,1> epotnp;

    // fill element arrays
    for (int i=0;i<iel;++i)
    {
      for (int k = 0; k< 2; ++k)
      {
        // split for each tranported scalar, insert into element arrays
        ephinp(i,k) = myphinp[k+(i*numdofpernode_)];
      }

      // get values for el. potential at element nodes
      epotnp(i) = myphinp[i*numdofpernode_+numscal_];
    } // for i

    CalErrorComparedToAnalytSolution(
        ele,
        params,
        ephinp,
        epotnp,
        elevec1_epetra,
        actmat);
  }
  else
    dserror("Unknown type of action for Condif3: %s",action.c_str());

return 0;
}


/*----------------------------------------------------------------------*
 |  calculate system matrix and rhs (public)                 g.bau 08/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::Sysmat(
    const DRT::ELEMENTS::Condif3*   ele, ///< the element those matrix is calculated
    const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ephinp,///< current scalar field
    const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ehist, ///< rhs from beginning of time step
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>& edens, ///< density*shc
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>& epotnp, ///< el. potential at element nodes
    Epetra_SerialDenseMatrix&       sys_mat,///< element matrix to calculate
    Epetra_SerialDenseVector&       residual, ///< element rhs to calculate
    Epetra_SerialDenseVector&       subgrdiff, ///< subgrid-diff.-scaling vector
    const struct _MATERIAL*         material, ///< material pointer
    const double                    time, ///< current simulation time
    const double                    timefac, ///< time discretization factor
    const LINALG::FixedSizeSerialDenseMatrix<3,iel>& evelnp,///< nodal velocities at t_{n+1}
    const bool                      temperature, ///< temperature flag
    const string                    fssgd, ///< subgrid-diff. flag
    const bool                      is_stationary, ///< flag indicating stationary formulation
    const double                    frt ///< factor F/RT needed for ELCH calculations
)
{
  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze_(0,i)=ele->Nodes()[i]->X()[0];
    xyze_(1,i)=ele->Nodes()[i]->X()[1];
    xyze_(2,i)=ele->Nodes()[i]->X()[2];
  }

  // dead load in element nodes
  BodyForce(ele,time);

  // get material constants
  GetMaterialParams(material,temperature);

  /*----------------------------------------------------------------------*/
  // calculation of stabilization parameter(s) tau
  /*----------------------------------------------------------------------*/
  CalTau(ele,subgrdiff,evelnp,epotnp,timefac,fssgd,is_stationary,false,frt);

  /*----------------------------------------------------------------------*/
  // integration loop for one condif3 element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool higher_order_ele = SCATRA::is3DHigherOrderElement<distype>();

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(SCATRA::get3DOptimalGaussrule<distype>());

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,higher_order_ele,ele);

    // density*specific heat capacity-weighted shape functions
    densfunct_.EMultiply(funct_,edens);

    // get (density*specific heat capacity-weighted) velocity at integration point
    velint_.Multiply(evelnp,funct_);

    //------------ get values of variables at integration point
    for (int k = 0;k<numdofpernode_;++k)     // loop of each transported sclar
    {
      // get history data at integration point (weighted by density)
      hist_[k] = densfunct_.Dot(ehist[k]);

      // get bodyforce in gausspoint (divided by shcacp for temperature eq.)
      rhs_[k] = bodyforce_[k].Dot(funct_) / shcacp_;
    }

    //----------- perform integration for entire matrix and rhs
    if (numdofpernode_-numscal_== 0) // 'standard' scalar transport
    {
      for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
      {
        if (not is_stationary)
          CalMat(sys_mat,residual,higher_order_ele,timefac,k);
        else
          CalMatStationary(sys_mat,residual,higher_order_ele,k);
      } // loop over each scalar
    }
    else  // ELCH problems
     CalMatElch(sys_mat,residual,ephinp,epotnp,higher_order_ele,frt,is_stationary,timefac);

  } // integration loop

  return;
}


/*----------------------------------------------------------------------*
 |  get the body force  (private)                              gjb 06/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::BodyForce(
    const DRT::ELEMENTS::Condif3* ele, 
    const double time
)
{
  vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique VolumeNeumann condition
  DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);

  if (myneumcond.size()>1)
    dserror("more than one VolumeNeumann cond on one node");

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
        curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
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

    // set this condition to the bodyforce array
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      for (int jnode=0; jnode<iel; jnode++)
      {
        (bodyforce_[idof])(jnode) = (*onoff)[idof]*(*val)[idof]*curvefac;
      }
    }
  }
  else
  {
    for(int idof=0;idof<numdofpernode_;idof++)
    {
      // we have no dead load
      bodyforce_[idof].Clear();
    }
  }

  return;

} //Condif3Impl::BodyForce


/*----------------------------------------------------------------------*
 |  get the material constants  (private)                      gjb 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::GetMaterialParams
(
    const struct _MATERIAL*   material,
    const bool&               temperature  
)
{
// get diffusivity / diffusivities
if (material->mattyp == m_matlist)
{
  for (int k = 0;k<numscal_;++k)
  {
    const int matid = material->m.matlist->matids[k];
    const _MATERIAL& singlemat =  DRT::Problem::Instance()->Material(matid-1);

    if (singlemat.mattyp == m_ion)
    {
      valence_[k]= singlemat.m.ion->valence;
      diffus_[k]= singlemat.m.ion->diffusivity;
    }
    else if (singlemat.mattyp == m_condif)
      diffus_[k]= singlemat.m.condif->diffusivity;
    else
      dserror("material type is not allowed");
  }
  // set specific heat capacity at constant pressure to 1.0
  shcacp_ = 1.0;
}
else if (material->mattyp == m_condif)
{
  dsassert(numdofpernode_==1,"more than 1 dof per node for condif material");

  // in case of a temperature equation, we get thermal conductivity instead of
  // diffusivity and have to divide by the specific heat capacity at constant
  // pressure; otherwise, it is the "usual" diffusivity
  if (temperature)
  {
    shcacp_ = material->m.condif->shc;
    diffus_[0] = material->m.condif->diffusivity/shcacp_;
  }
  else
  {
    // set specific heat capacity at constant pressure to 1.0, get diffusivity
    shcacp_ = 1.0;
    diffus_[0] = material->m.condif->diffusivity;
  }
}
else
  dserror("Material type is not supported");

return;
} //Condif3Impl::GetMaterialParams


/*----------------------------------------------------------------------*
 |  calculate stabilization parameter  (private)              gjb 06/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalTau(
    const DRT::ELEMENTS::Condif3*&          ele,
    Epetra_SerialDenseVector&               subgrdiff,
    const LINALG::FixedSizeSerialDenseMatrix<3,iel>& evelnp,
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>&         epot,
    const double&                           timefac,
    string                                  fssgd,
    const bool&                             is_stationary,
    const bool                              initial,
    const double&                           frt
  )
{
  /*------------------------------------------------------- initialize ---*/
  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D intrule_stabili = SCATRA::getIntegrationRuleForStabilization<distype>();

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D  intpoints_tau = getIntegrationPoints3D(intrule_stabili);

  // prepare the standard FE stuff for this single integration point
  // we do not need second derivatives for the calculation of tau
  // EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau,0,distype,false,ele);

  // shape functions and derivs at element center
  const double e1    = intpoints_tau.qxg[0][0];
  const double e2    = intpoints_tau.qxg[0][1];
  const double e3    = intpoints_tau.qxg[0][2];
  const double wquad = intpoints_tau.qwgt[0];

  // shape functions and their derivatives
  DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 0.0)
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
  if (abs(det) < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", ele->Id(), det);

  const double vol = wquad*det;

  // There exist different definitions for 'the' characteristic element length hk:
  // 1) get element length for tau_Mp/tau_C: volume-equival. diameter
  // const double hk = pow((6.*vol/PI),(1.0/3.0));

  // 2) streamlength (based on velocity vector at element centre)
  if (numdofpernode_-numscal_== 1) // ELCH
  {
    // compute global derivatives
    derxy_.Multiply(xij_,deriv_);

    // get "migration velocity" divided by D_k*z_k at element center
    migvelint_.Multiply(-frt,derxy_,epot);
    
  } // if ELCH

  
/*  double strle = 0.0;
  if (vel_norm>1e-6)
  {
    double val = 0;
    for (int i=0;i<iel;++i)
    {
      double sum = 0;
      for (int j=0;j<3;++j)
      {
        sum += velint_[j]*derxy_(j,i);
      }
      val+= abs(sum);
    }
    strle = 2.0*vel_norm/val; //this formula is not working in 3D in case of HEX8 elements!!
  }
  else
  //case: 'zero' velocity vector => tau will be very small in diffusion-dominated regions
  // => usage of arbitrary vector velint = (1 0 0)^T in above formula possible
  {
     double val = 0;
     for (int i=0;i<iel;++i)
     {
       val+=abs(derxy_(0,i));
     }
     strle = 2.0/val;
   }
   //const double hk = strle;*/

  // 3) use cubic root of the element volume as characteristic length
  const double hk = pow(vol,(1.0/3.0));

  // get element type constant for tau
  const double mk = SCATRA::MK<distype>();

  // get (density*specific heat capacity-weighted) velocity at element center
  velint_.Multiply(evelnp,funct_);

  // some necessary parameter definitions
  double vel_norm, epe1, epe2, xi1, xi2;

  for(int k = 0;k<numscal_;++k) // loop over all transported scalars
  {
    if (numdofpernode_ - numscal_ == 1) // ELCH
    {
      const double Dkzk = diffus_[k]*valence_[k];
      // get Euclidean norm of effective velocity at element center:
      // (weighted) convective velocity + individual migration velocity
      vel_norm = sqrt(DSQR(velint_(0)+Dkzk*migvelint_(0)) + DSQR(velint_(1)+Dkzk*migvelint_(1)) + DSQR(velint_(2)+Dkzk*migvelint_(2)));
    }
    else
    {
      // get Euclidean norm of (weighted) velocity at element center
      vel_norm = velint_.Norm2();
    }

  // stabilization parameter definition according to Franca and Valentin (2000)
  if (is_stationary == false)
  {
      /* parameter relating diffusive : reactive forces */
      epe1 = 2.0 * timefac * diffus_[k] / (mk * DSQR(hk)); 
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");
      /* parameter relating convective : diffusive forces */
      epe2 = mk * vel_norm * hk / diffus_[k];
      xi1 = DMAX(epe1,1.0);
      xi2 = DMAX(epe2,1.0);

      tau_[k] = DSQR(hk)/((DSQR(hk)*xi1)/timefac + (2.0*diffus_[k]/mk)*xi2);
  }
  else
  {
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");
      /* parameter relating convective : diffusive forces */
      epe2 = mk * vel_norm * hk / diffus_[k];
      xi2 = DMAX(epe2,1.0);

      tau_[k] = (DSQR(hk)*mk)/(2.0*diffus_[k]*xi2);
  }

  // compute artificial diffusivity kappa_art_[k]
  if (fssgd == "artificial_all" and (not initial))
  {
      if (diffus_[k] == 0.0) dserror("diffusivity is zero: Preventing division by zero at evaluation of stabilization parameter");
      /* parameter relating convective : diffusive forces */
      epe2 = mk * vel_norm * hk / diffus_[k];
      xi2 = DMAX(epe2,1.0);

      kart_[k] = (DSQR(hk)*mk*DSQR(vel_norm))/(2.0*diffus_[k]*xi2);

      for (int vi=0; vi<iel; ++vi)
      {
        subgrdiff(vi) = kart_[k]/ele->Nodes()[vi]->NumElement();
      }
  }
  else if (fssgd == "artificial_small" || fssgd == "Smagorinsky_all" ||
           fssgd == "Smagorinsky_small")
    dserror("only all-scale artficial diffusivity for convection-diffusion problems possible so far!\n");

  } // loop over scalars

#if 0
   // some debug output
   cout<<"hk (volume equiv diam)            = "<<pow((6.*vol/PI),(1.0/3.0))<<endl;
   cout<<"strle                             = "<<strle<<endl;
   cout<<"hk (cubic root of element volume) = "<<pow(vol,(1.0/3.0))<<endl;
   cout<<"tau = "<<tau<<endl;
#endif

  return;
} //Condif3Impl::Caltau


/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     gjb 08/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntegrationPoints3D& intpoints, ///< integration points
    const int&                             iquad,     ///< id of current Gauss point
    const bool&                            higher_order_ele,///< are second derivatives needed?
    const DRT::ELEMENTS::Condif3*          ele        ///< the element
)
{
  // coordinates of the current integration point
  const double e1 = intpoints.qxg[iquad][0];
  const double e2 = intpoints.qxg[iquad][1];
  const double e3 = intpoints.qxg[iquad][2];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function_3D(funct_,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv_,e1,e2,e3,distype);

  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  xjm_.MultiplyNT(deriv_,xyze_);
  const double det = xij_.Invert(xjm_);

  if (det < 0.0)
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
  if (abs(det) < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO JACOBIAN DETERMINANT: %f", ele->Id(), det);

  // set integration factor: fac = Gauss weight * det(J)
  fac_ = intpoints.qwgt[iquad]*det;

  // compute global derivatives
  derxy_.Multiply(xij_,deriv_);

  // compute second global derivatives (if needed)
  if (higher_order_ele) 
    CalSecondDeriv(e1,e2,e3);
  else 
    derxy2_.Clear();

  // say goodbye
  return;
}


/*----------------------------------------------------------------------*
 |  calculate second global derivatives w.r.t. x,y,z at point r,s,t
 |                                            (private)      gammi 07/07
 |
 | From the six equations
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dr^2     dr | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ------ = -- | --*-- + --*-- + --*-- |
 |  ds^2     ds | ds dx   ds dy   ds dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 |  ----   = -- | --*-- + --*-- + --*-- |
 |  dt^2     dt | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dr     ds | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | dt dr     dt | dr dx   dr dy   dr dz |
 |              +-                     -+
 |
 |              +-                     -+
 |  d^2N     d  | dx dN   dy dN   dz dN |
 | -----   = -- | --*-- + --*-- + --*-- |
 | ds dt     ds | dt dx   dt dy   dt dz |
 |              +-                     -+
 |
 | the matrix (jacobian-bar matrix) system
 |
 | +-                                                                                         -+   +-    -+
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dr/          \dr/           \dr/             dr dr           dr dr           dr dr     |   | dx^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \ds/          \ds/           \ds/             ds ds           ds ds           ds ds     |   | dy^2 |
 | |                                                                                           |   |      |
 | |   /dx\^2        /dy\^2         /dz\^2           dy dx           dz dx           dy dz     |   | d^2N |
 | |  | -- |        | ---|         | ---|          2*--*--         2*--*--         2*--*--     |   | ---- |
 | |   \dt/          \dt/           \dt/             dt dt           dt dt           dt dt     |   | dz^2 |
 | |                                                                                           | * |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr ds         dr ds          dr ds        dr ds   ds dr   dr ds   ds dr  dr ds   ds dr  |   | dxdy |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dr dt         dr dt          dr dt        dr dt   dt dr   dr dt   dt dr  dr dt   dt dr  |   | dxdz |
 | |                                                                                           |   |      |
 | |   dx dx         dy dy          dz dz        dx dy   dx dy   dx dz   dx dz  dy dz   dy dz  |   | d^2N |
 | |   --*--         --*--          --*--        --*-- + --*--   --*-- + --*--  --*-- + --*--  |   | ---- |
 | |   dt ds         dt ds          dt ds        dt ds   ds dt   dt ds   ds dt  dt ds   ds dt  |   | dydz |
 | +-                                                                                         -+   +-    -+
 |
 |                  +-    -+     +-                           -+
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dr^2 |     | dr^2 dx   dr^2 dy   dr^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | ds^2 |     | ds^2 dx   ds^2 dy   ds^2 dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dt^2 |     | dt^2 dx   dt^2 dy   dt^2 dz |
 |              =   |      |  -  |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drds |     | drds dx   drds dy   drds dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2y dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | drdt |     | drdt dx   drdt dy   drdt dz |
 |                  |      |     |                             |
 |                  | d^2N |     | d^2x dN   d^2y dN   d^2z dN |
 |                  | ---- |     | ----*-- + ----*-- + ----*-- |
 |                  | dtds |     | dtds dx   dtds dy   dtds dz |
 |                  +-    -+     +-                           -+
 |
 |
 | is derived. This is solved for the unknown global derivatives.
 |
 |
 |             jacobian_bar * derxy2 = deriv2 - xder2 * derxy
 |                                              |           |
 |                                              +-----------+
 |                                              'chainrulerhs'
 |                                     |                    |
 |                                     +--------------------+
 |                                          'chainrulerhs'
 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalSecondDeriv(
    const double&                            e1,      ///< first coordinate of GP
    const double&                            e2,      ///< second coordinate of GP
    const double&                            e3       ///< third coordinate of GP
    )
{
  /*--- get the second derivatives of standard element at current GP */
  DRT::UTILS::shape_function_3D_deriv2(deriv2_,e1,e2,e3,distype);
dserror("Please check CalSecondDeriv first");

  /*----------- now we have to compute the second global derivatives */
  static LINALG::FixedSizeSerialDenseMatrix<6,6> bm;

  /*------------------------------------------------- initialization */
  derxy2_.Clear();

  // calculate elements of jacobian_bar matrix
  bm(0,0) = xjm_(0,0)*xjm_(0,0);
  bm(1,0) = xjm_(1,0)*xjm_(1,0);
  bm(2,0) = xjm_(2,0)*xjm_(2,0);
  bm(3,0) = xjm_(0,0)*xjm_(1,0);
  bm(4,0) = xjm_(0,0)*xjm_(2,0);
  bm(5,0) = xjm_(2,0)*xjm_(1,0);

  bm(0,1) = xjm_(0,1)*xjm_(0,1);
  bm(1,1) = xjm_(1,1)*xjm_(1,1);
  bm(2,1) = xjm_(2,1)*xjm_(2,1);
  bm(3,1) = xjm_(0,1)*xjm_(1,1);
  bm(4,1) = xjm_(0,1)*xjm_(2,1);
  bm(5,1) = xjm_(2,1)*xjm_(1,1);

  bm(0,2) = xjm_(0,2)*xjm_(0,2);
  bm(1,2) = xjm_(1,2)*xjm_(1,2);
  bm(2,2) = xjm_(2,2)*xjm_(2,2);
  bm(3,2) = xjm_(0,2)*xjm_(1,2);
  bm(4,2) = xjm_(0,2)*xjm_(2,2);
  bm(5,2) = xjm_(2,2)*xjm_(1,2);

  bm(0,3) = 2.*xjm_(0,0)*xjm_(0,1);
  bm(1,3) = 2.*xjm_(1,0)*xjm_(1,1);
  bm(2,3) = 2.*xjm_(2,0)*xjm_(2,1);
  bm(3,3) = xjm_(0,0)*xjm_(1,1)+xjm_(1,0)*xjm_(0,1);
  bm(4,3) = xjm_(0,0)*xjm_(2,1)+xjm_(2,0)*xjm_(0,1);
  bm(5,3) = xjm_(1,0)*xjm_(2,1)+xjm_(2,0)*xjm_(1,1);

  bm(0,4) = 2.*xjm_(0,0)*xjm_(0,2);
  bm(1,4) = 2.*xjm_(1,0)*xjm_(1,2);
  bm(2,4) = 2.*xjm_(2,0)*xjm_(2,2);
  bm(3,4) = xjm_(0,0)*xjm_(1,2)+xjm_(1,0)*xjm_(0,2);
  bm(4,4) = xjm_(0,0)*xjm_(2,2)+xjm_(2,0)*xjm_(0,2);
  bm(5,4) = xjm_(1,0)*xjm_(2,2)+xjm_(2,0)*xjm_(1,2);

  bm(0,5) = 2.*xjm_(0,1)*xjm_(0,2);
  bm(1,5) = 2.*xjm_(1,1)*xjm_(1,2);
  bm(2,5) = 2.*xjm_(2,1)*xjm_(2,2);
  bm(3,5) = xjm_(0,1)*xjm_(1,2)+xjm_(1,1)*xjm_(0,2);
  bm(4,5) = xjm_(0,1)*xjm_(2,2)+xjm_(2,1)*xjm_(0,2);
  bm(5,5) = xjm_(1,1)*xjm_(2,2)+xjm_(2,1)*xjm_(1,2);

  /*------------------ determine 2nd derivatives of coord.-functions */

  /*
  |
  |         0 1 2              0...iel-1
  |        +-+-+-+             +-+-+-+-+        0 1 2
  |        | | | | 0           | | | | | 0     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | 0
  |        | | | | 1           | | | | | 1   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 2           | | | | | 2     +-+-+-+
  |        +-+-+-+       =     +-+-+-+-+       | | | | .
  |        | | | | 3           | | | | | 3     +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 4           | | | | | 4   * +-+-+-+ .
  |        +-+-+-+             +-+-+-+-+       | | | | .
  |        | | | | 5           | | | | | 5     +-+-+-+
  |        +-+-+-+             +-+-+-+-+       | | | | iel-1
  |                                            +-+-+-+
  |
  |        xder2               deriv2          xyze^T
  |
  |
  |                                     +-                  -+
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dr^2   dr^2   dr^2 |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | ds^2   ds^2   ds^2 |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dt^2   dt^2   dt^2 |
  |               yields    xder2  =    |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drds   drds   drds |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | drdt   drdt   drdt |
  |                                     |                    |
  |                                     | d^2x   d^2y   d^2z |
  |                                     | ----   ----   ---- |
  |                                     | dsdt   dsdt   dsdt |
  |                                     +-                  -+
  |
  |
  */

  //xder2_ = blitz::sum(deriv2_(i,k)*xyze_(j,k),k);
 /* for (int i = 0; i < 6; ++i)
  {
      for (int j = 0; j < 3; ++j)
      {
          for (int k = 0; k < iel; ++k)
          {
              xder2_(i,j) += deriv2_(i,k)*xyze_(j,k);
          }
      }
  }*/
  xder2_.MultiplyNT(deriv2_,xyze_);

  /*
  |        0...iel-1             0 1 2
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 0          | | | | 0
  |        +-+-+-+-+            +-+-+-+            0...iel-1
  |        | | | | | 1          | | | | 1         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 0
  |        | | | | | 2          | | | | 2         +-+-+-+-+
  |        +-+-+-+-+       =    +-+-+-+       *   | | | | | 1 * (-1)
  |        | | | | | 3          | | | | 3         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+           | | | | | 2
  |        | | | | | 4          | | | | 4         +-+-+-+-+
  |        +-+-+-+-+            +-+-+-+
  |        | | | | | 5          | | | | 5          derxy
  |        +-+-+-+-+            +-+-+-+
  |
  |       chainrulerhs          xder2
  */

  //derxy2_ = -blitz::sum(xder2_(i,k)*derxy_(k,j),k);
  //derxy2_ = deriv2 - blitz::sum(xder2(i,k)*derxy(k,j),k);
  for (int i = 0; i < 6; ++i)
  {
      for (int j = 0; j < iel; ++j)
      {
          derxy2_(i,j) += deriv2_(i,j);
          for (int k = 0; k < 3; ++k)
          {
              derxy2_(i,j) -= xder2_(i,k)*derxy_(k,j);
          }
      }
  }

  /*
  |        0...iel-1            0...iel-1         0...iel-1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 0          | | | | | 0       | | | | | 0
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 1          | | | | | 1       | | | | | 1
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 2          | | | | | 2       | | | | | 2
  |        +-+-+-+-+       =    +-+-+-+-+    +    +-+-+-+-+
  |        | | | | | 3          | | | | | 3       | | | | | 3
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 4          | | | | | 4       | | | | | 4
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |        | | | | | 5          | | | | | 5       | | | | | 5
  |        +-+-+-+-+            +-+-+-+-+         +-+-+-+-+
  |
  |       chainrulerhs         chainrulerhs        deriv2
  */

  //derxy2_ += deriv2_;

  /* make LR decomposition and solve system for all right hand sides
   * (i.e. the components of chainrulerhs)
  |
  |          0  1  2  3  4  5         i        i
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 0     | | 0    | | 0
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 1     | | 1    | | 1
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 2     | | 2    | | 2
  |        +--+--+--+--+--+--+    *  +-+   =  +-+      for i=0...iel-1
  |        |  |  |  |  |  |  | 3     | | 3    | | 3
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 4     | | 4    | | 4
  |        +--+--+--+--+--+--+       +-+      +-+
  |        |  |  |  |  |  |  | 5     | | 5    | | 5
  |        +--+--+--+--+--+--+       +-+      +-+
  |                                   |        |
  |                                   |        |
  |                                   derxy2[i]|
  |                                    |
  |                                    chainrulerhs[i]
  |
  |   yields
  |
  |                      0...iel-1
  |                      +-+-+-+-+
  |                      | | | | | 0 = drdr
  |                      +-+-+-+-+
  |                      | | | | | 1 = dsds
  |                      +-+-+-+-+
  |                      | | | | | 2 = dtdt
  |            derxy2 =  +-+-+-+-+
  |                      | | | | | 3 = drds
  |                      +-+-+-+-+
  |                      | | | | | 4 = drdt
  |                      +-+-+-+-+
  |                      | | | | | 5 = dsdt
  |                  +-+-+-+-+
  */

  LINALG::FixedSizeSerialDenseSolver<6,6,iel> solver;
  solver.SetMatrix(bm);

  // No need for a separate rhs. We assemble the rhs to the solution
  // vector. The solver will destroy the rhs and return the solution.
  solver.SetVectors(derxy2_,derxy2_);
  solver.Solve();

  return;
} //Condif3Impl::CalSecondDeriv


/*----------------------------------------------------------------------*
 |  evaluate instationary convection-diffusion matrix (private)gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif3 element are calculated for the instationary
case. The procedure is based on the Rothe method of first discretizing in
time. Hence the resulting terms include coefficients containing time
integration variables such as theta or delta t which are represented by
'timefac'.

The stabilization is based on the residuum:

R = rho * c_p * phi + timefac * rho * c_p * u * grad(phi)
                    - timefac * diffus * laplace(phi) - rhsint

The corresponding weighting operators are
L = timefac * rho * c_p * u * grad(w) +/- timefac * diffus * laplace(w)

'+': USFEM (default)
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\return void
------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalMat(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    const bool                higher_order_ele,
    const double&             timefac,
    const int&                dofindex
    )
{
static double             rhsint;           /* rhs at int. point     */

// stabilization parameter
const double taufac = tau_[dofindex]*fac_;

// integration factors and coefficients of single terms
const double timefacfac  = timefac * fac_;
const double timetaufac  = timefac * taufac;

/*-------------------------------- evaluate rhs at integration point ---*/
rhsint = hist_[dofindex] + rhs_[dofindex]*timefac;

/* convective part */
/* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
conv_.MultiplyTN(derxy_,velint_);


if (higher_order_ele)
{
  for (int i=0; i<iel; i++)
  {
    /* diffusive part */
    /* diffus * ( N,xx  +  N,yy +  N,zz ) */
    diff_(i) = diffus_[dofindex] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
  }
}

/*--------------------------------- now build single stiffness terms ---*/
const int numdof =numdofpernode_;
// -------------------------------------------System matrix
for (int vi=0; vi<iel; ++vi)
{
  for (int ui=0; ui<iel; ++ui)
  {
    /* Standard Galerkin terms: */
    /* transient term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += fac_*funct_(vi)*densfunct_(ui) ;

    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += timefacfac*funct_(vi)*conv_(ui) ;

    /* diffusive term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += timefacfac*diffus_[dofindex]*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi)) ;

    /* Stabilization terms: */
    /* 1) transient stabilization (USFEM assumed here, sign change necessary for GLS) */
    /* transient term */
    //estif(vi, ui) += -taufac*densfunct_(vi)*densfunct_[ui] ;

    /* convective term */
    //estif(vi, ui) += -timetaufac*densfunct_(vi)*conv_[ui] ;

    /* diffusive term */
    //if (higher_order_ele) estif(vi, ui) += timetaufac*densfunct_(vi)*diff[ui] ;

    /* 2) convective stabilization */
    /* transient term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*conv_(vi)*densfunct_(ui);

    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += timetaufac*conv_(vi)*conv_(ui) ;
  }
}

if (higher_order_ele)
{
  for (int vi=0; vi<iel; ++vi)
  {
    for (int ui=0; ui<iel; ++ui)
    {
      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += -timetaufac*conv_(vi)*diff_(ui) ;

      /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
      /* transient term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*diff_(vi)*densfunct_(ui) ;

      /* convective term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += timetaufac*diff_(vi)*conv_(ui) ;

      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += -timetaufac*diff_(vi)*diff_(ui) ;
    }
  }
}

// ----------------------------------------------RHS
for (int vi=0; vi<iel; ++vi)
{
  /* RHS source term */
  eforce[vi*numdof+dofindex] += fac_*funct_(vi)*rhsint ;

  /* transient stabilization of RHS source term */
  //eforce(vi) += -taufac*densfunct_(vi)*rhsint ;

  /* convective stabilization of RHS source term */
  eforce[vi*numdof+dofindex] += taufac*conv_(vi)*rhsint ;

  /* diffusive stabilization of RHS source term */
  if (higher_order_ele) eforce[vi*numdof+dofindex] += taufac*diff_(vi)*rhsint ;
}

return;
} //Condif3Impl::Condif3CalMat


/*----------------------------------------------------------------------*
 |  evaluate stationary convection-diffusion matrix (private)  gjb 06/08|
 *----------------------------------------------------------------------*/

/*
In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilized condif3 element are calculated for the stationary
case.

The stabilization is based on the residuum:

R = rho * c_p * u * grad(phi) - diffus *  laplace(phi) - rhsint

The corresponding weighting operators are
L = rho * c_p * u * grad(w) +/- diffus *  laplace(w)

'+': USFEM (default)
'-': GLS


The calculation proceeds as follows.
1) obtain single operators of R and L
2) build Galerkin terms from them
3) build stabilizing terms from them
4) build Galerkin and stabilizing terms of RHS

NOTE: Galerkin and stabilization matrices are calculated within one
      routine.


for further comments see comment lines within code.

</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\return void
------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalMatStationary(
    Epetra_SerialDenseMatrix& estif,
    Epetra_SerialDenseVector& eforce,
    const bool                higher_order_ele,
    const int&                dofindex
    )
{
static double rhsint;           /* rhs at int. point     */
const double  fac_diffus = fac_*diffus_[dofindex];

// stabilization parameter
const double taufac = tau_[dofindex]*fac_;

/*------------------------------------- set rhs at integration point ---*/
rhsint = rhs_[dofindex];

/* convective part */
/* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
conv_.MultiplyTN(derxy_,velint_);

if (higher_order_ele)
{
  for (int i=0; i<iel; i++)
  {
    /* diffusive part */
    /* diffus * ( N,xx  +  N,yy +  N,zz ) */
    diff_(i) = diffus_[dofindex] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
  }
}

/*--------------------------------- now build single stiffness terms ---*/
const int numdof = numdofpernode_;
// -------------------------------------------System matrix
for (int vi=0; vi<iel; ++vi)
{
  for (int ui=0; ui<iel; ++ui)
  {
    /* Standard Galerkin terms: */
    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += fac_*funct_(vi)*conv_(ui) ;

    /* diffusive term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += fac_diffus*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi)+ derxy_(2, ui)*derxy_(2, vi));

    /* Stabilization term: */
    /* 1) convective stabilization */

    /* convective term */
    estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*conv_(vi)*conv_(ui) ;
  }
}

if (higher_order_ele)
{
  for (int vi=0; vi<iel; ++vi)
  {
    for (int ui=0; ui<iel; ++ui)
    {
      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += -taufac*conv_(vi)*diff_(ui) ;

      /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

      /* convective term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) += taufac*diff_(vi)*conv_(ui) ;

      /* diffusive term */
      estif(vi*numdof+dofindex, ui*numdof+dofindex) -= taufac*diff_(vi)*diff_(ui) ;
    }
  }
}

// ----------------------------------------------RHS
for (int vi=0; vi<iel; ++vi)
{
  /* RHS source term */
  eforce[vi*numdof+dofindex] += fac_*funct_(vi)*rhsint ;

  /* convective stabilization of RHS source term */
  eforce[vi*numdof+dofindex] += taufac*conv_(vi)*rhsint ;

  /* diffusive stabilization of RHS source term */
  if (higher_order_ele) eforce[vi*numdof+dofindex] += taufac*diff_(vi)*rhsint ;
}

return;
} //Condif3Impl::Condif3CalMatStationary


/*----------------------------------------------------------------------*
 | calculate mass matrix and rhs for initializing OST          gjb 08/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::InitializeOST(
    const DRT::ELEMENTS::Condif3*   ele,
    const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ephi0,
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>& edens,
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>& epot0,
    Epetra_SerialDenseMatrix&       massmat,
    Epetra_SerialDenseVector&       rhs,
    Epetra_SerialDenseVector&       subgrdiff,
    const struct _MATERIAL*         material,
    const double                    time,
    const double                    timefac,
    const LINALG::FixedSizeSerialDenseMatrix<3,iel>& evel,
    const bool                      temperature,
    const string                    fssgd,
    const double                    frt
)
{
  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze_(0,i)=ele->Nodes()[i]->X()[0];
    xyze_(1,i)=ele->Nodes()[i]->X()[1];
    xyze_(2,i)=ele->Nodes()[i]->X()[2];
  }

  // dead load in element nodes
  BodyForce(ele,time);

  // get material constants
  GetMaterialParams(material,temperature);

  /*----------------------------------------------------------------------*/
  // calculation of instationary(!) stabilization parameter(s)
  /*----------------------------------------------------------------------*/
  CalTau(ele,subgrdiff,evel,epot0,timefac,fssgd,false,true,frt);

  /*----------------------------------------------------------------------*/
  // integration loop for one condif3 element
  /*----------------------------------------------------------------------*/

  // flag for higher order elements
  const bool higher_order_ele = SCATRA::is3DHigherOrderElement<distype>();

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(SCATRA::get3DOptimalGaussrule<distype>());

  // integration loop
  for (int iquad=0; iquad<intpoints.nquad; ++iquad)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,higher_order_ele,ele);

    // density*specific heat capacity-weighted shape functions
    densfunct_.EMultiply(funct_,edens);

    // get (density*specific heat capacity-weighted) velocity at element center
    velint_.Multiply(evel,funct_);

    //------------ get values of variables at integration point
    for (int k = 0;k<numdofpernode_;++k)     // loop of each transported sclar
    {
      // get bodyforce in gausspoint (divided by shcacp for temperature eq.)
      rhs_[k] = bodyforce_[k].Dot(funct_) / shcacp_;
    }

    // get values of all transported scalars at integration point
    for (int k=0; k<numscal_; ++k)
    {
      conint_[k] = funct_.Dot(ephi0[k]);
    }

    // get gradient of el. potential at integration point
    gradpot_.Multiply(derxy_,epot0);

    // convective part
    /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
        with  N .. form function matrix */
    conv_.MultiplyTN(derxy_,velint_);

    // migration part
    mig_.MultiplyTN(-frt,derxy_,gradpot_);

    /*-------------- perform integration for entire matrix and rhs ---*/
    for (int k=0;k<numscal_;++k) // deal with a system of transported scalars
    {
      static double             rhsint;            /* rhs at int. point     */

      // stabilization parameter
      const double taufac = tau_[k]*fac_;

      /*-------------------------------- evaluate rhs at integration point ---*/
      rhsint = rhs_[k];

      if (higher_order_ele)
      {
        for (int i=0; i<iel; i++) /* loop over nodes of element */
        {
          /* diffusive part */
          /* diffus * ( N,xx  +  N,yy +  N,zz ) */
          diff_(i) = diffus_[k] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
        } // end of loop over nodes of element
      }
      else
        diff_.Clear();

      /*--------------------------------- now build single stiffness terms ---*/
      // -------------------------------------------System matrix
      const double conv_ephi0_k = conv_.Dot(ephi0[k]);
      double diff_ephi0_k(0.0);
      if (higher_order_ele) diff_ephi0_k = diff_.Dot(ephi0[k]); // only necessary for higher order ele!

      for (int vi=0; vi<iel; ++vi)
      {
        for (int ui=0; ui<iel; ++ui)
        {
          /* Standard Galerkin terms: */
          /* transient term */
          massmat(vi*numdofpernode_+k, ui*numdofpernode_+k) += fac_*funct_(vi)*densfunct_(ui) ;
        }

        /* convective term */
        rhs[vi*numdofpernode_+k] += -(fac_*funct_(vi)*conv_ephi0_k) ;

        for (int ui=0; ui<iel; ++ui)
        {
          /* diffusive term */
          rhs[vi*numdofpernode_+k] += -(fac_*diffus_[k]*(derxy_(0, ui)*derxy_(0, vi) + derxy_(1, ui)*derxy_(1, vi) + derxy_(2, ui)*derxy_(2, vi))*(ephi0[k])(ui));
        }

        // nonlinear migration term
        rhs[vi*numdofpernode_+k] += conint_[k]*fac_*diffus_[k]*valence_[k]*mig_(vi);

        /* Stabilization terms: */
        /* 1) transient stabilization (USFEM assumed here, sign change necessary for GLS) */
        /* transient term */
        //massmat(vi, ui) += -taufac*densfunct_(vi)*densfunct_[ui] ;

        /* convective term */
        //massmat(vi, ui) += -timetaufac*densfunct_(vi)*conv[ui] ;

        /* diffusive term */
        //massmat(vi, ui) += timetaufac*densfunct_(vi)*diff[ui] ;

        /* 2) convective stabilization */
        /* transient term */
        for (int ui=0; ui<iel; ++ui)
        {
          massmat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*conv_(vi)*densfunct_(ui)*ephi0[k](ui);
        }

        /* convective term */
        rhs[vi*numdofpernode_+k] += -(taufac*conv_(vi)*conv_ephi0_k);
      }

      if (higher_order_ele)
      {
        for (int vi=0; vi<iel; ++vi)
        {
          /* diffusive term */
          rhs[vi*numdofpernode_+k] += -(-taufac*conv_(vi)*diff_ephi0_k);

          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
          /* transient term */
          for (int ui=0; ui<iel; ++ui)
          {
            massmat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*diff_(vi)*densfunct_(ui);
          }
          /* convective term */
          rhs[vi*numdofpernode_+k] += -(taufac*diff_(vi)*conv_ephi0_k); 

          /* diffusive term */
          rhs[vi*numdofpernode_+k] += -(-taufac*diff_(vi)*diff_ephi0_k);
        }
      }
      // ----------------------------------------------RHS
      for (int vi=0; vi<iel; ++vi)
      {
        /* RHS source term */
        rhs[vi*numdofpernode_+k] += fac_*funct_(vi)*rhsint ;

        /* transient stabilization of RHS source term */
        //eforce(vi) += -taufac*densfunct(vi)*rhsint ;

        /* convective stabilization of RHS source term */
        rhs[vi*numdofpernode_+k] += taufac*conv_(vi)*rhsint ;

        if (higher_order_ele)
        {
          /* diffusive stabilization of RHS source term */
          rhs[vi*numdofpernode_+k] += taufac*diff_(vi)*rhsint ;
        }

      }
    } // loop over each scalar

    if (numdofpernode_-numscal_== 1) // ELCH
    {
      // dof for el. potential have no 'acceleration'
      for (int vi=0; vi<iel; ++vi)
      {
        massmat(vi*numdofpernode_+numscal_, vi*numdofpernode_+numscal_) += 1.0;
      }
    }

  } // integration loop

  return;
} // Condif3Impl::InitializeOST


/*----------------------------------------------------------------------*
 | calculate normalized subgrid-diffusivity matrix              vg 10/08|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalcSubgridDiffMatrix(
    const DRT::ELEMENTS::Condif3*   ele,
    Epetra_SerialDenseMatrix&       sys_mat_sd,
    const double                    timefac,
    const bool                      is_stationary
    )
{
// get node coordinates
for (int i=0;i<iel;i++)
{
  xyze_(0,i)=ele->Nodes()[i]->X()[0];
  xyze_(1,i)=ele->Nodes()[i]->X()[1];
  xyze_(2,i)=ele->Nodes()[i]->X()[2];
}

/*----------------------------------------------------------------------*/
// integration loop for one condif2 element
/*----------------------------------------------------------------------*/
// gaussian points
const DRT::UTILS::IntegrationPoints3D intpoints(SCATRA::get3DOptimalGaussrule<distype>());

// integration loop
for (int iquad=0; iquad<intpoints.nquad; ++iquad)
{
  EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,false,ele);

  for (int k=0;k<numscal_;++k)
  {
    // parameter for artificial diffusivity (scaled to one here)
    double kartfac = fac_;
    if (not is_stationary) kartfac *= timefac;

    for (int vi=0; vi<iel; ++vi)
    {
      for (int ui=0; ui<iel; ++ui)
      {
        sys_mat_sd(vi*numdofpernode_+k,ui*numdofpernode_+k) += kartfac*(derxy_(0,vi)*derxy_(0,ui)+derxy_(1,vi)*derxy_(1,ui));

        /*subtract SUPG term */
        //sys_mat_sd(vi*numdofpernode_+k,ui*numdofpernode_+k) -= taufac*conv[vi]*conv[ui] ;
      }
    }
  }
} // integration loop

return;
} // Condif3Impl::CalcSubgridDiffMatrix


/*----------------------------------------------------------------------*
 | calculate matrix and rhs for electrochemistry problem      gjb 10/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalMatElch(
    Epetra_SerialDenseMatrix& emat,
    Epetra_SerialDenseVector& erhs,
    const vector<LINALG::FixedSizeSerialDenseMatrix<iel,1> >& ephinp,
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>& epotnp,
    const bool&               higher_order_ele,
    const double&             frt,
    const bool&               is_stationary,
    const double&             timefac
)
{
  // get values of all transported scalars at integration point
  for (int k=0; k<numscal_; ++k)
  {
    conint_[k] = funct_.Dot(ephinp[k]);

    // when concentration becomes zero, the coupling terms in the system matrix get lost!
    if (abs(conint_[k])<1e-18) dserror("concentration is nearly singular: %g",conint_[k]);
  }

  // get gradient of el. potential at integration point
  gradpot_.Multiply(derxy_,epotnp);

  // convective part
  /* rho * c_p * u_x * N,x  +  rho * c_p * u_y * N,y +  rho * c_p * u_z * N,z
      with  N .. form function matrix */
  conv_.MultiplyTN(derxy_,velint_);

  // migration part
  mig_.MultiplyTN(-frt,derxy_,gradpot_);

#if 0
  // DEBUG output
  cout<<endl<<"values at GP:"<<endl;
  cout<<"factor F/RT = "<<frt<<endl;
  for (int k=0;k<numscal_;++k)
  {cout<<"conint_["<<k<<"] = "<<conint_[k]<<endl;}
  for (int k=0;k<3;++k)
  {cout<<"gradpot_["<<k<<"] = "<<gradpot_[k]<<endl;}
#endif
  // some 'working doubles'
  static double diffus_valence_k(0.0);
  double rhsint(0.0);  // rhs at int. point
  // integration factors and coefficients of single terms
  static double timefacfac(0.0);
  static double timetaufac(0.0);
  static double taufac(0.0);

  for (int k = 0; k < numscal_;++k) // loop over all transported sclars
  {
    // stabilization parameters
    taufac = tau_[k]*fac_;

    // factor D_k * z_k
    diffus_valence_k = diffus_[k]*valence_[k];

    if (is_stationary)
    {
      timefacfac  = fac_;
      timetaufac  = taufac;
      rhsint = rhs_[k];     // set rhs at integration point
    }
    else
    {
      timefacfac  = timefac * fac_;
      timetaufac  = timefac * taufac;
      rhsint = hist_[k] + rhs_[k]*timefac;     // set rhs at integration point
    }

    // compute gradient of scalar k at integration point
    gradphi_.Multiply(derxy_,ephinp[k]);

    if (higher_order_ele)
    {
      for (int i=0; i<iel; i++)
      {
        // diffusive part
        /* diffus * ( N,xx  +  N,yy +  N,zz ) */
        diff_(i) = diffus_[k] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));

        /* reactive part of migration*/
        /* diffus * ( N,xx  +  N,yy +  N,zz ) */
        //migr_[i] = diffus_[k] * funct_[i] * (derxy2_(0,i) + derxy2_(1,i) + derxy2_(2,i));
      }
    }

    const double frt_timefacfac_diffus_valence_k_conint_k = frt*timefacfac*diffus_valence_k*conint_[k];

    // ----------------------------------------matrix entries
    for (int vi=0; vi<iel; ++vi)
    {
      const double timetaufac_conv_eff_vi = timetaufac*(conv_(vi)+diffus_valence_k*mig_(vi));
      const double timefacfac_funct_vi = timefacfac*funct_(vi);
      const double timefacfac_diffus_valence_k_mig_vi = timefacfac*diffus_valence_k*mig_(vi);
      const double valence_k_fac_funct_vi = valence_[k]*fac_*funct_(vi);

      for (int ui=0; ui<iel; ++ui)
      {
        /* Standard Galerkin terms: */
        /* convective term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timefacfac_funct_vi*conv_(ui) ;

        /* diffusive term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timefacfac*diffus_[k]*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui) + derxy_(2, vi)*derxy_(2, ui));

        /* migration term (directional derivatives) */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) -= timefacfac_diffus_valence_k_mig_vi*funct_(ui);
        emat(vi*numdofpernode_+k,ui*numdofpernode_+numscal_) += frt_timefacfac_diffus_valence_k_conint_k*(derxy_(0, vi)*derxy_(0, ui) + derxy_(1, vi)*derxy_(1, ui) + derxy_(2, vi)*derxy_(2, ui));

        /* electroneutrality condition */
        emat(vi*numdofpernode_+numscal_, ui*numdofpernode_+k) += valence_k_fac_funct_vi*densfunct_(ui);

        /* Stabilization term: */
        /* 0) transient stabilization */
        // not implemented

        /* 1) convective stabilization */

        /* convective term */
        emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timetaufac_conv_eff_vi*(conv_(ui)+diffus_valence_k*mig_(ui));

      } // for ui
    } // for vi

    if (higher_order_ele)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        for (int ui=0; ui<iel; ++ui)
        {
          /* diffusive term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += -timetaufac*(conv_(vi)+diffus_valence_k*mig_(vi))*diff_(ui) ;

          /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */

          /* convective term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += timetaufac*diff_(vi)*(conv_(ui)+diffus_valence_k*mig_(ui));

          /* diffusive term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) -= timetaufac*diff_(vi)*diff_(ui) ;
        } // for ui
      } // for vi

    } // higher_order_ele

    // ----------------------------------------------RHS
    const double conv_ephinp_k = conv_.Dot(ephinp[k]);
    const double Dkzk_mig_ephinp_k = diffus_valence_k*(mig_.Dot(ephinp[k]));
    const double conv_eff_k = conv_ephinp_k + Dkzk_mig_ephinp_k;
    const double densfunct_ephinp_k = densfunct_.Dot(ephinp[k]);
    double diff_ephinp_k(0.0);
    if (higher_order_ele) diff_ephinp_k = diff_.Dot(ephinp[k]); // only necessary for higher order ele!

    // compute residual of strong form for stabilization
    double taufacresidual = taufac*rhsint - timetaufac*(conv_eff_k + diff_ephinp_k);
    if (!is_stationary) // add transient term to the residual
      taufacresidual -= taufac*densfunct_ephinp_k;

    //------------residual formulation (Newton iteration)
    for (int vi=0; vi<iel; ++vi)
    {
      // RHS source term
      erhs[vi*numdofpernode_+k] += fac_*funct_(vi)*rhsint ;

      // nonlinear migration term
      erhs[vi*numdofpernode_+k] += conint_[k]*timefacfac*diffus_valence_k*mig_(vi);

      // convective term
      erhs[vi*numdofpernode_+k] -= timefacfac*funct_(vi)*conv_ephinp_k;

      // diffusive term
      erhs[vi*numdofpernode_+k] -= timefacfac*diffus_[k]*(gradphi_(0)*derxy_(0, vi) + gradphi_(1)*derxy_(1, vi)+ gradphi_(2)*derxy_(2, vi));

      // electroneutrality condition
      // for incremental formulation, there is the residuum on the rhs! : 0-ENC*phi_i
      erhs[vi*numdofpernode_+numscal_] -= valence_[k]*fac_*funct_(vi)*densfunct_ephinp_k;

      // Stabilization terms:

      // 0) transient stabilization
      // not implemented

      // 1) convective stabilization
      erhs[vi*numdofpernode_+k] += (conv_(vi)+diffus_valence_k*mig_(vi)) * taufacresidual;

    } // for vi

    if (higher_order_ele)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        dserror("higher order terms not yet tested");

        // 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS)
        erhs[vi*numdofpernode_+k] += diff_(vi)*taufacresidual ;

      } // for vi
    } // higher_order_ele

    // -----------------------------------INSTATIONARY TERMS
    if (!is_stationary)
    {
      for (int vi=0; vi<iel; ++vi)
      {
        const double fac_funct_vi = fac_*funct_(vi);
        for (int ui=0; ui<iel; ++ui)
        {
          /* Standard Galerkin terms: */
          /* transient term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += fac_funct_vi*densfunct_(ui) ;

          /* 1) convective stabilization */
          /* transient term */
          emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*(conv_(vi)+diffus_valence_k*mig_(vi))*densfunct_(ui);

          if (higher_order_ele)
          {
            /* 2) diffusive stabilization (USFEM assumed here, sign change necessary for GLS) */
            /* transient term */
            emat(vi*numdofpernode_+k, ui*numdofpernode_+k) += taufac*diff_(vi)*densfunct_(ui) ;
          }
        } // for ui

        // residuum on RHS:

        /* Standard Galerkin terms: */
        /* transient term */
        erhs[vi*numdofpernode_+k] -= fac_funct_vi*densfunct_ephinp_k;

      } // for vi
    } // instationary case

  } // loop over scalars

  return;
} // Condif3Impl::CalMatElch


/*---------------------------------------------------------------------*
 |  calculate error compared to analytical solution           gjb 10/08|
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Condif3Impl<distype>::CalErrorComparedToAnalytSolution(
    const DRT::ELEMENTS::Condif3*   ele,
    ParameterList& params,
    const LINALG::FixedSizeSerialDenseMatrix<iel,2>& ephinp,
    const LINALG::FixedSizeSerialDenseMatrix<iel,1>& epotnp,
    Epetra_SerialDenseVector& errors,
    struct _MATERIAL* material
)
{
  //at the moment, there is only one analytical test problem available!
  if (params.get<string>("action") != "calc_elch_kwok_error")
    dserror("Unknown analytical solution");

  //------------------------------------------------- Kwok et Wu,1995
  //   Reference:
  //   Kwok, Yue-Kuen and Wu, Charles C. K.
  //   "Fractional step algorithm for solving a multi-dimensional diffusion-migration equation"
  //   Numerical Methods for Partial Differential Equations
  //   1995, Vol 11, 389-397

  // get node coordinates
  for (int i=0;i<iel;i++)
  {
    xyze_(0,i)=ele->Nodes()[i]->X()[0];
    xyze_(1,i)=ele->Nodes()[i]->X()[1];
    xyze_(2,i)=ele->Nodes()[i]->X()[2];
  }

  // set constants for analytical solution
  const double t = params.get("total time",-1.0);
  dsassert (t >= 0.0, "no total time for error calculation");
  const double frt = params.get<double>("frt");

  // get material constants
  GetMaterialParams(material,false);


  // working arrays
  double                                  potint;
  LINALG::FixedSizeSerialDenseMatrix<2,1> conint;
  LINALG::FixedSizeSerialDenseMatrix<3,1> xint;
  LINALG::FixedSizeSerialDenseMatrix<2,1> c;
  double                                  deltapot;
  LINALG::FixedSizeSerialDenseMatrix<2,1> deltacon(true);

  // integration points
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point; // for cos/sin
  const DRT::UTILS::IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,false,ele);

    // get both concentration solutions at integration point
    conint.MultiplyTN(ephinp,funct_);

    // get el. potential solution at integration point
    potint = funct_.Dot(epotnp);

    // get global coordinate of integration point
    xint.Multiply(xyze_,funct_);

    // compute varous constants
    const double d= (frt*(diffus_[0]*valence_[0]-diffus_[1]*valence_[1]));
    if (abs(d) == 0.0) dserror("division by zero");
    const double D = frt*(valence_[0]*diffus_[0]*diffus_[1] - valence_[1]*diffus_[1]*diffus_[0])/d;

    // compute analytical concentrations for the problem of Kwok(1995)
    const double A0 = 5.0;
    const double m = 2.0;
    const double n = 2.0;
    c(0) = A0 + ((cos(m*PI*xint(0))*cos(n*PI*xint(1)))*exp((-D)*(m*m + n*n)*t*PI*PI));
    c(1) = (-valence_[0]/valence_[1])* c(0);

    // compute analytical el. potential
    const double c_0_0_t = A0 + exp((-D)*(m*m + n*n)*t*PI*PI);
    const double pot = ((diffus_[1]-diffus_[0])/d) * log(c(0)/c_0_0_t);

    // compute differences between analytical solution and numerical solution
    deltapot = potint - pot;

    deltacon = conint;
    deltacon -= c;

    // add square to L2 error
    errors[0] += deltacon(0)*deltacon(0)*fac_;
    errors[1] += deltacon(1)*deltacon(1)*fac_;
    errors[2] += deltapot*deltapot*fac_;

  } // end of loop over integration points

  return;
} // Condif3Impl::CalErrorComparedToAnalytSolution

#endif
#endif
