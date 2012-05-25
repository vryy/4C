/*!------------------------------------------------------------------------------------------------*
\file topopt_optimizer_ele_impl.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_optimizer_ele_impl.H"
#include "topopt_optimizer_ele_parameter.H"

#include "../drt_lib/drt_utils.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptImplInterface* DRT::ELEMENTS::TopOptImplInterface::Impl(
  const DRT::Element* ele
)
{
  switch (ele->Shape())
  {
  case DRT::Element::hex8:
  {
    return TopOptImpl<DRT::Element::hex8>::Instance();
  }
  case DRT::Element::hex20:
  {
    return TopOptImpl<DRT::Element::hex20>::Instance();
  }
  case DRT::Element::hex27:
  {
    return TopOptImpl<DRT::Element::hex27>::Instance();
  }
  case DRT::Element::tet4:
  {
    return TopOptImpl<DRT::Element::tet4>::Instance();
  }
  case DRT::Element::tet10:
  {
    return TopOptImpl<DRT::Element::tet10>::Instance();
  }
  case DRT::Element::quad4:
  {
    return TopOptImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return TopOptImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return TopOptImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return TopOptImpl<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return TopOptImpl<DRT::Element::tri6>::Instance();

  case DRT::Element::line2:
  {
    return TopOptImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return TopOptImpl<DRT::Element::line3>::Instance();
  }
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(ele->Shape()).c_str());
  }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptImpl<distype> * DRT::ELEMENTS::TopOptImpl<distype>::Instance(
    bool create)
{
  static TopOptImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new TopOptImpl<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptImpl<distype>::TopOptImpl()
: intpoints_( distype ),
xsi_(true),
det_(0.0),
fac_(0.0),
visc_(0.0),
reacoeff_(0.0),
dens_(0.0),
is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateObjective(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm
)
{
  return EvaluateObjective(
      ele,
      params,
      discretization,
      lm,
      intpoints_
  );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateObjective(
  DRT::Element*                 ele,
  ParameterList&                params,
  DRT::Discretization&          discretization,
  vector<int>&                  lm,
  DRT::UTILS::GaussIntegration& intpoints
)
{
  // TODO coming...

  // work is done
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradient(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{
  return EvaluateGradient(
      ele,
      params,
      discretization,
      lm,
      elevec1_epetra,
      intpoints_
  );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptImpl<distype>::EvaluateGradient(
  DRT::Element*                 ele,
  ParameterList&                params,
  DRT::Discretization&          discretization,
  vector<int>&                  lm,
  Epetra_SerialDenseVector&     elevec1_epetra,
  DRT::UTILS::GaussIntegration& intpoints

  )
{
  RCP<DRT::Discretization> fluiddis = params.get<RCP<DRT::Discretization> >("fluiddis");

  RCP<map<int,RCP<Epetra_Vector> > > fluidvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("fluidvel");
  RCP<map<int,RCP<Epetra_Vector> > > adjointvels = params.get<RCP<map<int,RCP<Epetra_Vector> > > >("adjointvel");

  map<int,LINALG::Matrix<nsd_,nen_> > fluidnodalvels;
  map<int,LINALG::Matrix<nsd_,nen_> > adjointnodalvels;

  LINALG::Matrix<nsd_,nen_> fluidnodalvel;
  LINALG::Matrix<nsd_,nen_> adjointnodalvel;

  vector<int> fluidlm;
  {
    vector<int> lmowner; // dummy for function call
    vector<int> lmstride; // dummy for function call
    ele->LocationVector(*fluiddis,fluidlm,lmowner,lmstride);
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=fluidvels->begin();
      i!=fluidvels->end();i++)
  {
    ExtractValuesFromGlobalVector(fluiddis,fluidlm,&fluidnodalvel,NULL,i->second);
    fluidnodalvels.insert(pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,fluidnodalvel));
  }

  for (map<int,RCP<Epetra_Vector> >::iterator i=adjointvels->begin();
      i!=adjointvels->end();i++)
  {
    ExtractValuesFromGlobalVector(fluiddis,fluidlm,&adjointnodalvel,NULL,i->second);
    adjointnodalvels.insert(pair<int,LINALG::Matrix<nsd_,nen_> >(i->first,fluidnodalvel));
  }

//  std::vector<double> nodalfluidvel(dim);
//  std::vector<double> nodaladjointvel(dim);
//
//  double dissipation_fac = 0.0;
//  if (fluidParams_->get<bool>("OBJECTIVE_DISSIPATION")) dissipation_fac = fluidParams_->get<double>("DISSIPATION_FAC");
//
//  for (int inode=0;inode<discret_->NumMyRowNodes();inode++)
//  {
//    node = discret_->lRowNode(inode);
//
//    vector<int> lm = discret_->Dof(node);
//    lm.pop_back(); // delete pressure dof
//
//    value = 0.0;
//
//    if (fluidParams_->get<int>("time int algo")==INPAR::FLUID::timeint_stationary)
//    {
//      fluidvel = vel_->find(1)->second;
//      adjointvel = adjointvel_->find(1)->second;
//
//      DRT::UTILS::ExtractMyValues(*fluidvel,nodalfluidvel,lm);
//      DRT::UTILS::ExtractMyValues(*adjointvel,nodaladjointvel,lm);
//
//      for (int idim=0;idim<dim;idim++)
//        value +=
//            dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
//            -nodalfluidvel[idim]*nodaladjointvel[idim]; // adjoint part
//    }
//    else
//    {
//      for (int timestep=0;timestep<=fluidParams_->get<int>("max number timesteps");timestep++)
//      {
//        fluidvel = vel_->find(timestep)->second;
//        adjointvel = adjointvel_->find(timestep)->second;
//
//        DRT::UTILS::ExtractMyValues(*fluidvel,nodalfluidvel,lm);
//        DRT::UTILS::ExtractMyValues(*adjointvel,nodaladjointvel,lm);
//        if (inode==6) cout << "initial value is " << value << endl;
//        if (timestep!=0 && timestep!=fluidParams_->get<int>("max number timesteps")) // default case, weight 1
//        {
//          for (int idim=0;idim<dim;idim++)
//          {
//            value +=
//                dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
//                -nodalfluidvel[idim]*nodaladjointvel[idim]; // adjoint part
//          }
//        }
//        else // first and last time step, weight 0.5
//        {
//          for (int idim=0;idim<dim;idim++)
//          {
//            value += 0.5*(
//                dissipation_fac*nodalfluidvel[idim]*nodalfluidvel[idim] // dissipation part
//                -nodalfluidvel[idim]*nodaladjointvel[idim]); // adjoint part
//          }
//        }
//      }
//      value *= dt; // scale with time step size
//    }
//
//    int err = obj_grad_->SumIntoMyValue(inode,0,value);
//    if (err)
//      dserror("error while adding value to gradient of objective");
//  }
  // work is done
  return 0;
}


///*----------------------------------------------------------------------*
//|  calculate system matrix and rhs (public)                 g.bau 08/08|
//*----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::TopOptImpl<distype>::Sysmat(
//  DRT::Element*                         ele, ///< the element those matrix is calculated
//  Epetra_SerialDenseMatrix&             emat,///< element matrix to calculate
//  Epetra_SerialDenseVector&             erhs, ///< element rhs to calculate
//  Epetra_SerialDenseVector&             subgrdiff, ///< subgrid-diff.-scaling vector
//  const double                          time, ///< current simulation time
//  const double                          dt, ///< current time-step length
//  const double                          timefac, ///< time discretization factor
//  const double                          alphaF, ///< factor for generalized-alpha time integration
//  const enum INPAR::SCATRA::AssgdType   whichassgd, ///< all-scale subgrid-diffusivity definition
//  const enum INPAR::SCATRA::FSSUGRDIFF  whichfssgd, ///< fine-scale subgrid-diffusivity definition
//  const bool                            assgd, ///< all-scale subgrid-diff. flag
//  const bool                            fssgd, ///< fine-scale subgrid-diff. flag
//  const double                          Cs, ///< Smagorinsky constant
//  const double                          tpn, ///< turbulent Prandtl number
//  const double                          Csgs_sgvel, ///< parameter of multifractal subgrid-scales
//  const double                          alpha, ///< grid-filter to test-filter ratio
//  const bool                            calc_N, ///< flag to activate calculation of N
//  const double                          N_vel, ///< value for N if not calculated
//  const enum INPAR::FLUID::RefVelocity  refvel, ///< reference velocity
//  const enum INPAR::FLUID::RefLength    reflength, ///< reference length
//  const double                          c_nu, ///< scaling for Re
//  const bool                            nwl, ///< flag to activate near-wall limit
//  const double                          Csgs_sgphi, ///< parameter of multifractal subgrid-scales
//  const double                          c_diff, ///< scaling for Re*Pr
//  const bool                            BD_gp, ///< evaluation of model coefficient at gp
//  const double                          frt, ///< factor F/RT needed for ELCH calculations
//  const enum INPAR::SCATRA::ScaTraType  scatratype ///< type of scalar transport problem
//  )
//{
//}



///*----------------------------------------------------------------------*
//  |  get the material constants  (private)                      gjb 10/08|
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::TopOptImpl<distype>::GetMaterialParams(
//  const DRT::Element*  ele
//)
//{
//  // get the material
//  RefCountPtr<MAT::Material> material = ele->Material();
//
//  if (material->MaterialType() == INPAR::MAT::m_scatra)
//  {
//    const MAT::ScatraMat* actmat = static_cast<const MAT::ScatraMat*>(material.get());
//
//    dsassert(numdofpernode_==1,"more than 1 dof per node for SCATRA material");
//
//    // get constant diffusivity
//    diffus_[0] = actmat->Diffusivity();
//
//    // in case of reaction with (non-zero) constant coefficient:
//    // read coefficient and set reaction flag to true
//    reacoeff_[0] = actmat->ReaCoeff();
//    if (reacoeff_[0] > EPS14) is_reactive_ = true;
//    if (reacoeff_[0] < -EPS14)
//      dserror("Reaction coefficient is not positive: %f",0, reacoeff_[0]);
//
//    reacoeffderiv_[0] = reacoeff_[0];
//
//    // set specific heat capacity at constant pressure to 1.0
//    shc_ = 1.0;
//
//    // set temperature rhs for reactive equation system to zero
//    reatemprhs_[0] = 0.0;
//
//    // set density at various time steps and density gradient factor to 1.0/0.0
//    densn_[0]       = 1.0;
//    densnp_[0]      = 1.0;
//    densam_[0]      = 1.0;
//    densgradfac_[0] = 0.0;
//
//    // in case of multifrcatal subgrid-scales, read Schmidt number
//    if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales or sgvel_)
//    {
//      double scnum = actmat->ScNum();
//      visc_ = scnum * diffus_[0];
//    }
//  }
//  else dserror("Material type is not supported");
//
//// check whether there is negative (physical) diffusivity
//  if (diffus_[0] < -EPS15) dserror("negative (physical) diffusivity");
//
//  return;
//} //TopOptImpl::GetMaterialParams



///*----------------------------------------------------------------------*
//  |  Integrate shape functions over domain (private)           gjb 07/09 |
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::TopOptImpl<distype>::IntegrateShapeFunctions(
//  const DRT::Element*             ele,
//  Epetra_SerialDenseVector&       elevec1,
//  const Epetra_IntSerialDenseVector& dofids
//)
//{
//  // integrations points and weights
//  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
//
//  // safety check
//  if (dofids.M() < numdofpernode_)
//    dserror("Dofids vector is too short. Received not enough flags");
//
//  // loop over integration points
//  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
//  {
//    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,gpid,ele->Id());
//
//    // compute integral of shape functions (only for dofid)
//    for (int k=0;k<numdofpernode_;k++)
//    {
//      if (dofids[k] >= 0)
//      {
//        for (int node=0;node<nen_;node++)
//        {
//          elevec1[node*numdofpernode_+k] += funct_(node) * fac;
//        }
//      }
//    }
//
//  } //loop over integration points
//
//  return;
//
//} //TopOptImpl<distype>::IntegrateShapeFunction
//
//
///*----------------------------------------------------------------------*
//  | evaluate shape functions and derivatives at int. point     gjb 08/08 |
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//double DRT::ELEMENTS::TopOptImpl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
//  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
//  const int                                    iquad,      ///< id of current Gauss point
//  const int                                    eleid       ///< the element id
//  )
//{
//  // coordinates of the current integration point
//  const double* gpcoord = (intpoints.IP().qxg)[iquad];
//  for (int idim=0;idim<nsd_;idim++)
//    xsi_(idim) = gpcoord[idim];
//
//  // shape functions and their first derivatives
//  DRT::UTILS::shape_function<distype>(xsi_,funct_);
//  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
//  if (use2ndderiv_)
//  {
//    // get the second derivatives of standard element at current GP
//    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
//  }
//
//  // compute Jacobian matrix and determinant
//  // actually compute its transpose....
//  /*
//    +-            -+ T      +-            -+
//    | dx   dx   dx |        | dx   dy   dz |
//    | --   --   -- |        | --   --   -- |
//    | dr   ds   dt |        | dr   dr   dr |
//    |              |        |              |
//    | dy   dy   dy |        | dx   dy   dz |
//    | --   --   -- |   =    | --   --   -- |
//    | dr   ds   dt |        | ds   ds   ds |
//    |              |        |              |
//    | dz   dz   dz |        | dx   dy   dz |
//    | --   --   -- |        | --   --   -- |
//    | dr   ds   dt |        | dt   dt   dt |
//    +-            -+        +-            -+
//  */
//
//  xjm_.MultiplyNT(deriv_,xyze_);
//  const double det = xij_.Invert(xjm_);
//
//  if (det < 1E-16)
//    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det);
//
//  // set integration factor: fac = Gauss weight * det(J)
//  const double fac = intpoints.IP().qwgt[iquad]*det;
//
//  // compute global derivatives
//  derxy_.Multiply(xij_,deriv_);
//
//  // compute second global derivatives (if needed)
//  if (use2ndderiv_)
//  {
//    // get global second derivatives
//    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
//  }
//  else
//    derxy2_.Clear();
//
//  // return integration factor for current GP: fac = Gauss weight * det(J)
//  return fac;
//
//} //TopOptImpl::CalcSubgrVelocity
//
//
//
///*---------------------------------------------------------------------*
//  |  calculate error compared to analytical solution           gjb 10/08|
//  *---------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::TopOptImpl<distype>::CalErrorComparedToAnalytSolution(
//  const DRT::Element*                   ele,
//  const enum INPAR::SCATRA::ScaTraType  scatratype,
//  ParameterList&                        params,
//  Epetra_SerialDenseVector&             errors
//  )
//{
//  return;
//} // TopOptImpl::CalErrorComparedToAnalytSolution
//
//
//
///*----------------------------------------------------------------------*
//  |  calculation of characteristic element length               vg 01/11 |
//  *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//double DRT::ELEMENTS::TopOptImpl<distype>::CalcCharEleLength(
//  const double  vol,
//  const double  vel_norm
//  )
//{
//  //---------------------------------------------------------------------
//  // various definitions for characteristic element length
//  //---------------------------------------------------------------------
//  // a) streamlength due to Tezduyar et al. (1992) -> default
//  // normed velocity vector
//  LINALG::Matrix<nsd_,1> velino;
//  if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,convelint_);
//  else
//  {
//    velino.Clear();
//    velino(0,0) = 1;
//  }
//
//  // get streamlength using the normed velocity at element centre
//  LINALG::Matrix<nen_,1> tmp;
//  tmp.MultiplyTN(derxy_,velino);
//  const double val = tmp.Norm1();
//  const double hk = 2.0/val; // h=streamlength
//
//  // b) volume-equivalent diameter (warning: 3-D formula!)
//  //hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
//
//  // c) cubic/square root of element volume/area or element length (3-/2-/1-D)
//  // cast dimension to a double varibale -> pow()
//  //const double dim = (double) nsd_;
//  //hk = pow(vol,1/dim);
//
//  return hk;
//}



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptImpl<distype>::ExtractValuesFromGlobalVector(
    Teuchos::RCP<DRT::Discretization> discretization, ///< discretization
    const vector<int>&           lm,                  ///<
    LINALG::Matrix<nsd_,nen_> *  matrixtofill,        ///< vector field
    LINALG::Matrix<nen_,1> *     vectortofill,        ///< scalar field
    RCP<Epetra_Vector>&          globalvector         ///< global vector
) const
{
  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*globalvector,mymatrix,lm);
cout << "here not" << endl;
  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*(nsd_+1))];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*(nsd_+1))];
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TopOptBoundaryImplInterface* DRT::ELEMENTS::TopOptBoundaryImplInterface::Impl(
    const DRT::Element* ele
)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return TopOptBoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return TopOptBoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return TopOptBoundaryImpl<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return TopOptBoundaryImpl<DRT::Element::tri3>::Instance();
  }
  /*  case DRT::Element::tri6:
  {
    return TopOptBoundaryImpl<DRT::Element::tri6>::Instance();
  }*/
  case DRT::Element::line2:
  {
    return TopOptBoundaryImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return TopOptBoundaryImpl<DRT::Element::line3>::Instance();
  }
  case DRT::Element::nurbs2:    // 1D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:    // 1D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:    // 2D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs9:    // 2D nurbs boundary element
  {
    return TopOptBoundaryImpl<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptBoundaryImpl<distype> * DRT::ELEMENTS::TopOptBoundaryImpl<distype>::Instance(
    bool create
)
{
  static TopOptBoundaryImpl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new TopOptBoundaryImpl<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::TopOptBoundaryImpl<distype>::TopOptBoundaryImpl
()
: intpoints_( distype ),
xsi_(true),
det_(0.0),
fac_(0.0),
visc_(0.0),
reacoeff_(0.0),
dens_(0.0),
is_higher_order_ele_(false)
{
  // pointer to class FluidEleParameter (access to the general parameter)
  optiparams_ = DRT::ELEMENTS::TopOptParam::Instance();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryObjective(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm
  )
{
  // TODO coming...

  // work is done
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvaluateBoundaryGradient(
  DRT::Element*              ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  Epetra_SerialDenseVector&  elevec1_epetra
  )
{
  // TODO coming...

  // work is done
  return 0;
}





///*----------------------------------------------------------------------*
// | evaluate shape functions and int. factor at int. point     gjb 01/09 |
// *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//double DRT::ELEMENTS::TopOptBoundaryImpl<distype>::EvalShapeFuncAndIntFac(
//    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
//    const int                                    iquad,      ///< id of current Gauss point
//    const int                                    eleid,      ///< the element id
//    LINALG::Matrix<1 + nsd_,1>*         normalvec ///< normal vector at Gauss point(optional)
//)
//{
//  // coordinates of the current integration point
//  const double* gpcoord = (intpoints.IP().qxg)[iquad];
//  for (int idim=0;idim<nsd_;idim++)
//  {xsi_(idim) = gpcoord[idim];}
//
//  if(not DRT::NURBS::IsNurbs(distype))
//  {
//    // shape functions and their first derivatives
//    DRT::UTILS::shape_function<distype>(xsi_,funct_);
//    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
//  }
//  else // nurbs elements are always somewhat special...
//  {
//    DRT::NURBS::UTILS::nurbs_get_funct_deriv(
//        funct_  ,
//        deriv_  ,
//        xsi_    ,
//        myknots_,
//        weights_,
//        distype );
//  }
//
//  // the metric tensor and the area of an infinitesimal surface/line element
//  // optional: get normal at integration point as well
//  // Note: this is NOT yet a unit normal. Its norm corresponds to the area/length of the element
//  double drs(0.0);
//  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor_,drs,normalvec);
//
//  // for nurbs elements the normal vector must be scaled with a special orientation factor!!
//  if(DRT::NURBS::IsNurbs(distype))
//  {
//    if (normalvec != NULL)
//      normal_.Scale(normalfac_);
//  }
//
//  // return the integration factor
//  return intpoints.IP().qwgt[iquad] * drs;
//}
//
//
//
///*----------------------------------------------------------------------*
// |  Integrate shapefunctions over surface (private)           gjb 02/09 |
// *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::IntegrateShapeFunctions(
//    const DRT::Element*        ele,
//    ParameterList&             params,
//    Epetra_SerialDenseVector&  elevec1,
//    const bool                 addarea
//)
//{
//  // access boundary area variable with its actual value
//  double boundaryint = params.get<double>("boundaryint");
//
//  // integrations points and weights
//  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
//
//  // loop over integration points
//  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
//  {
//    const double fac = EvalShapeFuncAndIntFac(intpoints,gpid,ele->Id());
//
//    // compute integral of shape functions
//    for (int node=0;node<nen_;++node)
//    {
//      for (int k=0; k< numscal_; k++)
//      {
//        elevec1[node*numdofpernode_+k] += funct_(node) * fac;
//      }
//    }
//
//    if (addarea)
//    {
//      // area calculation
//      boundaryint += fac;
//    }
//
//  } //loop over integration points
//
//  // add contribution to the global value
//  params.set<double>("boundaryint",boundaryint);
//
//  return;
//
//} //TopOptBoundaryImpl<distype>::IntegrateShapeFunction



/*---------------------------------------------------------------------------------*
 | extract element data from global vector                        winklmaier 03/12 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::TopOptBoundaryImpl<distype>::ExtractValuesFromGlobalVector(
    Teuchos::RCP<DRT::Discretization> discretization, ///< discretization
    const vector<int>&           lm,                  ///<
    LINALG::Matrix<nsd_,nen_> *  matrixtofill,        ///< vector field
    LINALG::Matrix<nen_,1> *     vectortofill,        ///< scalar field
    RCP<Epetra_Vector>&          globalvector         ///< global vector
) const
{
  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*globalvector,mymatrix,lm);

  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*(nsd_+1))];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*(nsd_+1))];
  }
}


