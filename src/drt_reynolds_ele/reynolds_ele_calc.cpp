/*--------------------------------------------------------------------------*/
/*!
\file reynolds_ele_calc.cpp

\brief main file containing routines for calculation of reynolds element

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "reynolds_ele_calc.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_scatra_ele/scatra_ele_calc_utils.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_mat/material.H"
#include "../drt_mat/scatra_mat.H"

#include "reynolds_ele_parameter.H"

/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::ReynoldsEleCalc(const std::string& disname)
  : reynoldspara_(DRT::ELEMENTS::ReynoldsEleParameter::Instance(disname)),            // standard parameter list
    eprenp_(true),  // initialized to zero
    xsi_(true),     // initialized to zero
    xyze_(true),    // initialized to zero
    funct_(true),   // initialized to zero
    deriv_(true),   // initialized to zero
    derxy_(true),   // initialized to zero
    xjm_(true),     // initialized to zero
    xij_(true),     // initialized to zero
    bodyforce_(true), // size of vector
    diffmanager_(Teuchos::rcp(new ReynoldsEleDiffManager())),           // diffusion manager for diffusivity
    reynoldsvarmanager_(Teuchos::rcp(new ReynoldsEleInternalVariableManager<nsd_,nen_>())),   // internal variable manager
    eid_(0),
    ele_(NULL)
{
  dsassert(nsd_ >= nsd_ele_,"problem dimension has to be equal or larger than the element dimension!");

  return;
}

/*----------------------------------------------------------------------*
 | singleton access method                                  wirtz 10/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>* DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::Instance(
    const std::string& disname,
    const ReynoldsEleCalc* delete_me
    )
{
  static std::map<std::string,ReynoldsEleCalc<distype,probdim>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ReynoldsEleCalc<distype,probdim>(disname);
  }
  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for( typename std::map<std::string,ReynoldsEleCalc<distype,probdim>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::Done()
{
  // delete instance
  Instance("",this);

  return;
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::Evaluate(
  DRT::Element*                 ele,
  Teuchos::ParameterList&       params,
  DRT::Discretization&          discretization,
  DRT::Element::LocationArray&  la,
  Epetra_SerialDenseMatrix&     elemat1_epetra,
  Epetra_SerialDenseMatrix&     elemat2_epetra,
  Epetra_SerialDenseVector&     elevec1_epetra,
  Epetra_SerialDenseVector&     elevec2_epetra,
  Epetra_SerialDenseVector&     elevec3_epetra
  )
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  if(SetupCalc(ele,discretization) == -1)
    return 0;

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractElementAndNodeValues(ele,params,discretization,la);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  Sysmat(ele,elemat1_epetra,elevec1_epetra);

  return 0;
}

/*----------------------------------------------------------------------*
 | setup element evaluation                                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::SetupCalc(
    DRT::Element*               ele,
    DRT::Discretization&        discretization
    )
{
  // get element coordinates
  ReadElementCoordinates(ele);

  // set element id
  eid_ = ele->Id();
  // set element
  ele_ = ele;

  return 0;
}

/*----------------------------------------------------------------------*
 | read element coordinates                                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::ReadElementCoordinates(
    const DRT::Element*     ele
    )
{
  // Directly copy the coordinates since in 3D the transformation is just the identity
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  return;
} //ReynoldsEleCalc::ReadElementCoordinates

/*----------------------------------------------------------------------*
 | extract element based or nodal values                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
)
{
  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> prenp = discretization.GetState("phinp");
  if (prenp==Teuchos::null)
    dserror("Cannot get state vector 'phinp'");

  //values of pressure field are always in first dofset
  const std::vector<int>&    lm = la[0].lm_;
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_,1> >(*prenp,eprenp_,lm);

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele);

  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 wirtz 10/15 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::Sysmat(
  DRT::Element*                         ele,        ///< the element whose matrix is calculated
  Epetra_SerialDenseMatrix&             emat,       ///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs        ///< element rhs to calculate
  )
{
  //----------------------------------------------------------------------
  // get material and stabilization parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

    // compute rhs containing bodyforce (divided by specific heat capacity) and,
    // for temperature equation, the time derivative of thermodynamic pressure,
    // if not constant, and for temperature equation of a reactive
    // equation system, the reaction-rate term
    double rhsint(0.0);
    GetRhsInt(rhsint,densnp);

    //----------------------------------------------------------------
    // standard Galerkin terms
    //----------------------------------------------------------------

    // stabilization parameter and integration factors
    const double timefacfac = fac; // ************************** works only for stationary problems!

    //----------------------------------------------------------------
    // 1) element matrix: stationary terms
    //----------------------------------------------------------------

    // calculation of diffusive element matrix
    CalcMatDiff(emat,timefacfac);

    //----------------------------------------------------------------
    // 5) element right hand side
    //----------------------------------------------------------------
    //----------------------------------------------------------------
    // computation of bodyforce (and potentially history) term,
    // residual, integration factors and standard Galerkin transient
    // term (if required) on right hand side depending on respective
    // (non-)incremental stationary or time-integration scheme
    //----------------------------------------------------------------
    double rhsfac    = fac; // ************************** works only for stationary problems!

    //----------------------------------------------------------------
    // standard Galerkin transient, old part of rhs and bodyforce term
    //----------------------------------------------------------------
    CalcRHSHistAndSource(erhs,fac,rhsint);

    //----------------------------------------------------------------
    // standard Galerkin terms on right hand side
    //----------------------------------------------------------------

    // diffusive term
    CalcRHSDiff(erhs,rhsfac);

  }// end loop Gauss points

  return;
}


/*------------------------------------------------------------------------------*
 | set internal variables                                           wirtz 10/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::SetInternalVariablesForMatAndRHS()
{
  reynoldsvarmanager_->SetInternalVariables(funct_,derxy_,eprenp_);
  return;
}

/*----------------------------------------------------------------------*
 |  get the material constants  (private)                     wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
// get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  Materials(material,densn,densnp,densam,visc,iquad);

  return;
} //ReynoldsEleCalc::GetMaterialParams

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                   wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point

  )
{

  // todo: reynolds material
  switch(material->MaterialType())
  {
  case INPAR::MAT::m_scatra:
    MatScaTra(material,densn,densnp,densam,visc,iquad);
    break;
  default:
    dserror("Material type %i is not supported",material->MaterialType());
    break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                         wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad   //!< id of current gauss point (default = -1)
  )
{

  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  // get constant diffusivity
  diffmanager_->SetIsotropicDiff(actmat->Diffusivity());

  return;
} // ReynoldsEleCalc<distype>::MatScaTra

/*-----------------------------------------------------------------------------*
 | compute rhs containing bodyforce                                wirtz 10/15 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::GetRhsInt(
  double&      rhsint,  //!< rhs containing bodyforce at Gauss point
  const double densnp  //!< density at t_(n+1)
  )
{
  // compute rhs containing bodyforce (divided by specific heat capacity) and,
  // for temperature equation, the time derivative of thermodynamic pressure,
  // if not constant, and for temperature equation of a reactive
  // equation system, the reaction-rate term
  rhsint = bodyforce_.Dot(funct_);

  return;
} // GetRhsInt

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix               wirtz 10/15 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::CalcMatDiff(
  Epetra_SerialDenseMatrix&     emat,
  const double                  timefacfac
  )
{
  // diffusive term
  const double fac_diffus = timefacfac * diffmanager_->GetIsotropicDiff();
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int ui=0; ui<nen_; ++ui)
    {
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf,ui,vi);
      emat(vi,ui) += fac_diffus*laplawf;
    }
  }
  return;
}

/*-------------------------------------------------------------------------------------- *
 |  standard Galerkin transient, old part of rhs and source term             wirtz 10/15 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::CalcRHSHistAndSource(
  Epetra_SerialDenseVector&     erhs,
  const double                  fac,
  const double                  rhsint
  )
{
  double vrhs = fac*rhsint;
  for (int vi=0; vi<nen_; ++vi)
  {
    erhs[vi] += vrhs*funct_(vi);
  }

  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side    wirtz 10/15 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::CalcRHSDiff(
  Epetra_SerialDenseVector&     erhs,
  const double                  rhsfac
  )
{
  const LINALG::Matrix<nsd_,1>& gradpre = reynoldsvarmanager_->GradPre();

  double vrhs = rhsfac*diffmanager_->GetIsotropicDiff();

  for (int vi=0; vi<nen_; ++vi)
  {
    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf,gradpre,vi);
    erhs[vi] -= vrhs*laplawf;
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point   wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
double DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::EvalShapeFuncAndDerivsAtIntPoint(
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_>& intpoints,  ///< integration points
  const int                                           iquad       ///< id of current Gauss point
  )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_ele_;idim++)
    xsi_(idim) = gpcoord[idim];

  const double det = EvalShapeFuncAndDerivsInParameterSpace();

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

  // compute global spatial derivatives
  derxy_.Multiply(xij_,deriv_);

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad]*det;

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

} //ReynoldsImpl::EvalShapeFuncAndDerivsAtIntPoint

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives in parameter space wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
double DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::EvalShapeFuncAndDerivsInParameterSpace()
{
  double det=0.0;

  if(nsd_==nsd_ele_) //standard case
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);


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
    det = xij_.Invert(xjm_);
  }
  else //element dimension is smaller than problem dimension -> mannifold
  {
    static LINALG::Matrix<nsd_ele_,nen_> deriv_red;

    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_,funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_red);

    //! metric tensor at integration point
    static LINALG::Matrix<nsd_ele_,nsd_ele_>  metrictensor;
    static LINALG::Matrix<nsd_,        1> normalvec;

    // the metric tensor and the area of an infinitesimal surface/line element
    // optional: get unit normal at integration point as well
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype,nsd_>(xyze_,deriv_red,metrictensor,det,&normalvec);

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

    //transform the derivatives and Jacobians to the higher dimensional coordinates(problem dimension)
    static LINALG::Matrix<nsd_ele_,nsd_> xjm_red;
    xjm_red.MultiplyNT(deriv_red,xyze_);

    for(int i=0;i<nsd_;i++)
    {
      for(int j=0;j<nsd_ele_;j++)
        xjm_(j,i)=xjm_red(j,i);
      xjm_(nsd_ele_,i)=normalvec(i,0);
    }

    for(int i=0;i<nen_;i++)
    {
      for(int j=0;j<nsd_ele_;j++)
        deriv_(j,i)=deriv_red(j,i);
      deriv_(nsd_ele_,i)=0.0;
    }

    //special case: 1D element embedded in 3D problem
    if ( nsd_ele_ == 1 and nsd_ == 3 )
    {
      // compute second unit normal
      const double normalvec2_0 = xjm_red(0,1)*normalvec(2,0)-normalvec(1,0)*xjm_red(0,2);
      const double normalvec2_1 = xjm_red(0,2)*normalvec(0,0)-normalvec(2,0)*xjm_red(0,0);
      const double normalvec2_2 = xjm_red(0,0)*normalvec(1,0)-normalvec(0,0)*xjm_red(0,1);

      // norm
      const double norm2 = std::sqrt(normalvec2_0*normalvec2_0+normalvec2_1*normalvec2_1+normalvec2_2*normalvec2_2);

      xjm_(2,0) = normalvec2_0/norm2;
      xjm_(2,1) = normalvec2_1/norm2;
      xjm_(2,2) = normalvec2_2/norm2;

      for(int i=0;i<nen_;i++)
        deriv_(2,i)=0.0;
    }

    xij_.Invert(xjm_);
  }

  return det;
}

/*-----------------------------------------------------------------------*
  |  get the body force  (private)                           wirtz 10/15 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::BodyForce(
  const DRT::Element*    ele
  )
{
  std::vector<DRT::Condition*> myneumcond;

  // check whether all nodes have a unique Neumann condition
  switch(nsd_ele_)
  {
  case 3:
    DRT::UTILS::FindElementConditions(ele, "VolumeNeumann", myneumcond);
    break;
  case 2:
    DRT::UTILS::FindElementConditions(ele, "SurfaceNeumann", myneumcond);
    break;
  case 1:
    DRT::UTILS::FindElementConditions(ele, "LineNeumann", myneumcond);
    break;
  default:
    dserror("Illegal number of spatial dimensions: %d",nsd_ele_);
    break;
  }

  if (myneumcond.size()>1)
    dserror("More than one Neumann condition on one node!");

  if (myneumcond.size()==1)
  {

    // (SPATIAL) FUNCTION BUSINESS
    const std::vector<int>* funct = myneumcond[0]->Get<std::vector<int> >("funct");

    // check for potential time curve
    const std::vector<int>* curve  = myneumcond[0]->Get<std::vector<int> >("curve");
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0];

    // initialization of time-curve factor
    double curvefac(0.0);
    const double time = reynoldspara_->Time();

    // compute potential time curve or set time-curve factor to one
    if (curvenum >= 0)
    {
      // time factor (negative time indicating error)
      if (time >= 0.0)
        curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      else dserror("Negative time in bodyforce calculation: time = %f",time);
    }
    else curvefac = 1.0;

    // get values and switches from the condition
    const std::vector<int>*    onoff = myneumcond[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = myneumcond[0]->Get<std::vector<double> >("val"  );

    // function evaluation
    const int functnum = (funct) ? (*funct)[0] : -1;
    for (int jnode = 0; jnode < nen_; jnode++)
    {
      const double functfac =
          (functnum > 0) ?
              DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(0,
                  (ele->Nodes()[jnode])->X(), time,
                  NULL) :
              1.0;
      bodyforce_(jnode) = (*onoff)[0]*(*val)[0]*curvefac*functfac;
    }
  }
  else
  {
      // no bodyforce
      bodyforce_.Clear();
  }

  return;

} //ReynoldsEleCalc::BodyForce

/*----------------------------------------------------------------------*
 | evaluate service routine                                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::EvaluateService(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la,
    Epetra_SerialDenseMatrix&     elemat1_epetra,
    Epetra_SerialDenseMatrix&     elemat2_epetra,
    Epetra_SerialDenseVector&     elevec1_epetra,
    Epetra_SerialDenseVector&     elevec2_epetra,
    Epetra_SerialDenseVector&     elevec3_epetra
    )
{
  // setup
  if(SetupCalc(ele,discretization) == -1)
    return 0;

  // check for the action parameter
  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");

  // evaluate action
  EvaluateAction(
      ele,
      params,
      discretization,
      action,
      la,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );

  return 0;
}

/*----------------------------------------------------------------------*
 | evaluate action                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
int DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::EvaluateAction(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    const SCATRA::Action&         action,
    DRT::Element::LocationArray&  la,
    Epetra_SerialDenseMatrix&     elemat1_epetra,
    Epetra_SerialDenseMatrix&     elemat2_epetra,
    Epetra_SerialDenseVector&     elevec1_epetra,
    Epetra_SerialDenseVector&     elevec2_epetra,
    Epetra_SerialDenseVector&     elevec3_epetra
    )
{
  //(for now) only first dof set considered
  const std::vector<int> &    lm = la[0].lm_;
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::calc_error:
  {
    // check if length suffices
    if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

    // need current solution
    Teuchos::RCP<const Epetra_Vector> prenp = discretization.GetState("phinp");
    if (prenp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_,1> >(*prenp,eprenp_,lm);

    CalErrorComparedToAnalytSolution(
      ele,
      params,
      elevec1_epetra);

    break;
  }

  default:
  {
    dserror("Not acting on this action. Forgot implementation?");
    break;
  }
  } // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
  |  calculate error compared to analytical solution        wirtz 10/15 |
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::CalErrorComparedToAnalytSolution(
  const DRT::Element*                   ele,
  Teuchos::ParameterList&               params,
  Epetra_SerialDenseVector&             errors
  )
{
  if (DRT::INPUT::get<SCATRA::Action>(params,"action") != SCATRA::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // set constants for analytical solution
  const double t = reynoldspara_->Time();

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcError errortype = DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag");
  switch(errortype)
  {

  case INPAR::SCATRA::calcerror_byfunction:
  {
    const int errorfunctno = params.get<int>("error function number");

    // analytical solution
    double  pre_exact(0.0);
    double  deltapre(0.0);
    //! spatial gradient of current scalar value
    LINALG::Matrix<nsd_,1>  gradpre(true);
    LINALG::Matrix<nsd_,1>  gradpre_exact(true);
    LINALG::Matrix<nsd_,1>  deltagradpre(true);

    // start loop over integration points
    for (int iquad=0;iquad<intpoints.IP().nquad;iquad++)
    {
      const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // get coordinates at integration point
      // gp reference coordinates
      LINALG::Matrix<nsd_,1> xyzint(true);
      xyzint.Multiply(xyze_,funct_);

      // function evaluation requires a 3D position vector!!
      double position[3]={0.0,0.0,0.0};

      for (int dim=0; dim<nsd_; ++dim)
        position[dim] = xyzint(dim);

      // scalar at integration point at time step n+1
      const double prenp = funct_.Dot(eprenp_);
      // spatial gradient of current scalar value
      gradpre.Multiply(derxy_, eprenp_);

      pre_exact = DRT::Problem::Instance()->Funct(errorfunctno - 1).Evaluate(0,
          position, t, NULL);

      std::vector<std::vector<double> > gradpre_exact_vec =
          DRT::Problem::Instance()->Funct(errorfunctno - 1).FctDer(0, position,
              t, NULL);

      if (gradpre_exact_vec.size())
      {
        if (nsd_ == nsd_ele_)
          for (int dim = 0; dim < nsd_; ++dim)
            gradpre_exact(dim) = gradpre_exact_vec[0][dim];
        else
        {
          // Todo: calc gradpre correctly
//          std::cout
//              << "Warning: Gradient of analytical solution cannot be evaluated correctly for transport on curved surfaces!"
//              << std::endl;
          gradpre_exact.Clear();
        }
      }
      else
      {
        std::cout
            << "Warning: Gradient of analytical solution was not evaluated!"
            << std::endl;
        gradpre_exact.Clear();
      }

      // error at gauss point
      deltapre = prenp - pre_exact;
      deltagradpre.Update(1.0, gradpre, -1.0, gradpre_exact);

      // 0: delta scalar for L2-error norm
      // 1: delta scalar for H1-error norm
      // 2: analytical scalar for L2 norm
      // 3: analytical scalar for H1 norm

      // the error for the L2 and H1 norms are evaluated at the Gauss point

      // integrate delta scalar for L2-error norm
      errors(0) += deltapre * deltapre * fac;
      // integrate delta scalar for H1-error norm
      errors(1) += deltapre * deltapre * fac;
      // integrate analytical scalar for L2 norm
      errors(2) += pre_exact * pre_exact * fac;
      // integrate analytical scalar for H1 norm
      errors(3) += pre_exact * pre_exact * fac;

      // integrate delta scalar derivative for H1-error norm
      errors(1) += deltagradpre.Dot(deltagradpre) * fac;
      // integrate analytical scalar derivative for H1 norm
      errors(3) += gradpre_exact.Dot(gradpre_exact) * fac;
    } // loop over integration points

  }
  break;
  default: dserror("Unknown analytical solution!"); break;
  } //switch(errortype)

  return;
} // DRT::ELEMENTS::ReynoldsEleCalc<distype,probdim>::CalErrorComparedToAnalytSolution

// template classes

// 1D elements
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::line2,1>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::line2,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::line2,3>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::line3,1>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::line3,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::tri3,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::tri3,3>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::tri6,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::tri6,3>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::quad4,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::quad4,3>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::quad8,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::quad8,3>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::quad9,2>;
template class DRT::ELEMENTS::ReynoldsEleCalc<DRT::Element::quad9,3>;

