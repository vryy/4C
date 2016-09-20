/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_calc.cpp

 \brief implementation of the evaluation routines of the porofluidmultiphase element

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#include "porofluidmultiphase_ele_calc.H"
#include "porofluid_phasemanager.H"
#include "porofluid_variablemanager.H"
#include "porofluidmultiphase_ele_parameter.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/structporo.H"

#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::PoroFluidMultiPhaseEleCalc(
    const int numdofpernode,
    const std::string& disname):
    ele_(NULL),
    numdofpernode_(numdofpernode),
    para_(DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(disname)),            // standard parameter list
    xsi_(true),
    xyze0_(true),
    funct_(true),
    deriv_(true),
    deriv2_(true),
    derxy_(true),
    derxy2_(true),
    xjm_(true),
    xij_(true),
    det_(0.0),
    J_(0.0),
    structmat_(Teuchos::null),
    phasemanager_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*
 | singleton access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>* DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Instance(
    const int numdofpernode,
    const std::string& disname,
    const PoroFluidMultiPhaseEleCalc* delete_me
    )
{
  static std::map<std::pair<std::string,int>,PoroFluidMultiPhaseEleCalc<distype>* > instances;

  std::pair<std::string,int> key(disname,numdofpernode);

  if(delete_me == NULL)
  {
    if(instances.find(key) == instances.end())
      instances[key] = new PoroFluidMultiPhaseEleCalc<distype>(numdofpernode,disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for( typename std::map<std::pair<std::string,int>,PoroFluidMultiPhaseEleCalc<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}

/*----------------------------------------------------------------------*
 | singleton destruction                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Done()
{
  // delete instance
  Instance(0,"",this);

  return;
}

/*----------------------------------------------------------------------*
 | setup element evaluation                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::SetupCalc(
    DRT::Element*               ele,
    DRT::Discretization&        discretization,
    const POROFLUIDMULTIPHASE::Action& action
    )
{
  // get element coordinates
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze0_);

  // set current coordinates to initial coordinates
  // the displacements will be added later in ExtractElementAndNodeValues() for the moving mesh case
  xyze_=xyze0_;

  // set element
  ele_ = ele;

  // get poro structure material
  GetStructMaterial(ele);

  if(phasemanager_ == Teuchos::null)
  {
    //build the phase manager
    phasemanager_ =
        DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::CreatePhaseManager(
            *para_,
            ele->Material()->MaterialType(),
            action,
            numdofpernode_);
  }
  else
  {
    // just unwrap the manager and rebuild the extensions (we save one setup call)
    Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface > core = phasemanager_->UnWrap();
    core->ClearGPState();
    phasemanager_ =
        DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::WrapPhaseManager(
            *para_,
            ele->Material()->MaterialType(),
            action,
            core);
  }
  // setup the manager
  phasemanager_->Setup(*ele->Material(0));

  //just to be sure: rebuild the phase manager
  variablemanager_ = Teuchos::null;
  variablemanager_ =
      DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<nsd_,nen_>::CreateVariableManager(
          *para_,
          action,
          numdofpernode_);

  return 0;
}

/*----------------------------------------------------------------------*
 *    Get structure material                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GetStructMaterial(DRT::Element* ele)
{
  //get solid material
  {
    //access second material in structure element
    if (ele->NumMaterial() > 1)
    {
      structmat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(ele->Material(1));
      if(structmat_->MaterialType() != INPAR::MAT::m_structporo)
        dserror("invalid structure material for poroelasticity");
    }
    else
      dserror("no second material defined for element %i",ele_->Id());
  }

  return;
}

/*----------------------------------------------------------------------*
 * Evaluate element                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Evaluate(
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

  if(SetupCalc(ele,discretization,POROFLUIDMULTIPHASE::calc_mat_and_rhs) == -1)
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
 | extract element based or nodal values                     vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::ExtractElementAndNodeValues(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    DRT::Element::LocationArray&  la
)
{
  variablemanager_->ExtractElementAndNodeValues(*ele,discretization,la,xyze_);
  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::Sysmat(
  DRT::Element*                         ele,        ///< the element whose matrix is calculated
  Epetra_SerialDenseMatrix&             emat,       ///< element matrix to calculate
  Epetra_SerialDenseVector&             erhs       ///< element rhs to calculate
  )
{

  // prepare gauss point evaluation
  PrepareGaussPointLoop(ele);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(POROFLUIDMULTIPHASE::ELEUTILS::DisTypeToOptGaussRule<distype>::rule);

  // start loop over gauss points
  GaussPointLoop(intpoints,ele,emat,erhs);

  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::PrepareGaussPointLoop(DRT::Element* ele)
{

  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 vuong 08/16|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GaussPointLoop(
    const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs       //!< element rhs to calculate
)
{
  // start the loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    variablemanager_->EvaluateGPVariables(funct_,derxy_);

    //diffusion tensor
    static LINALG::Matrix<nsd_,nsd_> difftensor;
    difftensor.Clear();

    // set the current gauss point state in the manager
    phasemanager_->EvaluateGPState(*ele->Material(0),*structmat_,J_,*variablemanager_);

    // get porosity form structure material
    const double porosity = phasemanager_->Porosity();


    // loop over phases
    for (int k=0;k<numdofpernode_;++k) // deal with a system of transported scalars
    {

      //----------------------------------------------------------------
      // standard Galerkin terms
      //----------------------------------------------------------------

      // stabilization parameter and integration factors
      //const double taufac     = tau[k]*fac;
      const double timefacfac = para_->TimeFac()*fac;
      //const double timetaufac = para_->TimeFac()*taufac;

      //----------------------------------------------------------------
      // 1) element matrix: stationary terms
      //----------------------------------------------------------------

      // calculation of convective element matrix in convective form
//      if(k!=numdofpernode_-1)
//        CalcMatConv(
//            emat,
//            k,
//            *phasemanager_,
//            porosity,
//            timefacfac,
//            conv);

      // add conservative contributions
      if(para_->IsAle())
        if(k!=numdofpernode_-1)
          CalcMatConvCons(
              emat,
              k,
              *phasemanager_,
              timefacfac,
              variablemanager_->DivConVelnp());

      // compute prefactor
      ComputeDiffTensor(
          *phasemanager_,
          k,
          *ele->Material(0),
          *variablemanager_->Phinp(),
          difftensor);

      // calculation of diffusive element matrix
      CalcMatDiff(
          emat,
          k,
          k,
          *phasemanager_,
          *ele->Material(0),
          timefacfac,
          difftensor);

      // add the terms also into the last equation
      if(k!=numdofpernode_-1)
        CalcMatDiff(
            emat,
            k,
            numdofpernode_-1,
            *phasemanager_,
            *ele->Material(0),
            timefacfac,
            difftensor);

      // calculation of reactive element matrix
      if( phasemanager_->IsReactive(k) )
      {
        CalcMatReac(
            emat,
            k,
            k,
            *phasemanager_,
            *ele->Material(0),
            timefacfac,
            true // do scaling with saturation
            );
        // add the terms also into the last equation
        if(k!=numdofpernode_-1)
          CalcMatReac(
              emat,
              k,
              numdofpernode_-1,
              *phasemanager_,
              *ele->Material(0),
              timefacfac,
              false // do NOT scaling with saturation
              );
      }

      //----------------------------------------------------------------
      // 2) element matrix: instationary terms
      //----------------------------------------------------------------


      if (not para_->IsStationary())
      {
        CalcMatMassPressure(
            emat,
            k,
            k,
            *phasemanager_,
            *ele->Material(0),
            fac,
            timefacfac,
            porosity,
            (*variablemanager_->Hist())[k],
            *variablemanager_->Phinp(),
            *variablemanager_->Phidtnp());

        if(k!=numdofpernode_-1)
          CalcMatMassPressure(
              emat,
              k,
              numdofpernode_-1,
              *phasemanager_,
              *ele->Material(0),
              fac,
              timefacfac,
              porosity,
              0.0, // no history value! it has already been included
              *variablemanager_->Phinp(),
              *variablemanager_->Phidtnp());

        CalcMatMassSolidPressure(
            emat,
            k,
            *phasemanager_,
            *ele->Material(0),
            fac,
            timefacfac,
            porosity,
            (*variablemanager_->Hist())[k],
            *variablemanager_->Phinp(),
            *variablemanager_->Phidtnp(),
            k!=numdofpernode_-1);

        if(k!=numdofpernode_-1)
          CalcMatMassSaturation(
              emat,
              k,
              *phasemanager_,
              *ele->Material(0),
              fac,
              porosity,
              *variablemanager_->Phinp());
      }

      //----------------------------------------------------------------
      // 5) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac    = para_->TimeFacRhs() * fac;

      if (not para_->IsStationary())
      {
        CalcRHSMassPressure(
            erhs,
            k,
            k,
            *phasemanager_,
            *ele->Material(0),
            rhsfac,
            fac,
            (*variablemanager_->Hist())[k],
            porosity,
            *variablemanager_->Phinp(),
            *variablemanager_->Phidtnp());

        if(k!=numdofpernode_-1)
          CalcRHSMassPressure(
              erhs,
              k,
              numdofpernode_-1,
              *phasemanager_,
              *ele->Material(0),
              rhsfac,
              fac,
              0.0,
              porosity,
              *variablemanager_->Phinp(),
              *variablemanager_->Phidtnp());

        CalcRHSMassSolidPressure(
            erhs,
            k,
            k,
            *phasemanager_,
            *ele->Material(0),
            rhsfac,
            fac,
            (*variablemanager_->Hist())[k],
            porosity,
            *variablemanager_->Phinp(),
            *variablemanager_->Phidtnp(),
            k!=numdofpernode_-1);

        if(k!=numdofpernode_-1)
          CalcRHSMassSaturation(
              erhs,
              k,
              k,
              *phasemanager_,
              *ele->Material(0),
              rhsfac,
              fac,
              (*variablemanager_->Hist())[k],
              porosity,
              *variablemanager_->Phinp(),
              *variablemanager_->Phidtnp());
      }

      //----------------------------------------------------------------
      // standard Galerkin terms on right hand side
      //----------------------------------------------------------------

      // convective term
//      if(k!=numdofpernode_-1)
//        CalcRHSConv(
//            erhs,
//            k,
//            *phasemanager_,
//            porosity,
//            rhsfac,
//            conv_phi);

      // convective term conservative contribution
      if(para_->IsAle())
        CalcRHSConvCons(
            erhs,
            k,
            *phasemanager_,
            rhsfac,
            variablemanager_->DivConVelnp(),
            k!=numdofpernode_-1);

      // diffusive term
      CalcRHSDiff(
          erhs,
          k,
          k,
          *phasemanager_,
          rhsfac,
          *variablemanager_->GradPhinp(),
          difftensor);
      if(k!=numdofpernode_-1)
        CalcRHSDiff(
            erhs,
            k,
            numdofpernode_-1,
            *phasemanager_,
            rhsfac,
            *variablemanager_->GradPhinp(),
            difftensor);

      // calculation of reactive element matrix
      if( phasemanager_->IsReactive(k) )
      {
        // reactive term
        CalcRHSReac(
            erhs,
            k,
            k,
            *phasemanager_,
            *ele->Material(0),
            rhsfac,
            true //scale with saturation
            );
        if(k!=numdofpernode_-1)
          CalcRHSReac(
              erhs,
              k,
              numdofpernode_-1,
              *phasemanager_,
              *ele->Material(0),
              rhsfac,
              false //scale with saturation
              );
      }

    }// end loop all scalars

    // clear current gauss point data for safety
    phasemanager_->ClearGPState();
  }

  return;

}

/*------------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point     vuong 08/16 |
 *------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvalShapeFuncAndDerivsAtIntPoint(
  const DRT::UTILS::IntPointsAndWeights<nsd_>& intpoints,  ///< integration points
  const int                                    iquad       ///< id of current Gauss point
  )
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim=0;idim<nsd_;idim++)
    xsi_(idim) = gpcoord[idim];

  det_ = EvalShapeFuncAndDerivsInParameterSpace();

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", ele_->Id(), det_);

  // compute global spatial derivatives
  derxy_.Multiply(xij_,deriv_);

  // compute second global derivatives (if needed)
  if (use2ndderiv_)
  {
    // get global second derivatives
    DRT::UTILS::gder2<distype,nen_>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else
    derxy2_.Clear();

  //------------------------get determinant of Jacobian dX / ds
  // transposed jacobian "dX/ds"
  LINALG::Matrix<nsd_,nsd_> xjm0;
  xjm0.MultiplyNT(deriv_,xyze0_);

  // inverse of transposed jacobian "ds/dX"
  const double det0= xjm0.Determinant();

  // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  J_ = det_/det0;

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad]*det_;

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

} //ScaTraImpl::EvalShapeFuncAndDerivsAtIntPoint

/*--------------------------------------------------------------------------*
 | evaluate shape functions and derivatives in parameter space   vuong 08/16 |
 *--------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvalShapeFuncAndDerivsInParameterSpace()
{
  double det=0.0;

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  if (use2ndderiv_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
  }

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

  return det;
}

/*-----------------------------------------------------------------------------*
 |  calculate divergence of vector field (e.g., velocity)  (private) vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GetDivergence(
  double&                          vdiv,
  const LINALG::Matrix<nsd_,nen_>& evel)
{
  LINALG::Matrix<nsd_,nsd_> vderxy;
  vderxy.MultiplyNT(evel,derxy_);

  vdiv = 0.0;
  // compute vel x,x  + vel y,y +  vel z,z at integration point
  for (int j = 0; j<nsd_; ++j)
  {
    vdiv += vderxy(j,j);
  }
  return;
} // ScaTraEleCalc<distype>::GetDivergence


///*-----------------------------------------------------------------------------*
// | compute rhs containing bodyforce                                 vuong 08/16 |
// *-----------------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::GetRhsInt(
//  double&      rhsint,  //!< rhs containing bodyforce at Gauss point
//  const double densnp,  //!< density at t_(n+1)
//  const int    k        //!< index of current scalar
//  )
//{
//  // compute rhs containing bodyforce (divided by specific heat capacity) and,
//  // for temperature equation, the time derivative of thermodynamic pressure,
//  // if not constant, and for temperature equation of a reactive
//  // equation system, the reaction-rate term
//  rhsint = bodyforce_[k].Dot(funct_);
//
//  return;
//} // GetRhsInt

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix in convective form     vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatConv(
    Epetra_SerialDenseMatrix&          emat,
    const int                          k,
    const POROFLUIDMANAGER::PhaseManagerInterface&       phasemanager,
    const double                       porosity,
    const double                       timefacfac,
    const LINALG::Matrix<nen_,1>&      conv
  )
{

  // convective term in convective form
  /*
       /                               \
      |                                 |
      | prefac * v * nabla * Dphi  , q  |
      |                                 |
       \                               /
  */
  const double prefac = timefacfac*porosity;

  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = prefac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        const int fui = ui*numdofpernode_+idof;

        emat(fvi,fui) += v*conv(ui)*phasemanager.SaturationDeriv(k,idof);
      }
  }
  return;
} // ScaTraEleCalc<distype>::CalcMatConv


/*------------------------------------------------------------------------------------------*
 |  calculation of convective element matrix: add conservative contributions     vuong 08/16 |
 *------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatConvCons(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const double                  timefacfac,
  const double                  vdiv
  )
{

  const double consfac = timefacfac*vdiv;
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = consfac*funct_(vi);
    const int fvi = vi*numdofpernode_+k;

    for (int ui=0; ui<nen_; ++ui)
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        const int fui = ui*numdofpernode_+idof;

        emat(fvi,fui) += v*phasemanager.SaturationDeriv(k,idof)*funct_(ui);
      }
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                vuong 08/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatDiff(
  Epetra_SerialDenseMatrix&         emat,
  const int                         curphase,
  const int                         phasetoadd,
  const POROFLUIDMANAGER::PhaseManagerInterface&      phasemanager,
  const MAT::Material&              material,
  const double                      timefacfac,
  const LINALG::Matrix<nsd_,nsd_>&  difftensor
  )
{
  static LINALG::Matrix<nsd_,nen_> diffflux(true);
  diffflux.Multiply(difftensor,derxy_);
  // diffusive term
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+phasetoadd;

    for (int ui=0; ui<nen_; ++ui)
    {
      double laplawf(0.0);
      for (int j = 0; j<nsd_; j++)
        laplawf += derxy_(j, vi)*diffflux(j, ui);
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        const int fui = ui*numdofpernode_+idof;
        emat(fvi,fui) += timefacfac*laplawf*phasemanager.PressureDeriv(curphase,idof);
      }
    }
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of reactive element matrix                vuong 08/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatReac(
  Epetra_SerialDenseMatrix&          emat,
  const int                          curphase,
  const int                          phasetoadd,
  const POROFLUIDMANAGER::PhaseManagerInterface&       phasemanager,
  const MAT::Material&               material,
  const double                       timefacfac,
  bool                               scalewithsaturation
  )
{
  // TODO a constant density is assumed here
  double scaledtimefacfac =timefacfac/phasemanager.Density(material,curphase);

  if(scalewithsaturation)
  {
    // saturation
    scaledtimefacfac *= phasemanager.Saturation(curphase);
  }

  //----------------------------------------------------------------
  // reaction terms
  //----------------------------------------------------------------
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = scaledtimefacfac*funct_(vi);
    const int fvi = vi*numdofpernode_+phasetoadd;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double vfunct = v*funct_(ui);
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        const int fui = ui*numdofpernode_+idof;

        emat(fvi,fui) += vfunct*phasemanager.ReacDeriv(curphase,idof);
      }
    }
  }

  //----------------------------------------------------------------
  // linearization of saturation w.r.t. dof
  //----------------------------------------------------------------
  if(scalewithsaturation)
  {
    double facfacreac = timefacfac*phasemanager.ReacTerm(curphase);

    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = facfacreac*funct_(vi);
      const int fvi = vi*numdofpernode_+phasetoadd;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double vfunct = v*funct_(ui);
        for (int idof=0; idof<numdofpernode_; ++idof)
        {
          const int fui = ui*numdofpernode_+idof;

          emat(fvi,fui) +=
              vfunct*phasemanager.SaturationDeriv(curphase,idof);
        }
      }
    }
  }

  return;
}


/*------------------------------------------------------------------- *
 |  calculation of mass element matrix                    vuong 08/16  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatMassPressure(
  Epetra_SerialDenseMatrix&     emat,
  const int                     curphase,           //!< index of current phase
  const int                     phasetoadd,           //!< index of current phase
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const MAT::Material&          material,
  const double                  fac,
  const double                  timefacfac,
  const double                  porosity,
  const double                  hist,
  const std::vector<double>&    phinp,
  const std::vector<double>&    phidtnp
  )
{
  // saturation
  const double saturation = phasemanager.Saturation(curphase);

  // bulkmodulus of phase
  const double bulkmodulus = phasemanager.Bulkmodulus(material,curphase);

  // pre factor
  const double facfacmass = fac*porosity*saturation/bulkmodulus;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = facfacmass*funct_(vi);
    const int fvi = vi*numdofpernode_+phasetoadd;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double vfunct = v*funct_(ui);
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        const int fui = ui*numdofpernode_+idof;

        emat(fvi,fui) += vfunct*phasemanager.PressureDeriv(curphase,idof);
      }
    }
  }

  //----------------------------------------------------------------
  // linearization of saturation w.r.t. dof
  //----------------------------------------------------------------
  {
    double facfacmass2 = fac*phasemanager.PressureDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

    for (int idof=0; idof<numdofpernode_; ++idof)
    {
      if(idof!=phasetoadd)
        facfacmass2 += timefacfac*phasemanager.PressureDeriv(curphase,idof)*phidtnp[idof];
    }

    facfacmass2 *= porosity/bulkmodulus;

    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = facfacmass2*funct_(vi);
      const int fvi = vi*numdofpernode_+phasetoadd;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double vfunct = v*funct_(ui);
        for (int idof=0; idof<numdofpernode_; ++idof)
        {
          const int fui = ui*numdofpernode_+idof;

          emat(fvi,fui) +=
              vfunct*phasemanager.SaturationDeriv(curphase,idof);
        }
      }
    }
  }


  return;
}

/*------------------------------------------------------------------- *
 |  calculation of mass element matrix                    vuong 08/16  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatMassSolidPressure(
  Epetra_SerialDenseMatrix&     emat,
  const int                     curphase,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const MAT::Material&          material,
  const double                  fac,
  const double                  timefacfac,
  const double                  porosity,
  const double                  hist,
  const std::vector<double>&    phinp,
  const std::vector<double>&    phidtnp,
  bool                          scalewithsaturation
  )
{
  double scale =1.0;

  if(scalewithsaturation)
  {
    // saturation
    scale = phasemanager.Saturation(curphase);
  }

  //  get inverse bulkmodulus (=compressiblity)
  // TODO linearization of bulkmodulus
  const double invsolidbulkmodulus = structmat_->InvBulkmodulus();

  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  {
    const double facfacmass = fac*(1.0-porosity)*scale*invsolidbulkmodulus;
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = facfacmass*funct_(vi);
      const int fvi = vi*numdofpernode_+curphase;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double vfunct = v*funct_(ui);
        for (int idof=0; idof<numdofpernode_; ++idof)
        {
          const int fui = ui*numdofpernode_+idof;

          emat(fvi,fui) += vfunct*phasemanager.SolidPressureDeriv(idof);
        }
      }
    }
  }

  //----------------------------------------------------------------
  // linearization of saturation w.r.t. dof
  //----------------------------------------------------------------
  if(scalewithsaturation)
  {
    double facfacmass2 = fac*phasemanager.SolidPressureDeriv(curphase)*(phinp[curphase]-hist);

    for (int idof=0; idof<numdofpernode_; ++idof)
    {
      if(idof!=curphase)
        facfacmass2 += timefacfac*phasemanager.SolidPressureDeriv(idof)*phidtnp[idof];
    }

    facfacmass2 *= (1.0-porosity)*invsolidbulkmodulus;

    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = facfacmass2*funct_(vi);
      const int fvi = vi*numdofpernode_+curphase;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double vfunct = v*funct_(ui);
        for (int idof=0; idof<numdofpernode_; ++idof)
        {
          const int fui = ui*numdofpernode_+idof;

          emat(fvi,fui) +=
              vfunct*phasemanager.SaturationDeriv(curphase,idof);
        }
      }
    }
  }

  //----------------------------------------------------------------
  // linearization of solid pressure derivative w.r.t. dof
  //----------------------------------------------------------------
  {
    double facfacmass3 = scale*(1.0-porosity)*invsolidbulkmodulus;

    std::vector<double> val(numdofpernode_,0.0);

    for (int idof=0; idof<numdofpernode_; ++idof)
    {
      if(idof==curphase)
        for (int jdof=0; jdof<numdofpernode_; ++jdof)
          val[jdof]+= fac*phasemanager.SolidPressureDerivDeriv(idof,jdof)*(phinp[idof]-hist);
      else
        for (int jdof=0; jdof<numdofpernode_; ++jdof)
          val[jdof]+= timefacfac*phasemanager.SolidPressureDerivDeriv(idof,jdof)*(phidtnp[idof]);
    }

    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = facfacmass3*funct_(vi);
      const int fvi = vi*numdofpernode_+curphase;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double vfunct = v*funct_(ui);
        for (int idof=0; idof<numdofpernode_; ++idof)
        {
          const int fui = ui*numdofpernode_+idof;

          emat(fvi,fui) +=vfunct*val[idof];
        }
      }
    }
  }

}

/*------------------------------------------------------------------- *
 |  calculation of mass element matrix                    vuong 08/16  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcMatMassSaturation(
  Epetra_SerialDenseMatrix&     emat,
  const int                     curphase,           //!< index of current scalar
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const MAT::Material&          material,
  const double                  fac,
  const double                  porosity,
  const std::vector<double>&    phinp
  )
{
  const double facfacmass = fac*porosity;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (int vi=0; vi<nen_; ++vi)
  {
    const double v = facfacmass*funct_(vi);
    const int fvi = vi*numdofpernode_+curphase;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double vfunct = v*funct_(ui);
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        const int fui = ui*numdofpernode_+idof;

        emat(fvi,fui) += vfunct*phasemanager.SaturationDeriv(curphase,idof);
      }
    }
  }
}

/*------------------------------------------------------------------- *
 |  calculation of linearized mass rhs vector              vuong 08/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSMassPressure(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,
  const int                     phasetoadd,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const MAT::Material&          material,
  const double                  rhsfac,
  const double                  fac,
  const double                  hist,
  const double                  porosity,
  const std::vector<double>&    phinp,
  const std::vector<double>&    phidtnp
  )
{
  // saturation
  const double saturation = phasemanager.Saturation(curphase);

  // bulk modulus of phase
  const double bulkmodulus = phasemanager.Bulkmodulus(material,curphase);

  double vtrans = 0.0;

  //TODO check for Genalpha
  // compute scalar at integration point
  vtrans = fac*phasemanager.PressureDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);
  for (int idof=0; idof<numdofpernode_; ++idof)
    if(idof!=phasetoadd)
      vtrans += rhsfac*phasemanager.PressureDeriv(curphase,idof)*phidtnp[idof];

  vtrans *= porosity*saturation/bulkmodulus;

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+phasetoadd;

    erhs[fvi] -= vtrans*funct_(vi);
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of linearized mass rhs vector              vuong 08/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSMassSolidPressure(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,
  const int                     phasetoadd,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const MAT::Material&          material,
  const double                  rhsfac,
  const double                  fac,
  const double                  hist,
  const double                  porosity,
  const std::vector<double>&    phinp,
  const std::vector<double>&    phidtnp,
  bool                          scalewithsaturation
  )
{
  double scale=1.0;

  if(scalewithsaturation)
  {
    // saturation
    scale = phasemanager.Saturation(curphase);
  }

  //  get inverse bulkmodulus (=compressiblity)
  const double invsolidbulkmodulus = structmat_->InvBulkmodulus();

  //TODO check genalpha
  // compute scalar at integration point
  double vtrans = fac*phasemanager.SolidPressureDeriv(phasetoadd)*(phinp[phasetoadd]-hist);

  for (int idof=0; idof<numdofpernode_; ++idof)
  {
    if(idof!=phasetoadd)
    {
      vtrans += rhsfac*phasemanager.SolidPressureDeriv(idof)*phidtnp[idof];
    }
  }

  vtrans *= (1.0-porosity)*scale*invsolidbulkmodulus;

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+phasetoadd;

    erhs[fvi] -= vtrans*funct_(vi);
  }
}

/*------------------------------------------------------------------- *
 |  calculation of linearized mass rhs vector              vuong 08/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSMassSaturation(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,
  const int                     phasetoadd,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const MAT::Material&          material,
  const double                  rhsfac,
  const double                  fac,
  const double                  hist,
  const double                  porosity,
  const std::vector<double>&    phinp,
  const std::vector<double>&    phidtnp
  )
{
  //TODO genalpha
  // compute scalar at integration point
  double vtrans = fac*phasemanager.SaturationDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

  for (int idof=0; idof<numdofpernode_; ++idof)
    if(phasetoadd!=idof)
      vtrans += rhsfac*phasemanager.SaturationDeriv(curphase,idof)*phidtnp[idof];

  vtrans *= porosity;

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+phasetoadd;

    erhs[fvi] -= vtrans*funct_(vi);
  }

  return;
}

/*-------------------------------------------------------------------------------------- *
 |  standard Galerkin transient, old part of rhs and source term              vuong 08/16 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSHistAndSource(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,
  const double                  fac,
  const double                  rhsint
  )
{
  double vrhs = fac*rhsint;
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+curphase;

    erhs[fvi] += vrhs*funct_(vi);
  }

  return;
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin convective term on right hand side    vuong 08/16 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSConv(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const double                  porosity,
  const double                  rhsfac,
  const std::vector<double>&    conv_phi
  )
{
  double conv_sat = 0.0;
  for (int idof=0; idof<numdofpernode_; ++idof)
    conv_sat += rhsfac*porosity*phasemanager.SaturationDeriv(curphase,idof)*conv_phi[idof];
  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+curphase;

    erhs[fvi] -= conv_sat*funct_(vi);
  }

  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin convective term (conservative contribution       |
 |  on right hand side                                      vuong 08/16 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSConvCons(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,
  const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
  const double                  rhsfac,
  const double                  vdiv,
  bool                          scalewithsaturation
  )
{
  double vrhs=0.0;
  if(scalewithsaturation)
    vrhs = rhsfac*phasemanager.Saturation(curphase)*vdiv;
  else
    vrhs = rhsfac*vdiv;

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+curphase;

    erhs[fvi] -= vrhs*funct_(vi);
  }

  return;
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     vuong 08/16 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSDiff(
  Epetra_SerialDenseVector&     erhs,
  const int                     curphase,           //!< index of current phase
  const int                     phasetoadd,           //!< index of current phase
  const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
  const double                  rhsfac,
  const std::vector<LINALG::Matrix<nsd_,1> >& gradphi,
  const LINALG::Matrix<nsd_,nsd_>&     difftensor       //!< pre factor of diffusive term
  )
{
  // current pressure gradient
  LINALG::Matrix<nsd_,1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof=0; idof<numdofpernode_; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase,idof),gradphi[idof],1.0);

  // diffusive flux
  static LINALG::Matrix<nsd_,1> diffflux(true);
  diffflux.Multiply(difftensor,gradpres);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+phasetoadd;

    // laplacian in weak form
    double laplawf(0.0);
    for (int j = 0; j<nsd_; j++)
      laplawf += derxy_(j,vi)*diffflux(j);
    erhs[fvi] -= rhsfac*laplawf;
  }

  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin reaction term on right hand side     vuong 08/16 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcRHSReac(
  Epetra_SerialDenseVector&            erhs,
  const int                            curphase,
  const int                            phasetoadd,
  const POROFLUIDMANAGER::PhaseManagerInterface&         phasemanager,
  const MAT::Material&                 material,
  const double                         rhsfac,
  bool                                 scalewithsaturation
  )
{
  // TODO a constant density is assumed here
  double scale = 1.0/phasemanager.Density(material,curphase);

  if(scalewithsaturation)
  {
    // saturation
    scale *= phasemanager.Saturation(curphase);
  }

  double vrhs = scale*rhsfac*phasemanager.ReacTerm(curphase);

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi = vi*numdofpernode_+phasetoadd;
    erhs[fvi] -= vrhs*funct_(vi);
  }

  return;
}


/*------------------------------------------------------------------- *
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::ComputeDiffTensor(
    const POROFLUIDMANAGER::PhaseManagerInterface&  phasemanager,
    const int                     phase,
    const MAT::Material&          material,
    const std::vector<double>&    phinp,
    LINALG::Matrix<nsd_,nsd_>&    difftensor) const
{

  if(material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
     material.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
    dserror("only poro multiphase material valid");

  const MAT::FluidPoroMultiPhase& multiphasemat =
      static_cast<const MAT::FluidPoroMultiPhase&>(material);

  // TODO only isotropic, constant permeability for now
  difftensor.Clear();
  const double permeability=multiphasemat.Permeability();
  for (int i=0;i<nsd_;i++)
    difftensor(i,i) = permeability;

  // relative diffusivity (permeability)
  const double reldiffusivity = phasemanager.EvaluateRelDiffusivity(material,phase);

  difftensor.Scale(reldiffusivity);
  return;
}

/*----------------------------------------------------------------------*
 | evaluate service routine                                  vuong 08/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvaluateService(
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
  // check for the action parameter
  const POROFLUIDMULTIPHASE::Action action =
      DRT::INPUT::get<POROFLUIDMULTIPHASE::Action>(params,"action");

  // setup
  if(SetupCalc(ele,discretization,action) == -1)
    return 0;

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractElementAndNodeValues(ele,params,discretization,la);

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
 | evaluate action                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::EvaluateAction(
    DRT::Element*                          ele,
    Teuchos::ParameterList&                params,
    DRT::Discretization&                   discretization,
    const POROFLUIDMULTIPHASE::Action&     action,
    DRT::Element::LocationArray&           la,
    Epetra_SerialDenseMatrix&              elemat1_epetra,
    Epetra_SerialDenseMatrix&              elemat2_epetra,
    Epetra_SerialDenseVector&              elevec1_epetra,
    Epetra_SerialDenseVector&              elevec2_epetra,
    Epetra_SerialDenseVector&              elevec3_epetra
    )
{

  // determine and evaluate action
  switch(action)
  {
  // calculate true pressures and saturation
  case POROFLUIDMULTIPHASE::calc_pres_and_sat:
  {
    // calculate matrix and rhs
    CalcPressureAndSaturation(
      ele,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra,
      params,
      discretization,
      la
      );
    break;
  }
  // calculate solid pressure
  case POROFLUIDMULTIPHASE::calc_solidpressure:
  {
    // calculate matrix and rhs
    CalcSolidPressure(
      ele,
      elevec1_epetra,
      elevec2_epetra,
      params,
      discretization,
      la
      );
    break;
  }
  // reconstruct velocities
  case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
  {
    // calculate matrix and rhs
    ReconFlux(
      ele,
      elemat1_epetra,
      elemat2_epetra,
      params,
      discretization,
      la
      );
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

/*-----------------------------------------------------------------------------*
 | calculate pressures and saturation                               vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcPressureAndSaturation(
    DRT::Element*                 ele,              //!< current element
    Epetra_SerialDenseVector&     pressure,         //!< pressure vector
    Epetra_SerialDenseVector&     saturation,       //!< saturation vector
    Epetra_SerialDenseVector&     counter,          //!< counter
    Teuchos::ParameterList&       params,           //!< parameter list
    DRT::Discretization&          discretization,   //!< discretization
    DRT::Element::LocationArray&  la                //!< location array
    )
{
  // Note: this is an evaluation at every node (for post processing)
  //       Hence, the shape function matrix funct and the derivative matrix derxy
  //       in such a way that the standard routines can be used

  // dummy derivative matrix (initialized to zero)
  LINALG::Matrix<nsd_,nen_> derxy(true);

  for (int inode = 0; inode < nen_; ++inode)
  {
    // dummy shape function matrix (initialized to zero)
    LINALG::Matrix<nen_,1>   funct(true);
    // set value at current node to 1
    funct(inode) = 1.0;

    // evaluate needed variables
    variablemanager_->EvaluateGPVariables(funct,derxy);

    // set the current node state as GP state in the manager
    phasemanager_->EvaluateGPState(*ele->Material(0),*structmat_,1.0,*variablemanager_);

    // loop over phases
    for (int k=0;k<numdofpernode_;++k) // deal with a system of transported scalars
    {
      // save the pressure value
      pressure[inode*numdofpernode_+k] = phasemanager_->Pressure(k);
      // save the saturation value
      saturation[inode*numdofpernode_+k] = phasemanager_->Saturation(k);
      // mark the evaluated node
      counter[inode*numdofpernode_+k] = 1.0;
    }

    phasemanager_->ClearGPState();
  }
}

/*-----------------------------------------------------------------------------*
 | calculate solid pressure                                         vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::CalcSolidPressure(
    DRT::Element*                 ele,              //!< current element
    Epetra_SerialDenseVector&     solidpressure,    //!< pressure vector
    Epetra_SerialDenseVector&     counter,          //!< counter
    Teuchos::ParameterList&       params,           //!< parameter list
    DRT::Discretization&          discretization,   //!< discretization
    DRT::Element::LocationArray&  la                //!< location array
    )
{
  // Note: this is an evaluation at every node (for post processing)
  //       Hence, the shape function matrix funct and the derivative matrix derxy
  //       in such a way that the standard routines can be used

  // dummy derivative matrix (initialized to zero)
  LINALG::Matrix<nsd_,nen_> derxy(true);

  for (int inode = 0; inode < nen_; ++inode)
  {
    // dummy shape function matrix (initialized to zero)
    LINALG::Matrix<nen_,1>   funct(true);
    // set value at current node to 1
    funct(inode) = 1.0;

    // evaluate needed variables
    variablemanager_->EvaluateGPVariables(funct,derxy);

    // set the current node state as GP state in the manager
    phasemanager_->EvaluateGPState(*ele->Material(0),*structmat_,1.0,*variablemanager_);

    // save the pressure value
    solidpressure[inode] = phasemanager_->SolidPressure();
    // mark the evaluated node
    counter[inode] = 1.0;

    phasemanager_->ClearGPState();
  }
}

/*-----------------------------------------------------------------------------*
 | reconstruct flux field                                          vuong 08/16 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<distype>::ReconFlux(
    DRT::Element*                 ele,              //!< current element
    Epetra_SerialDenseMatrix&     elemat1,          //!< linearizations
    Epetra_SerialDenseMatrix&     elemat2,          //!< rhs
    Teuchos::ParameterList&       params,           //!< parameter list
    DRT::Discretization&          discretization,   //!< discretization
    DRT::Element::LocationArray&  la                //!< location array
    )
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(POROFLUIDMULTIPHASE::ELEUTILS::DisTypeToOptGaussRule<distype>::rule);

   // Loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
  {
    // Get integration factor and evaluate shape func and its derivatives at the integration points.
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,gpid);

    variablemanager_->EvaluateGPVariables(funct_,derxy_);

    //diffusion tensor
    static LINALG::Matrix<nsd_,nsd_> difftensor;
    difftensor.Clear();

    // set the current gauss point state in the manager
    phasemanager_->EvaluateGPState(*ele->Material(0),*structmat_,J_,*variablemanager_);

    // Compute element matrix. For L2-projection
    for (int vi=0; vi<nen_; ++vi)
    {
      const double v = fac*funct_(vi);
      for (int ui=0; ui<nen_; ++ui)
      {
        elemat1(vi,ui) += v*funct_(ui);
      }
    }

    // Loop over degrees of freedom
    for (int k=0; k<numdofpernode_; k++)
    {
      // compute prefactor
      ComputeDiffTensor(
          *phasemanager_,
          k,
          *ele->Material(0),
          *variablemanager_->Phinp(),
          difftensor);

      // current pressure gradient
      LINALG::Matrix<nsd_,1> gradpres(true);
      gradpres.Clear();

      // compute the pressure gradient from the phi gradients
      for (int idof=0; idof<numdofpernode_; ++idof)
      {
        gradpres.Update(phasemanager_->PressureDeriv(k,idof),(*variablemanager_->GradPhinp())[idof],1.0);
      }

      // diffusive flux
      static LINALG::Matrix<nsd_,1> diffflux(true);
      diffflux.Multiply(-1.0,difftensor,gradpres);

      // Compute element vectors. For L2-Projection
      for (int node_i=0;node_i<nen_;node_i++)
      {
        for (int j=0;j<nsd_;j++)
        {
          elemat2(node_i,nsd_*k+j) += funct_(node_i) * fac * diffflux(j);
        }
      }

    } //loop over degrees of freedom

    phasemanager_->ClearGPState();
  } //loop over integration points
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line2>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line2,2>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line2,3>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line3>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line3,2>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tri3>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad4>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad9>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::quad9,3>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::hex8>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::tet10>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::PoroFluidMultiPhaseEleCalc<DRT::Element::nurbs27>;
