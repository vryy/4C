/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_electrode.cpp

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge within isothermal electrodes

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_electrode.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_utils_elch_electrode.H"

#include "../drt_mat/material.H"


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>* DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcElchElectrode* delete_me )
{
  static std::map<std::string,ScaTraEleCalcElchElectrode<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcElchElectrode<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcElchElectrode<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::Done()
{
  // delete singleton
  Instance( 0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::ScaTraEleCalcElchElectrode(const int numdofpernode,const int numscal,const std::string& disname) :
  myelch::ScaTraEleCalcElch(numdofpernode,numscal,disname)
{
  // replace elch diffusion manager by diffusion manager for electrodes
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchElectrode(my::numscal_));

  // replace elch internal variable manager by internal variable manager for electrodes
  my::scatravarmanager_ = Teuchos::rcp(new ScaTraEleInternalVariableManagerElchElectrode<my::nsd_, my::nen_>(my::numscal_,myelch::elchparams_));

  // replace elch utility class by utility class for electrodes
  myelch::utils_ = DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(numdofpernode,numscal,disname);

  return;
}


/*----------------------------------------------------------------------------------------------------*
 | calculate contributions to element matrix and residual (inside loop over all scalars)   fang 02/15 |
 *----------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcMatAndRhs(
    Epetra_SerialDenseMatrix&     emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&     erhs,         //!< element rhs to calculate+
    const int                     k,            //!< index of current scalar
    const double                  fac,          //!< domain-integration factor
    const double                  timefacfac,   //!< domain-integration factor times time-integration factor
    const double                  rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    const double                  taufac,       //!< tau times domain-integration factor
    const double                  timetaufac,   //!< domain-integration factor times tau times time-integration factor
    const double                  rhstaufac,    //!< time-integration factor for rhs times tau times domain-integration factor
    LINALG::Matrix<my::nen_,1>&   tauderpot,    //!< derivatives of stabilization parameter w.r.t. electric potential
    double&                       rhsint        //!< rhs at Gauss point
    )
{
  //----------------------------------------------------------------------
  // 1) element matrix: instationary terms arising from transport equation
  //----------------------------------------------------------------------

  if (not my::scatraparatimint_->IsStationary())
    // 1a) element matrix: standard Galerkin mass term
    my::CalcMatMass(emat,k,fac,1.);

  //--------------------------------------------------------------------
  // 2) element matrix: stationary terms arising from transport equation
  //--------------------------------------------------------------------

  // 2a) element matrix: standard Galerkin diffusive term
  my::CalcMatDiff(emat,k,timefacfac);

  // 2b) element matrix: additional term arising from concentration dependency of diffusion coefficient
  CalcMatDiffCoeffLin(emat,k,timefacfac,VarManager()->GradPhi(k),1.);

  //----------------------------------------------------------------------------
  // 3) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from transport equation
  //----------------------------------------------------------------------------

  // 3a) element rhs: standard Galerkin contributions from non-history part of instationary term if needed
  if(not my::scatraparatimint_->IsStationary())
    my::CalcRHSLinMass(erhs,k,rhsfac,fac,1.,1.);

  // 3b) element rhs: standard Galerkin contributions from rhsint vector (contains body force vector and history vector)
  // need to adapt rhsint vector to time integration scheme first
  my::ComputeRhsInt(rhsint,1.,1.,VarManager()->Hist(k));
  my::CalcRHSHistAndSource(erhs,k,fac,rhsint);

  // 3c) element rhs: standard Galerkin diffusion term
  my::CalcRHSDiff(erhs,k,rhsfac);

  //----------------------------------------------------------------------------
  // 4) element matrix: stationary terms arising from potential equation
  // 5) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from potential equation
  //----------------------------------------------------------------------------
  // see function CalcMatAndRhsOutsideScalarLoop()

  return;
}


/*-----------------------------------------------------------------------------------------------------*
 | calculate contributions to element matrix and residual (outside loop over all scalars)   fang 02/15 |
 *-----------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcMatAndRhsOutsideScalarLoop(
    Epetra_SerialDenseMatrix&   emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&   erhs,         //!< element rhs to calculate
    const double                fac,          //!< domain-integration factor
    const double                timefacfac,   //!< domain-integration factor times time-integration factor
    const double                rhsfac        //!< time-integration factor for rhs times domain-integration factor
  )
{
  //--------------------------------------------------------------------
  // 4) element matrix: stationary terms arising from potential equation
  //--------------------------------------------------------------------

  // element matrix: standard Galerkin terms from potential equation
  CalcMatPotEquDiviOhm(emat,timefacfac,VarManager()->InvF(),VarManager()->GradPot(),1.);

  //----------------------------------------------------------------------------
  // 5) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from potential equation
  //----------------------------------------------------------------------------

  // element rhs: standard Galerkin terms from potential equation
  CalcRhsPotEquDiviOhm(erhs,rhsfac,VarManager()->InvF(),VarManager()->GradPot(),1.);

  return;
}


/*--------------------------------------------------------------------------------*
 | CalcMat: linearization of diffusion coefficient in diffusion term   fang 02/15 |
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcMatDiffCoeffLin(
    Epetra_SerialDenseMatrix&                          emat,         //!< element matrix to be filled
    const int                                          k,            //!< index of current scalar
    const double                                       timefacfac,   //!< domain-integration factor times time-integration factor
    const LINALG::Matrix<my::nsd_,1>&                  gradphi,      //!< gradient of concentration at GP
    const double                                       scalar        //!< scaling factor for element matrix contributions
)
{
  // linearization of diffusion coefficient in ionic diffusion term (transport equation):
  //
  // (nabla w, D(D(c)) nabla c)
  //
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double laplawfrhs_gradphi(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi,gradphi,vi);

      emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k) += scalar*timefacfac*DiffManager()->GetDerivIsoDiffCoef(k,k)*laplawfrhs_gradphi*my::funct_(ui);
    }
  }

  return;
}


/*--------------------------------------------------------------------------------------------*
 | CalcMat: potential equation div i with inserted current - ohmic overpotential   fang 02/15 |
 *--------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcMatPotEquDiviOhm(
    Epetra_SerialDenseMatrix&                          emat,         //!< element matrix to be filled
    const double                                       timefacfac,   //!< domain-integration factor times time-integration factor
    const double                                       invf,         //!< 1/F
    const LINALG::Matrix<my::nsd_,1>&                  gradpot,      //!< gradient of potenial at GP
    const double                                       scalar        //!< scaling factor for element matrix contributions
)
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    for (int ui=0; ui<my::nen_; ++ui)
    {
      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf,ui,vi);

      // linearization of the ohmic term
      //
      // (grad w, 1/F kappa D(grad pot))
      //
      emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += scalar*timefacfac*invf*DiffManager()->GetCond()*laplawf;

      for(int iscal=0; iscal<my::numscal_; ++iscal)
      {
        double laplawfrhs_gradpot(0.);
        my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);

        // linearization of the ohmic term with respect to conductivity
        //
        // (grad w, 1/F kappa D(grad pot))
        //
        emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+iscal) += scalar*timefacfac*invf*DiffManager()->GetDerivCond(iscal)*my::funct_(ui)*laplawfrhs_gradpot;
      }
    }
  }

  return;
}


/*--------------------------------------------------------------------------------------------*
 | CalcRhs: potential equation div i with inserted current - ohmic overpotential   fang 02/15 |
 *--------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcRhsPotEquDiviOhm(
    Epetra_SerialDenseVector&                          erhs,      //!< element vector to be filled
    const double                                       rhsfac,    //!< time-integration factor for rhs times domain-integration factor
    const double                                       invf,      //!< 1./F
    const LINALG::Matrix<my::nsd_,1>&                  gradpot,   //!< gradient of potenial at GP
    const double                                       scalar     //!< scaling factor for element residual contributions
  )
{
  for (int vi=0; vi<my::nen_; ++vi)
  {
    double laplawfrhs_gradpot(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot,gradpot,vi);

    erhs[vi*my::numdofpernode_+my::numscal_] -= scalar*rhsfac*invf*DiffManager()->GetCond()*laplawfrhs_gradpot;
  }

  return;
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::GetMaterialParams(
    const DRT::Element*                  ele,           //!< the element we are dealing with
    std::vector<double>&                 densn,         //!< density at t_(n)
    std::vector<double>&                 densnp,        //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>&                 densam,        //!< density at t_(n+alpha_M)
    double&                              visc,          //!< fluid viscosity
    const int                            iquad          //!< id of current gauss point (default = -1)
    )
{
  // get material
  Teuchos::RCP<const MAT::Material> material = ele->Material();

  // evaluate electrode material
  if(material->MaterialType() == INPAR::MAT::m_electrode)
    Utils()->MatElectrode(material,VarManager()->Phinp(0),DiffManager());
  else
    dserror("Material type not supported!");

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::GetMaterialParams


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::nurbs27>;
