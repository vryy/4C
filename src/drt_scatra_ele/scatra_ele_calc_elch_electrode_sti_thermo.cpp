/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch_electrode_sti_thermo.cpp

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within thermodynamic electrodes

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_electrode_sti_thermo.H"

#include "scatra_ele_parameter_timint.H"

#include "../drt_mat/material.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>*
DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname,
    const ScaTraEleCalcElchElectrodeSTIThermo* delete_me)
{
  static std::map<std::string, ScaTraEleCalcElchElectrodeSTIThermo<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleCalcElchElectrodeSTIThermo<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcElchElectrodeSTIThermo<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::ExtractElementAndNodeValues(
    DRT::Element* ele,                    //!< current element
    Teuchos::ParameterList& params,       //!< parameter list
    DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la       //!< location array
)
{
  // call base class routine to extract scatra-related quantities
  myelch::ExtractElementAndNodeValues(ele, params, discretization, la);

  // call base class routine to extract thermo-related quantitites
  mythermo::ExtractElementAndNodeValues(ele, params, discretization, la);

  return;
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::GetMaterialParams(
    const DRT::Element* ele,      //!< current element
    std::vector<double>& densn,   //!< density at t_(n)
    std::vector<double>& densnp,  //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,  //!< density at t_(n+alpha_M)
    double& visc,                 //!< fluid viscosity
    const int iquad               //!< ID of current integration point
)
{
  // call base class routine to get parameters of primary, isothermal electrode material
  myelectrode::GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

  // get parameters of secondary, thermodynamic electrolyte material
  Teuchos::RCP<const MAT::Material> material = ele->Material(1);
  if (material->MaterialType() == INPAR::MAT::m_soret)
    mythermo::MatSoret(material);
  else
    dserror("Invalid electrode material!");

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::GetMaterialParams


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::CalcMatAndRhs(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
    const int k,                     //!< index of current scalar
    const double fac,                //!< domain integration factor
    const double timefacfac,         //!< domain integration factor times time integration factor
    const double rhsfac,      //!< domain integration factor times time integration factor for
                              //!< right-hand side vector
    const double taufac,      //!< domain integration factor times stabilization parameter
    const double timetaufac,  //!< domain integration factor times stabilization parameter times
                              //!< time integration factor
    const double rhstaufac,  //!< domain integration factor times stabilization parameter times time
                             //!< integration factor for right-hand side vector
    LINALG::Matrix<my::nen_, 1>&
        tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
    double& rhsint  //!< body force value
)
{
  // call base class routine for isothermal problems
  myelectrode::CalcMatAndRhs(
      emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac, tauderpot, rhsint);

  // matrix and vector contributions arising from additional, thermodynamic term for Soret effect
  mythermo::CalcMatSoret(emat, timefacfac, VarManager()->Phinp(0),
      myelectrode::DiffManager()->GetIsotropicDiff(0),
      myelectrode::DiffManager()->GetDerivIsoDiffCoef(0, 0), VarManager()->Temp(),
      VarManager()->GradTemp(), my::funct_, my::derxy_);
  mythermo::CalcRHSSoret(erhs, VarManager()->Phinp(0),
      myelectrode::DiffManager()->GetIsotropicDiff(0), rhsfac, VarManager()->Temp(),
      VarManager()->GradTemp(), my::derxy_);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::EvaluateActionOD(
    DRT::Element* ele,                         //!< current element
    Teuchos::ParameterList& params,            //!< parameter list
    DRT::Discretization& discretization,       //!< discretization
    const SCATRA::Action& action,              //!< action parameter
    DRT::Element::LocationArray& la,           //!< location array
    Epetra_SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    Epetra_SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    Epetra_SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    Epetra_SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    Epetra_SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::calc_scatra_mono_odblock_scatrathermo:
    {
      SysmatODScatraThermo(ele, elemat1_epetra);

      break;
    }

    default:
    {
      // call base class routine
      my::EvaluateActionOD(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return 0;
}


/*------------------------------------------------------------------------------------------------------*
 | fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs   fang
 11/15 |
 *------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::SysmatODScatraThermo(
    DRT::Element* ele,              //!< current element
    Epetra_SerialDenseMatrix& emat  //!< element matrix
)
{
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions, their derivatives, and domain integration factor at current
    // integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

    // evaluate internal variables at current integration point
    SetInternalVariablesForMatAndRHS();

    // evaluate material parameters at current integration point
    double dummy(0.);
    std::vector<double> dummyvec(my::numscal_, 0.);
    GetMaterialParams(ele, dummyvec, dummyvec, dummyvec, dummy, iquad);

    // provide element matrix with linearizations of Soret term in discrete scatra residuals w.r.t.
    // thermo dofs
    mythermo::CalcMatSoretOD(emat, timefacfac, VarManager()->Phinp(0),
        myelectrode::DiffManager()->GetIsotropicDiff(0), VarManager()->Temp(),
        VarManager()->GradTemp(), my::funct_, my::derxy_);
  }

  return;
}


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables for element evaluation
  VarManager()->SetInternalVariables(my::funct_, my::derxy_, my::ephinp_, my::ephin_,
      mythermo::etempnp_, my::econvelnp_, my::ehist_);

  return;
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::ScaTraEleCalcElchElectrodeSTIThermo(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructors of base classes
      ScaTraEleCalcElchElectrode<distype>::ScaTraEleCalcElchElectrode(
          numdofpernode, numscal, disname),
      ScaTraEleSTIThermo<distype>::ScaTraEleSTIThermo(numscal)
{
  // safety check
  if (numscal != 1 or numdofpernode != 2)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // replace internal variable manager for isothermal electrodes by internal variable manager for
  // thermodynamic electrodes
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<my::nsd_, my::nen_>(
          my::numscal_, myelch::elchparams_));

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<DRT::Element::nurbs27>;
