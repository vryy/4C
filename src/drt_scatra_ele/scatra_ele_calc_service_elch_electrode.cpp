/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_electrode.cpp

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge within electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/material.H"

#include "scatra_ele.H"
#include "scatra_ele_calc_elch_electrode.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::EvaluateAction(
    DRT::ELEMENTS::Transport*   ele,
    Teuchos::ParameterList&     params,
    DRT::Discretization&        discretization,
    const SCATRA::Action&       action,
    const std::vector<int> &    lm,
    Epetra_SerialDenseMatrix&   elemat1_epetra,
    Epetra_SerialDenseMatrix&   elemat2_epetra,
    Epetra_SerialDenseVector&   elevec1_epetra,
    Epetra_SerialDenseVector&   elevec2_epetra,
    Epetra_SerialDenseVector&   elevec3_epetra
    )
{
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::calc_elch_electrode_soc:
  {
    CalculateElectrodeSOC(ele,params,discretization,lm);

    break;
  }
  default:
  {
    myelch::EvaluateAction(
        ele,
        params,
        discretization,
        action,
        lm,
        elemat1_epetra,
        elemat2_epetra,
        elevec1_epetra,
        elevec2_epetra,
        elevec3_epetra
        );

    break;
  }
  } // switch(action)

  return 0;
}


/*----------------------------------------------------------------------------------------------------------*
 | validity check with respect to input parameters, degrees of freedom, number of scalars etc.   fang 02/15 |
 *----------------------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CheckElchElementParameter(
    DRT::ELEMENTS::Transport*   ele   //!< current element
    )
{
  // safety checks
  if(ele->Material()->MaterialType() != INPAR::MAT::m_electrode)
    dserror("Invalid material type!");

  if(my::numscal_ != 1)
    dserror("Invalid number of transported scalars!");

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CheckElchElementParameter


/*----------------------------------------------------------------------*
 | get conductivity                                          fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::GetConductivity(
    const enum INPAR::ELCH::EquPot   equpot,      //!< type of closing equation for electric potential
    double&                          sigma_all,   //!< conductivity of electrolyte solution
    Epetra_SerialDenseVector&        sigma        //!< conductivity or a single ion + overall electrolyte solution
    )
{
  // use precomputed conductivity
  sigma_all = DiffManager()->GetCond();

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::GetConductivity


/*----------------------------------------------------------------------*
 | calculate electrode state of charge                       fang 01/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateElectrodeSOC(
    const DRT::Element*               ele,              //!< the element we are dealing with
    Teuchos::ParameterList&           params,           //!< parameter list
    DRT::Discretization&              discretization,   //!< discretization
    const std::vector<int>&           lm                //!< location vector
  )
{
  // consider only unghosted elements for integration
  if(ele->Owner() == discretization.Comm().MyPID())
  {
    // safety check
    if(my::numscal_ != 1)
      dserror("Electrode state of charge can only be computed for one transported scalar!");

    // get global state vector
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if(phinp == Teuchos::null)
      dserror("Cannot get state vector \"phinp\"!");

    // extract local nodal values from global state vector
    std::vector<double> ephinpvec(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinpvec,lm);
    for(int inode = 0; inode < my::nen_; ++inode)
      my::ephinp_[0](inode,0) = ephinpvec[inode*my::numdofpernode_];

    // get current values of concentration and domain integrals
    double intconcentration = params.get<double>("intconcentration");
    double intdomain = params.get<double>("intdomain");

    // integration points and weights
    const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    // loop over integration points
    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      // evaluate values of shape functions and domain integration factor at current integration point
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      // calculate concentration and domain integrals
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double funct_vi_fac = my::funct_(vi)*fac;

        // concentration integral
        intconcentration += funct_vi_fac*my::ephinp_[0](vi,0);

        // domain integral
        intdomain += funct_vi_fac;
      }
    } // loop over integration points

    // write updated values of concentration and domain integrals back to parameter list
    params.set<double>("intconcentration",intconcentration);
    params.set<double>("intdomain",intdomain);
  } // if(ele->Owner() == discretization.Comm().MyPID())

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateElectrodeSOC


/*---------------------------------------------------------------------------------------*
 | compute element matrix and residual for initial electric potential field   fang 02/15 |
 *---------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcMatAndRhsElectricPotentialField(
    const enum INPAR::ELCH::EquPot   equpot,   //!< type of closing equation for electric potential
    Epetra_SerialDenseMatrix&        emat,     //!< element matrix
    Epetra_SerialDenseVector&        erhs,     //!< element rhs
    const double                     fac,      //!< integration factor
    const double                     scalar    //!< scaling factor for element matrix and residual contributions
    )
{
  // provide standard mass matrix for degrees of freedom associated with concentration
  for (int k=0; k<my::numscal_; ++k)
  {
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const int fvi = vi*my::numdofpernode_+k;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;

        emat(fvi,fui) += fac*my::funct_(vi)*my::funct_(ui);
      }
    }
  }

  // provide standard Galerkin terms for degrees of freedom associated with electric potential
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+my::numscal_;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+my::numscal_;

      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf,ui,vi);

      emat(fvi,fui) += scalar*fac*VarManager()->InvF()*DiffManager()->GetCond()*laplawf;
    }

    double laplawf(0.);
    my::GetLaplacianWeakFormRHS(laplawf,VarManager()->GradPot(),vi);

    erhs[fvi] -= scalar*fac*VarManager()->InvF()*DiffManager()->GetCond()*laplawf;
  }

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalcMatAndRhsElectricPotentialField


/*---------------------------------------------------------------------*
 | calculate weighted mass flux (no reactive flux so far)   fang 02/15 |
 *---------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateFlux(
    LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac         //!< integration factor
    )
{
  /*
    Actually, we compute here a weighted (and integrated) form of the fluxes!
    On time integration level, these contributions are then used to calculate
    an L2-projected representation of fluxes.
    Thus, this method here DOES NOT YET provide flux values that are ready to use!
  */

  // add convective flux contribution
  switch (fluxtype)
  {
  case INPAR::SCATRA::flux_total_domain:
  {
    q.Update(VarManager()->ConInt(k),VarManager()->ConVelInt());
    break;
  }

  default:
  {
    dserror("received illegal flag inside flux evaluation for whole domain");
    break;
  }
  };

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateFlux


/*------------------------------------------------------------------------------*
 | set internal variables for electrodes                             fang 02/15 |
 *------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchElectrode(my::funct_,my::derxy_,my::ephinp_,my::ephin_,myelch::epotnp_,my::econvelnp_,my::ehist_,DiffManager());

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::SetInternalVariablesForMatAndRHS()


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri6>;
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
