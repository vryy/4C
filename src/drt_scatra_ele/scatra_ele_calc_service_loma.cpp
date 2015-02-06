/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_loma.cpp

\brief evaluation of scatra elements for loma

<pre>
Maintainer: Ursula Rasthofer / Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15236/45
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_loma.H"

#include "scatra_ele.H"
#include "scatra_ele_action.H"

#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::EvaluateAction(
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
  if (action == SCATRA::calc_dissipation or action == SCATRA::calc_mean_Cai)
  {
    if (my::scatraparatimint_->IsGenAlpha())
    {
      // extract additional local values from global vector
      Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
      if (phiam==Teuchos::null) dserror("Cannot get state vector 'phiam'");
      std::vector<double> myphiam(lm.size());
      DRT::UTILS::ExtractMyValues(*phiam,myphiam,lm);

      // fill element array
      for (int i=0;i<my::nen_;++i)
      {
        for (int k = 0; k< my::numscal_; ++k)
        {
          // split for each transported scalar, insert into element arrays
          ephiam_[k](i,0) = myphiam[k+(i*my::numdofpernode_)];
        }
       } // for i
    }
  }

  if (action == SCATRA::calc_dissipation or action == SCATRA::calc_mean_Cai)
  {
    // get thermodynamic pressure (for CalcDissipation)
    thermpressnp_ = params.get<double>("thermodynamic pressure");
    thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
    if (my::scatraparatimint_->IsGenAlpha())
      thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");
  }

  switch (action)
  {
  case SCATRA::calc_domain_and_bodyforce:
  {
    // NOTE: add integral values only for elements which are NOT ghosted!
    if (ele->Owner() == discretization.Comm().MyPID())
      // calculate domain and bodyforce integral
      CalculateDomainAndBodyforce(elevec1_epetra,ele);

    break;
  }
  default:
  {
    my::EvaluateAction(ele,
                       params,
                       discretization,
                       action,
                       lm,
                       elemat1_epetra,
                       elemat2_epetra,
                       elevec1_epetra,
                       elevec2_epetra,
                       elevec3_epetra);
    break;
  }
  } // switch(action)

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate domain integral                                   vg 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::CalculateDomainAndBodyforce(
Epetra_SerialDenseVector&  scalars,
const DRT::Element*        ele
)
{
  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::BodyForce(ele);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // get bodyforce in gausspoint
    const double rhs = my::bodyforce_[0].Dot(my::funct_);

    // calculate integrals of domain and bodyforce
    for (int i=0; i<my::nen_; i++)
    {
      scalars[0] += fac*my::funct_(i);
    }
    scalars[1] += fac*rhs;

  } // loop over integration points

  return;
} // ScaTraEleCalcLoma::CalculateDomain


/*-----------------------------------------------------------------------------------------*
 | extract element based or nodal values and return extracted values of phinp   fang 02/15 |
 *-----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
std::vector<double>  DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::ExtractElementAndNodeValues(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm
)
{
  // add loma-specific values
  if(my::scatraparatimint_->IsGenAlpha())
  {
    // extract local values from global vector
    Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
    if (phiam==Teuchos::null) dserror("Cannot get state vector 'phiam'");
    std::vector<double> myphiam(lm.size());
    DRT::UTILS::ExtractMyValues(*phiam,myphiam,lm);

    // fill element array
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0; k< my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        ephiam_[k](i,0) = myphiam[k+(i*my::numdofpernode_)];
      }
    } // for i
  }

  // get thermodynamic pressure
  thermpressnp_ = params.get<double>("thermodynamic pressure");
  thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
  if (my::scatraparatimint_->IsGenAlpha())
    thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");

  // call base class routine
  return my::ExtractElementAndNodeValues(ele,params,discretization,lm);
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcLoma<DRT::Element::nurbs27>;
