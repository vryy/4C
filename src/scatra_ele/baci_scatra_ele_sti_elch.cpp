/*--------------------------------------------------------------------------*/
/*! \file

\brief supplementary element calculation class providing general utility for evaluation of heat
transport within electrochemical substances

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_sti_elch.H"

#include "baci_lib_discret.H"
#include "baci_lib_utils.H"

/*-------------------------------------------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from thermal source terms in
 discrete thermo residuals   fang 11/15 |
 *-------------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleSTIElch<distype>::CalcMatAndRhsSource(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
    CORE::LINALG::SerialDenseVector& erhs,  //!< element right-hand side vector
    const double& timefacfac,  //!< domain integration factor times time integration factor
    const double& rhsfac       //!< domain integration factor times time integration factor for
                               //!< right-hand side vector
)
{
  // matrix and vector contributions arising from Joule's heat
  CalcMatAndRhsJoule(emat, erhs, timefacfac, rhsfac);

  // matrix and vector contributions arising from heat of mixing
  CalcMatAndRhsMixing(emat, erhs, timefacfac, rhsfac);

  // matrix and vector contributions arising from Soret effect
  CalcMatAndRhsSoret(emat, erhs, timefacfac, rhsfac);

  return;
};


/*-------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of source terms in discrete thermo residuals w.r.t.
 scatra dofs   fang 11/15 |
 *-------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleSTIElch<distype>::CalcMatSourceOD(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix
    const double& timefacfac  //!< domain integration factor times time integration factor
)
{
  // matrix contributions arising from Joule's heat
  CalcMatJouleOD(emat, timefacfac);

  // matrix contributions arising from heat of mixing
  CalcMatMixingOD(emat, timefacfac);

  // matrix contributions arising from Soret effect
  CalcMatSoretOD(emat, timefacfac);

  return;
};


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleSTIElch<distype>::ExtractElementAndNodeValues(
    DRT::Element* ele,                    //!< current element
    Teuchos::ParameterList& params,       //!< parameter list
    DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la       //!< location array
)
{
  // extract electrochemistry state vector from discretization
  const Teuchos::RCP<const Epetra_Vector> elchnp = discretization.GetState(2, "scatra");
  if (elchnp == Teuchos::null)
    dserror("Cannot extract electrochemistry state vector from discretization!");

  // extract local nodal values of concentration and electric potential from global state vector
  const std::vector<int>& lm = la[2].lm_;
  std::vector<double> myelchnp(lm.size());
  DRT::UTILS::ExtractMyValues(*elchnp, myelchnp, lm);
  for (int inode = 0; inode < nen_; ++inode)
  {
    econcnp_(inode) = myelchnp[inode * 2];
    epotnp_(inode) = myelchnp[inode * 2 + 1];
  }

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleSTIElch<distype>::ScaTraEleSTIElch(
    const int numdofpernode, const int numscal, const std::string& disname)
    : econcnp_(true), epotnp_(true)
{
  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::line2>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::tri3>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::tri6>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::quad8>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::quad9>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::hex20>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::hex27>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::tet4>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<DRT::Element::DiscretizationType::nurbs27>;
