/*--------------------------------------------------------------------------*/
/*! \file

\brief supplementary element calculation class providing general utility for evaluation of heat
transport within electrochemical substances

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_sti_elch.hpp"

#include "baci_discretization_fem_general_extract_values.hpp"
#include "baci_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------------------------------------------------------------------*
 | element matrix and right-hand side vector contributions arising from thermal source terms in
 discrete thermo residuals   fang 11/15 |
 *-------------------------------------------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
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
template <CORE::FE::CellType distype>
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
template <CORE::FE::CellType distype>
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
  CORE::FE::ExtractMyValues(*elchnp, myelchnp, lm);
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
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleSTIElch<distype>::ScaTraEleSTIElch(
    const int numdofpernode, const int numscal, const std::string& disname)
    : econcnp_(true), epotnp_(true)
{
  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleSTIElch<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
