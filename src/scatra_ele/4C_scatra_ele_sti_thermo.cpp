/*--------------------------------------------------------------------------*/
/*! \file

\brief supplementary element calculation class providing general utility for thermodynamic scalar
transport

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_sti_thermo.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_mat_soret.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::extract_element_and_node_values(
    Core::Elements::Element* ele,               //!< current element
    Teuchos::ParameterList& params,             //!< parameter list
    Discret::Discretization& discretization,    //!< discretization
    Core::Elements::Element::LocationArray& la  //!< location array
)
{
  // extract thermo state vector from discretization
  Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(2, "thermo");
  if (tempnp == Teuchos::null)
    FOUR_C_THROW("Cannot extract thermo state vector from discretization!");

  // extract local nodal temperature values from global state vector
  Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_, 1>>(*tempnp, etempnp_, la[2].lm_);
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of Soret effect term in discrete scatra residuals
 w.r.t. scatra dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::calc_mat_soret(
    Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
    const double& timefacfac,      //!< domain integration factor times time integration factor
    const double& conc,            //!< concentration
    const double& diffcoeff,       //!< diffusion coefficient
    const double& diffcoeffderiv,  //!< derivative of diffusion coefficient w.r.t. concentration
    const double& temp,            //!< temperature
    const Core::LinAlg::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
    const Core::LinAlg::Matrix<nen_, 1>& funct,     //!< shape functions
    const Core::LinAlg::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
)
{
  for (int vi = 0; vi < nen_; ++vi)
  {
    // recurring index
    const int rowconc(vi * 2);

    // gradient of test function times temperature gradient
    double laplawfrhs_gradtemp(0.);
    get_laplacian_weak_form_rhs(laplawfrhs_gradtemp, vi, gradtemp, derxy);

    for (int ui = 0; ui < nen_; ++ui)
    {
      // linearizations of Soret effect term in concentration residuals w.r.t. concentration dofs
      emat(rowconc, ui * 2) += timefacfac * funct(ui) * laplawfrhs_gradtemp / temp *
                               diffmanagerstithermo_->GetSoret() *
                               (diffcoeff + conc * diffcoeffderiv);
    }

    // linearizations of Soret effect term in concentration residuals w.r.t. electric potential dofs
    // are zero Soret effect term does not appear in electric potential residuals
  }
}


/*------------------------------------------------------------------------------------------------------------------------------*
 | provide element matrix with linearizations of Soret effect term in discrete scatra residuals
 w.r.t. thermo dofs   fang 11/15 |
 *------------------------------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::calc_mat_soret_od(
    Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
    const double& timefacfac,     //!< time integration factor times domain integration factor
    const double& concentration,  //!< concentration
    const double& diffcoeff,      //!< diffusion coefficient
    const double& temperature,    //!< temperature
    const Core::LinAlg::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
    const Core::LinAlg::Matrix<nen_, 1>& funct,     //!< shape functions
    const Core::LinAlg::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
)
{
  for (int vi = 0; vi < nen_; ++vi)
  {
    // recurring index
    const int rowconc = vi * 2;

    // gradient of test function times temperature gradient
    double laplawfrhs_gradtemp(0.);
    get_laplacian_weak_form_rhs(laplawfrhs_gradtemp, vi, gradtemp, derxy);

    for (int ui = 0; ui < nen_; ++ui)
    {
      // gradient of test function times gradient of shape function
      double laplawf(0.);
      get_laplacian_weak_form(laplawf, vi, ui, derxy);

      // linearizations of Soret effect term in concentration residuals w.r.t. thermo dofs
      emat(rowconc, ui) +=
          timefacfac * diffcoeff * diffmanagerstithermo_->GetSoret() * concentration *
          (laplawf / temperature - funct(ui) * laplawfrhs_gradtemp / (pow(temperature, 2)));

      // Soret effect term does not appear in electric potential residuals
    }
  }
}


/*--------------------------------------------------------------------------------------------------------------------------*
 | provide element right-hand side vector with contributions of Soret effect term to discrete scatra
 residuals   fang 11/15 |
 *--------------------------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::calc_rhs_soret(
    Core::LinAlg::SerialDenseVector& erhs,  //!< element right-hand side vector
    const double& rhsfac,     //!< domain integration factor times time integration factor for
                              //!< right-hand side vector
    const double& conc,       //!< concentration
    const double& diffcoeff,  //!< diffusion coefficient
    const double& temp,       //!< temperature
    const Core::LinAlg::Matrix<nsd_, 1>& gradtemp,  //!< gradient of temperature
    const Core::LinAlg::Matrix<nsd_, nen_>& derxy   //!< spatial derivatives of shape functions
)
{
  for (int vi = 0; vi < nen_; ++vi)
  {
    // gradient of test function times temperature gradient
    double laplawfrhs_gradtemp(0.);
    get_laplacian_weak_form_rhs(laplawfrhs_gradtemp, vi, gradtemp, derxy);

    // contributions of Soret effect term to concentration residuals
    erhs[vi * 2] -=
        rhsfac * laplawfrhs_gradtemp * diffcoeff * conc * diffmanagerstithermo_->GetSoret() / temp;

    // Soret effect term does not appear in electric potential residuals
  }
}


/*----------------------------------------------------------------------*
 | evaluate Soret material                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::mat_soret(
    const Teuchos::RCP<const Core::Mat::Material> material  //!< Soret material
)
{
  // extract material parameters from Soret material
  const Teuchos::RCP<const Mat::Soret> matsoret =
      Teuchos::rcp_static_cast<const Mat::Soret>(material);
  diffmanagerstithermo_->SetIsotropicDiff(matsoret->Conductivity(), 0);
  diffmanagerstithermo_->SetSoret(matsoret->SoretCoefficient());
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::ScaTraEleSTIThermo(
    const int& numscal  //!< number of transported scalars
    )
    : etempnp_(true),

      // initialize thermo diffusion manager
      diffmanagerstithermo_(Teuchos::rcp(new ScaTraEleDiffManagerSTIThermo(numscal)))
{
  // safety check
  if (numscal != 1)
  {
    FOUR_C_THROW(
        "Thermodynamic scalar transport only works for exactly one transported scalar at the "
        "moment!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleSTIThermo<distype>::calc_mat_diff_thermo_od(
    Core::LinAlg::SerialDenseMatrix& emat, const int& numdofpernode, const double& timefacfac,
    const double& invF, const Core::LinAlg::Matrix<nsd_, 1>& gradconc,
    const Core::LinAlg::Matrix<nsd_, 1>& gradpot, const double& tempderivisodiffcoef,
    const double& tempderivcond, const Core::LinAlg::Matrix<nen_, 1>& funct,
    const Core::LinAlg::Matrix<nsd_, nen_>& derxy, const double& scalefac)
{
  for (int vi = 0; vi < static_cast<int>(nen_); ++vi)
  {
    const int rowconc = vi * 2;
    const int rowpot = rowconc + 1;

    for (int ui = 0; ui < static_cast<int>(nen_); ++ui)
    {
      double laplawfrhs(0.0);
      get_laplacian_weak_form_rhs(laplawfrhs, vi, gradconc, derxy);
      emat(rowconc, ui) += timefacfac * tempderivisodiffcoef * laplawfrhs * funct(ui);

      get_laplacian_weak_form_rhs(laplawfrhs, vi, gradpot, derxy);
      emat(rowpot, ui) += timefacfac * invF * tempderivcond * laplawfrhs * funct(ui) * scalefac;
    }
  }
}


// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::quad4>;
// template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::hex8>;
// template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::tet10>;
// template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::pyramid5>;
// template class Discret::ELEMENTS::ScaTraEleSTIThermo<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
