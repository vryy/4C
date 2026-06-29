// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_scatra_ele_boundary_calc_elch_NP.hpp"

#include "4C_mat_ion.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>*
Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>::instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::Utils::make_singleton_map<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchNP<distype, probdim>>(
            new ScaTraEleBoundaryCalcElchNP<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].instance(
      Core::Utils::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>::ScaTraEleBoundaryCalcElchNP(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      myelch::ScaTraEleBoundaryCalcElch(numdofpernode, numscal, disname)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>::evaluate_action(
    Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, ScaTra::BoundaryAction action,
    Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1,
    Core::LinAlg::SerialDenseMatrix& elemat2, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseVector& elevec2, Core::LinAlg::SerialDenseVector& elevec3)
{
  // determine and evaluate action
  switch (action)
  {
    case ScaTra::BoundaryAction::calc_elch_boundary_kinetics:
    {
      myelch::calc_elch_boundary_kinetics(ele, params, discretization, la, elemat1, elevec1, 1.);

      break;
    }

    default:
    {
      myelch::evaluate_action(
          ele, params, discretization, action, la, elemat1, elemat2, elevec1, elevec2, elevec3);

      break;
    }
  }  // switch action

  return 0;
}  // Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>::evaluate_action

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>::evaluate_neumann(
    Core::Elements::FaceElement* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, const Core::Conditions::Condition& condition,
    Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseVector& elevec1,
    const double scalar)
{
  // call base class routine
  my::evaluate_neumann(ele, params, discretization, condition, la, elevec1, scalar);

  // add boundary flux contributions to potential equation
  switch (myelch::elchparams_->equ_pot())
  {
    case ElCh::equpot_enc_pde:
    case ElCh::equpot_enc_pde_elim:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        // get valence
        const double valence_k = get_valence(ele->parent_element()->material(), k);

        for (int vi = 0; vi < nen_; ++vi)
          elevec1[vi * my::numdofpernode_ + my::numscal_] +=
              valence_k * elevec1[vi * my::numdofpernode_ + k];
      }

      break;
    }

    case ElCh::equpot_enc:
    case ElCh::equpot_poisson:
    case ElCh::equpot_laplace:
      // do nothing in these cases
      break;

    default:
    {
      FOUR_C_THROW("Closing equation for electric potential not recognized!");
    }
  }  // switch(myelch::elchparams_->EquPot())

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype,
    probdim>::evaluate_elch_boundary_kinetics(const Core::Elements::Element* ele,
    Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs,
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ephinp,
    const std::vector<Core::LinAlg::Matrix<nen_, 1>>& ehist, double timefac,
    std::shared_ptr<const Core::Mat::Material> material, const Core::Conditions::Condition& cond,
    const int nume, const std::vector<int>& stoich, const int kinetics, const double pot0,
    const double frt, const double scalar)
{
  // call base class routine
  myelch::evaluate_elch_boundary_kinetics(ele, emat, erhs, ephinp, ehist, timefac, material, cond,
      nume, stoich, kinetics, pot0, frt, scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch (myelch::elchparams_->equ_pot())
  {
    case ElCh::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case ElCh::equpot_enc_pde:
    case ElCh::equpot_enc_pde_elim:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (int vi = 0; vi < nen_; ++vi)
        {
          for (int ui = 0; ui < nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] += nume * erhs[vi * my::numdofpernode_ + k];
        }
      }

      break;
    }

    // need special treatment for Laplace equation due to missing scaling with inverse of Faraday
    // constant
    case ElCh::equpot_laplace:
    {
      const double faraday = myelch::elchparams_->faraday();
      for (int k = 0; k < my::numscal_; ++k)
      {
        for (int vi = 0; vi < nen_; ++vi)
        {
          for (int ui = 0; ui < nen_; ++ui)
          {
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                faraday * nume * emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k);
            emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
                faraday * nume *
                emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + my::numscal_);
          }

          erhs[vi * my::numdofpernode_ + my::numscal_] +=
              faraday * nume * erhs[vi * my::numdofpernode_ + k];
        }
      }

      break;
    }

    case ElCh::equpot_poisson:
    {
      FOUR_C_THROW("Poisson equation combined with electrode boundary conditions not implemented!");
    }

    default:
    {
      FOUR_C_THROW("Unknown closing equation for electric potential!");
    }
  }  // switch(myelch::elchparams_->EquPot())

}  // Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype,
   // probdim>::evaluate_elch_boundary_kinetics


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
double Discret::Elements::ScaTraEleBoundaryCalcElchNP<distype, probdim>::get_valence(
    const std::shared_ptr<const Core::Mat::Material>& material, const int k) const
{
  double valence(0.);

  if (material->material_type() == Core::Materials::m_matlist)
  {
    const std::shared_ptr<const Mat::MatList> matlist =
        std::static_pointer_cast<const Mat::MatList>(material);

    const std::shared_ptr<const Core::Mat::Material> species =
        matlist->material_by_id(matlist->mat_id(k));

    if (species->material_type() == Core::Materials::m_ion)
    {
      valence = std::static_pointer_cast<const Mat::Ion>(species)->valence();
      FOUR_C_ASSERT_ALWAYS(std::abs(valence) > 1.0e-16, "Received zero valence!");
    }
    else
    {
      FOUR_C_THROW("Material species is not an ion!");
    }
  }

  else
  {
    FOUR_C_THROW("Unknown material!");
  }

  return valence;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::quad4, 3>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::quad8, 3>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::quad9, 3>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::tri3, 3>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::tri6, 3>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::line2, 2>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::line2, 3>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::line3, 2>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::nurbs3, 2>;
template class Discret::Elements::ScaTraEleBoundaryCalcElchNP<Core::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
