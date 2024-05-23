/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra boundary elements for Nernst-Planck formulation

\level 2

 */
/*----------------------------------------------------------------------*/
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
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>*
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleBoundaryCalcElchNP<distype, probdim>>(
            new ScaTraEleBoundaryCalcElchNP<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::ScaTraEleBoundaryCalcElchNP(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      myelch::ScaTraEleBoundaryCalcElch(numdofpernode, numscal, disname)
{
  return;
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 08/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele,                            //!< boundary element
    Teuchos::ParameterList& params,                   //!< parameter list
    DRT::Discretization& discretization,              //!< discretization
    SCATRA::BoundaryAction action,                    //!< action
    DRT::Element::LocationArray& la,                  //!< location array
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_elch_boundary_kinetics:
    {
      myelch::calc_elch_boundary_kinetics(
          ele, params, discretization, la, elemat1_epetra, elevec1_epetra, 1.);

      break;
    }

    default:
    {
      myelch::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch action

  return 0;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::EvaluateAction


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                       fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::EvaluateNeumann(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    CORE::Conditions::Condition& condition, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& elevec1, const double scalar)
{
  // call base class routine
  my::EvaluateNeumann(ele, params, discretization, condition, la, elevec1, scalar);

  // add boundary flux contributions to potential equation
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    {
      for (int k = 0; k < my::numscal_; ++k)
      {
        // get valence
        const double valence_k = GetValence(ele->ParentElement()->Material(), k);

        for (int vi = 0; vi < nen_; ++vi)
          elevec1[vi * my::numdofpernode_ + my::numscal_] +=
              valence_k * elevec1[vi * my::numdofpernode_ + k];
      }

      break;
    }

    case INPAR::ELCH::equpot_enc:
    case INPAR::ELCH::equpot_poisson:
    case INPAR::ELCH::equpot_laplace:
      // do nothing in these cases
      break;

    default:
    {
      FOUR_C_THROW("Closing equation for electric potential not recognized!");
      break;
    }
  }  // switch(myelch::elchparams_->EquPot())

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition         fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::evaluate_elch_boundary_kinetics(
    const DRT::Element* ele,                ///< current element
    CORE::LINALG::SerialDenseMatrix& emat,  ///< element matrix
    CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
    double timefac,                                           ///< time factor
    Teuchos::RCP<const CORE::MAT::Material> material,         ///< material
    Teuchos::RCP<CORE::Conditions::Condition> cond,  ///< electrode kinetics boundary condition
    const int nume,                                  ///< number of transferred electrons
    const std::vector<int> stoich,                   ///< stoichiometry of the reaction
    const int kinetics,                              ///< desired electrode kinetics model
    const double pot0,                               ///< electrode potential on metal side
    const double frt,                                ///< factor F/RT
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // call base class routine
  myelch::evaluate_elch_boundary_kinetics(ele, emat, erhs, ephinp, ehist, timefac, material, cond,
      nume, stoich, kinetics, pot0, frt, scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch (myelch::elchparams_->EquPot())
  {
    case INPAR::ELCH::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
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
    case INPAR::ELCH::equpot_laplace:
    {
      const double faraday = myelch::elchparams_->Faraday();
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

    case INPAR::ELCH::equpot_poisson:
    {
      FOUR_C_THROW("Poisson equation combined with electrode boundary conditions not implemented!");
      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown closing equation for electric potential!");
      break;
    }
  }  // switch(myelch::elchparams_->EquPot())

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::evaluate_elch_boundary_kinetics


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<distype, probdim>::GetValence(
    const Teuchos::RCP<const CORE::MAT::Material>& material,  // element material
    const int k                                               // species number
) const
{
  double valence(0.);

  if (material->MaterialType() == CORE::Materials::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList> matlist =
        Teuchos::rcp_static_cast<const MAT::MatList>(material);

    const Teuchos::RCP<const CORE::MAT::Material> species =
        matlist->MaterialById(matlist->MatID(k));

    if (species->MaterialType() == CORE::Materials::m_ion)
    {
      valence = Teuchos::rcp_static_cast<const MAT::Ion>(species)->Valence();
      if (abs(valence) < 1.e-14) FOUR_C_THROW("Received zero valence!");
    }
    else
      FOUR_C_THROW("Material species is not an ion!");
  }

  else
    FOUR_C_THROW("Unknown material!");

  return valence;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchNP<CORE::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
