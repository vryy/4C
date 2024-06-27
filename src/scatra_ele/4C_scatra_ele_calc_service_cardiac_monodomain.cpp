/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for cardiac monodomain problems

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_myocard.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_calc_cardiac_monodomain.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::evaluate_action(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, const ScaTra::Action& action,
    Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  //(for now) only first dof set considered
  const std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  switch (action)
  {
    case ScaTra::Action::time_update_material:
    {
      std::vector<Teuchos::RCP<Mat::Myocard>> updatemat;
      updatemat.reserve(my::numscal_);

      // access the general material
      Teuchos::RCP<Core::Mat::Material> material = ele->Material();

      // first, determine the materials which need a time update, i.e. myocard materials
      if (material->MaterialType() == Core::Materials::m_matlist)
      {
        const Teuchos::RCP<Mat::MatList> actmat = Teuchos::rcp_dynamic_cast<Mat::MatList>(material);
        if (actmat->NumMat() < my::numscal_) FOUR_C_THROW("Not enough materials in MatList.");

        for (int k = 0; k < my::numscal_; ++k)
        {
          const int matid = actmat->MatID(k);
          Teuchos::RCP<Core::Mat::Material> singlemat = actmat->MaterialById(matid);

          if (singlemat->MaterialType() == Core::Materials::m_myocard)
          {
            // reference to Teuchos::rcp not possible here, since the material
            // is required to be not const for this application
            updatemat.push_back(Teuchos::rcp_dynamic_cast<Mat::Myocard>(singlemat));
          }
        }
      }

      if (material->MaterialType() == Core::Materials::m_myocard)
      {  // reference to Teuchos::rcp not possible here, since the material is required to be
        // not const for this application
        updatemat.push_back(Teuchos::rcp_dynamic_cast<Mat::Myocard>(material));
      }

      if (updatemat.size() > 0)  // found at least one material to be updated
      {
        // all materials in the matlist should be of the same kind
        if (updatemat.size() != (unsigned)my::numscal_)
          FOUR_C_THROW("Number of materials to be updated is not equal to number of scalars!");

        // get time-step length
        const double dt = my::scatraparatimint_->Dt();

        // extract local values from the global vectors
        Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
        if (phinp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'phinp'");
        Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_, 1>>(*phinp, my::ephinp_, lm);

        my::eval_shape_func_and_derivs_at_ele_center();

        for (unsigned i = 0; i < updatemat.size(); i++)
        {
          const double csnp = my::funct_.dot(my::ephinp_[i]);  // be careful, we assume k==i here
          updatemat[i]->update(csnp, dt);
        }
      }

      break;
    }

    case ScaTra::Action::get_material_internal_state:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // access the general material
        Teuchos::RCP<Core::Mat::Material> material = ele->Material();
        Teuchos::RCP<Epetra_MultiVector> material_internal_state =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");

        if (material->MaterialType() == Core::Materials::m_myocard)
        {
          Teuchos::RCP<Mat::Myocard> material =
              Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
          for (int k = 0; k < material_internal_state->NumVectors(); ++k)
          {
            int err = material_internal_state->ReplaceGlobalValue(
                ele->Id(), k, material->GetInternalState(k));
            if (err != 0) FOUR_C_THROW("%i", err);
          }
        }
        params.set<Teuchos::RCP<Epetra_MultiVector>>(
            "material_internal_state", material_internal_state);
      }

      break;
    }

    case ScaTra::Action::set_material_internal_state:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // access the general material
        Teuchos::RCP<Core::Mat::Material> material = ele->Material();
        Teuchos::RCP<Epetra_Vector> material_internal_state_component =
            params.get<Teuchos::RCP<Epetra_Vector>>("material_internal_state_component");

        if (material->MaterialType() == Core::Materials::m_myocard)
        {
          Teuchos::RCP<Mat::Myocard> material =
              Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
          int k = params.get<int>("k");
          material->SetInternalState(k, (*material_internal_state_component)[ele->Id()]);
        }
      }
    }

    break;

    case ScaTra::Action::get_material_ionic_currents:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // access the general material
        Teuchos::RCP<Core::Mat::Material> material = ele->Material();
        Teuchos::RCP<Epetra_MultiVector> material_ionic_currents =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("material_ionic_currents");

        if (material->MaterialType() == Core::Materials::m_myocard)
        {
          Teuchos::RCP<Mat::Myocard> material =
              Teuchos::rcp_dynamic_cast<Mat::Myocard>(ele->Material());
          for (int k = 0; k < material_ionic_currents->NumVectors(); ++k)
          {
            int err = material_ionic_currents->ReplaceGlobalValue(
                ele->Id(), k, material->GetIonicCurrents(k));
            if (err != 0) FOUR_C_THROW("%i", err);
          }
        }
        params.set<Teuchos::RCP<Epetra_MultiVector>>(
            "material_ionic_currents", material_ionic_currents);
      }

      break;
    }

    default:
    {
      my::evaluate_action(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch(action)

  return 0;
}


// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad4, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::quad9, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::hex8, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::tet10, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::pyramid5, 3>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcCardiacMonodomain<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
