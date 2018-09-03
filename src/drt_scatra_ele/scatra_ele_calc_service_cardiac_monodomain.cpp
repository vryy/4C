/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_cardiac_monodomain.cpp

\brief evaluation of scatra elements for cardiac monodomain problems

\level 2

<pre>
\maintainer Lasse Jagschies
            lasse.jagschies@tum.de
            http://www.lnm.mw.tum.de
            089-289-10365
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_cardiac_monodomain.H"

#include "scatra_ele.H"
#include "scatra_ele_action.H"

#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_mat/myocard.H"
#include "../drt_mat/matlist.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::EvaluateAction(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  //(for now) only first dof set considered
  const std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  switch (action)
  {
    case SCATRA::time_update_material:
    {
      std::vector<Teuchos::RCP<MAT::Myocard>> updatemat;
      updatemat.reserve(my::numscal_);

      // access the general material
      Teuchos::RCP<MAT::Material> material = ele->Material();

      // first, determine the materials which need a time update, i.e. myocard materials
      if (material->MaterialType() == INPAR::MAT::m_matlist)
      {
        const Teuchos::RCP<MAT::MatList> actmat = Teuchos::rcp_dynamic_cast<MAT::MatList>(material);
        if (actmat->NumMat() < my::numscal_) dserror("Not enough materials in MatList.");

        for (int k = 0; k < my::numscal_; ++k)
        {
          const int matid = actmat->MatID(k);
          Teuchos::RCP<MAT::Material> singlemat = actmat->MaterialById(matid);

          if (singlemat->MaterialType() == INPAR::MAT::m_myocard)
          {
            // reference to Teuchos::rcp not possible here, since the material
            // is required to be not const for this application
            updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(singlemat));
          }
        }
      }

      if (material->MaterialType() == INPAR::MAT::m_myocard)
      {  // reference to Teuchos::rcp not possible here, since the material is required to be
        // not const for this application
        updatemat.push_back(Teuchos::rcp_dynamic_cast<MAT::Myocard>(material));
      }

      if (updatemat.size() > 0)  // found at least one material to be updated
      {
        // all materials in the matlist should be of the same kind
        if (updatemat.size() != (unsigned)my::numscal_)
          dserror("Number of materials to be updated is not equal to number of scalars!");

        // get time-step length
        const double dt = my::scatraparatimint_->Dt();

        // extract local values from the global vectors
        Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
        if (phinp == Teuchos::null) dserror("Cannot get state vector 'phinp'");
        DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, lm);

        my::EvalShapeFuncAndDerivsAtEleCenter();

        for (unsigned i = 0; i < updatemat.size(); i++)
        {
          const double csnp = my::funct_.Dot(my::ephinp_[i]);  // be careful, we assume k==i here
          updatemat[i]->Update(csnp, dt);
        }
      }

      break;
    }

    case SCATRA::get_material_internal_state:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // access the general material
        Teuchos::RCP<MAT::Material> material = ele->Material();
        Teuchos::RCP<Epetra_MultiVector> material_internal_state =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("material_internal_state");

        if (material->MaterialType() == INPAR::MAT::m_myocard)
        {
          Teuchos::RCP<MAT::Myocard> material =
              Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
          for (int k = 0; k < material_internal_state->NumVectors(); ++k)
          {
            int err = material_internal_state->ReplaceGlobalValue(
                ele->Id(), k, material->GetInternalState(k));
            if (err != 0) dserror("%i", err);
          }
        }
        params.set<Teuchos::RCP<Epetra_MultiVector>>(
            "material_internal_state", material_internal_state);
      }

      break;
    }

    case SCATRA::set_material_internal_state:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // access the general material
        Teuchos::RCP<MAT::Material> material = ele->Material();
        Teuchos::RCP<Epetra_Vector> material_internal_state_component =
            params.get<Teuchos::RCP<Epetra_Vector>>("material_internal_state_component");

        if (material->MaterialType() == INPAR::MAT::m_myocard)
        {
          Teuchos::RCP<MAT::Myocard> material =
              Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
          int k = params.get<int>("k");
          material->SetInternalState(k, (*material_internal_state_component)[ele->Id()]);
        }
      }
    }

    break;

    case SCATRA::get_material_ionic_currents:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
      {
        // access the general material
        Teuchos::RCP<MAT::Material> material = ele->Material();
        Teuchos::RCP<Epetra_MultiVector> material_ionic_currents =
            params.get<Teuchos::RCP<Epetra_MultiVector>>("material_ionic_currents");

        if (material->MaterialType() == INPAR::MAT::m_myocard)
        {
          Teuchos::RCP<MAT::Myocard> material =
              Teuchos::rcp_dynamic_cast<MAT::Myocard>(ele->Material());
          for (int k = 0; k < material_ionic_currents->NumVectors(); ++k)
          {
            int err = material_ionic_currents->ReplaceGlobalValue(
                ele->Id(), k, material->GetIonicCurrents(k));
            if (err != 0) dserror("%i", err);
          }
        }
        params.set<Teuchos::RCP<Epetra_MultiVector>>(
            "material_ionic_currents", material_ionic_currents);
      }

      break;
    }

    default:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch(action)

  return 0;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<DRT::Element::nurbs27>;
