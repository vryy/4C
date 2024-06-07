/*! \file

\brief PreEvaluator of Shell7p-ScaTra elements

\level 1
*/


#include "4C_shell7p_ele_scatra_preevaluator.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_shell7p_ele_calc_lib.hpp"

FOUR_C_NAMESPACE_OPEN


void Discret::ELEMENTS::Shell::PreEvaluateScatraByElement(Core::Elements::Element& ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization,
    Core::Elements::Element::LocationArray& dof_index_array)
{
  switch (ele.Shape())
  {
    case Core::FE::CellType::quad4:
      return PreEvaluateScatra<Core::FE::CellType::quad4>(
          ele, params, discretization, dof_index_array);
    case Core::FE::CellType::quad8:
      return PreEvaluateScatra<Core::FE::CellType::quad8>(
          ele, params, discretization, dof_index_array);
    case Core::FE::CellType::quad9:
      return PreEvaluateScatra<Core::FE::CellType::quad9>(
          ele, params, discretization, dof_index_array);
    case Core::FE::CellType::tri3:
      return PreEvaluateScatra<Core::FE::CellType::tri3>(
          ele, params, discretization, dof_index_array);
    case Core::FE::CellType::tri6:
      return PreEvaluateScatra<Core::FE::CellType::tri6>(
          ele, params, discretization, dof_index_array);
    default:
      FOUR_C_THROW(
          "The discretization type you are trying to pre-evaluate for shell7p scatra is not yet "
          "implemented.");
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Shell::PreEvaluateScatra(Core::Elements::Element& ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization,
    Core::Elements::Element::LocationArray& dof_index_array)
{
  Core::FE::IntegrationPoints2D intpoints_midsurface_ =
      CreateGaussIntegrationPoints<distype>(get_gauss_rule<distype>());

  if (dof_index_array.Size() > 1)
  {
    // ask for the number of dofs of second dofset (scatra)
    const int numscal = discretization.NumDof(1, ele.Nodes()[0]);

    if (dof_index_array[1].Size() != Shell::DETAIL::num_node<distype> * numscal)
      FOUR_C_THROW("location vector length does not match!");

    // name of scalarfield
    std::string scalarfield = "scalarfield";

    if (discretization.HasState(1, scalarfield))
    {
      // get the scalar state
      Teuchos::RCP<const Epetra_Vector> scalarnp = discretization.GetState(1, scalarfield);

      if (scalarnp == Teuchos::null)
        FOUR_C_THROW("can not get state vector %s", scalarfield.c_str());

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double>> myscalar =
          Teuchos::rcp(new std::vector<double>(dof_index_array[1].lm_.size(), 0.0));

      Core::FE::ExtractMyValues(*scalarnp, *myscalar, dof_index_array[1].lm_);

      // element vector for k-th scalar
      std::vector<Core::LinAlg::Matrix<Shell::DETAIL::num_node<distype>, 1>> elescalar(numscal);
      for (int k = 0; k < numscal; ++k)
      {
        for (int i = 0; i < Shell::DETAIL::num_node<distype>; ++i)
        {
          (elescalar.at(k))(i, 0) = myscalar->at(numscal * i + k);
        }
      }

      // create vector of gauss point values to be set in params list
      Teuchos::RCP<std::vector<std::vector<double>>> gpscalar =
          Teuchos::rcp(new std::vector<std::vector<double>>(
              intpoints_midsurface_.NumPoints(), std::vector<double>(numscal, 0.0)));

      // allocate vector for shape functions and matrix for derivatives at gp
      Core::LinAlg::Matrix<Shell::DETAIL::num_node<distype>, 1> shapefunctions(true);

      // loop over gauss points
      for (int gp = 0; gp < intpoints_midsurface_.NumPoints(); ++gp)
      {
        // get gauss points from integration rule
        double xi_gp = intpoints_midsurface_.qxg[gp][0];
        double eta_gp = intpoints_midsurface_.qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        Core::FE::shape_function_2D(shapefunctions, xi_gp, eta_gp, distype);

        // scalar at current gp
        std::vector<double> scalar_curr_gp(numscal, 0.0);

        for (int k = 0; k < numscal; ++k)
        {
          // identical shapefunctions for displacements and scalar fields
          scalar_curr_gp.at(k) = shapefunctions.Dot(elescalar.at(k));
        }

        gpscalar->at(gp) = scalar_curr_gp;
      }

      // set scalar states at gp to params list
      params.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc", gpscalar);
    }
  }
  Core::LinAlg::Matrix<2, 1> center(true);
  params.set("elecenter_coords_ref", center);
}

FOUR_C_NAMESPACE_CLOSE
