/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element evaluation with ScaTra coupling

*----------------------------------------------------------------------*/

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_lib_discret.hpp"
#include "4C_membrane_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  pre-evaluate the element (public)                      sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::MembraneScatra<distype>::pre_evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la)
{
  if (la.Size() > 1)
  {
    // ask for the number of dofs of second dofset (scatra)
    const int numscal = discretization.NumDof(1, Nodes()[0]);

    if (la[1].Size() != Membrane<distype>::numnod_ * numscal)
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
          Teuchos::rcp(new std::vector<double>(la[1].lm_.size(), 0.0));

      Core::FE::ExtractMyValues(*scalarnp, *myscalar, la[1].lm_);

      // element vector for k-th scalar
      std::vector<Core::LinAlg::Matrix<Membrane<distype>::numnod_, 1>> elescalar(numscal);
      for (int k = 0; k < numscal; ++k)
      {
        for (int i = 0; i < Membrane<distype>::numnod_; ++i)
        {
          (elescalar.at(k))(i, 0) = myscalar->at(numscal * i + k);
        }
      }

      // number of gauss points
      const int numgp = (Membrane<distype>::intpoints_).nquad;

      // create vector of gauss point values to be set in params list
      Teuchos::RCP<std::vector<std::vector<double>>> gpscalar = Teuchos::rcp(
          new std::vector<std::vector<double>>(numgp, std::vector<double>(numscal, 0.0)));

      // allocate vector for shape functions and matrix for derivatives at gp
      Core::LinAlg::Matrix<Membrane<distype>::numnod_, 1> shapefcts(true);

      // loop over gauss points
      for (int gp = 0; gp < numgp; ++gp)
      {
        // get gauss points from integration rule
        double xi_gp = (Membrane<distype>::intpoints_).qxg[gp][0];
        double eta_gp = (Membrane<distype>::intpoints_).qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        Core::FE::shape_function_2D(shapefcts, xi_gp, eta_gp, Shape());

        // scalar at current gp
        std::vector<double> scalar_curr_gp(numscal, 0.0);

        for (int k = 0; k < numscal; ++k)
        {
          // identical shapefunctions for displacements and scalar fields
          scalar_curr_gp.at(k) = shapefcts.Dot(elescalar.at(k));
        }

        gpscalar->at(gp) = scalar_curr_gp;
      }

      // set scalar states at gp to params list
      params.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_scalar", gpscalar);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                          sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::MembraneScatra<distype>::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // in some cases we need to write/change some data before evaluating
  pre_evaluate(params, discretization, la);

  Membrane<distype>::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra);

  return 0;
}

template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::MembraneScatra<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
