/*----------------------------------------------------------------------*/
/*! \file

\level 3


\brief Nonlinear Membrane Finite Element evaluation with ScaTra coupling

*----------------------------------------------------------------------*/

#include "membrane_scatra.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

/*----------------------------------------------------------------------*
 |  pre-evaluate the element (public)                      sfuchs 05/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::MembraneScatra<distype>::PreEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  if (la.Size() > 1)
  {
    // ask for the number of dofs of second dofset (scatra)
    const int numscal = discretization.NumDof(1, Nodes()[0]);

    if (la[1].Size() != Membrane<distype>::numnod_ * numscal)
      dserror("location vector length does not match!");

    // name of scalarfield
    std::string scalarfield = "scalarfield";

    if (discretization.HasState(1, scalarfield))
    {
      // get the scalar state
      Teuchos::RCP<const Epetra_Vector> scalarnp = discretization.GetState(1, scalarfield);

      if (scalarnp == Teuchos::null) dserror("can not get state vector %s", scalarfield.c_str());

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double>> myscalar =
          Teuchos::rcp(new std::vector<double>(la[1].lm_.size(), 0.0));

      DRT::UTILS::ExtractMyValues(*scalarnp, *myscalar, la[1].lm_);

      // element vector for k-th scalar
      std::vector<LINALG::Matrix<Membrane<distype>::numnod_, 1>> elescalar(numscal);
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
      LINALG::Matrix<Membrane<distype>::numnod_, 1> shapefcts(true);

      // loop over gauss points
      for (int gp = 0; gp < numgp; ++gp)
      {
        // get gauss points from integration rule
        double xi_gp = (Membrane<distype>::intpoints_).qxg[gp][0];
        double eta_gp = (Membrane<distype>::intpoints_).qxg[gp][1];

        // get shape functions and derivatives in the plane of the element
        DRT::UTILS::shape_function_2D(shapefcts, xi_gp, eta_gp, Shape());

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
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::MembraneScatra<distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // in some cases we need to write/change some data before evaluating
  PreEvaluate(params, discretization, la);

  Membrane<distype>::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra);

  return 0;
}

template class DRT::ELEMENTS::MembraneScatra<DRT::Element::tri3>;
template class DRT::ELEMENTS::MembraneScatra<DRT::Element::tri6>;
template class DRT::ELEMENTS::MembraneScatra<DRT::Element::quad4>;
template class DRT::ELEMENTS::MembraneScatra<DRT::Element::quad9>;
