/*! \file
 *
\brief Nonlinear Shell 7-Parameter Model Finite Element line evaluation

\level 3
*/

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_integration.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_shell7p_line.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

int DRT::ELEMENTS::Shell7pLine::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& dof_index_array, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // set the interface ptr in the parent element
  ParentElement()->SetParamsInterfacePtr(params);

  // we need the displacement at the previous step
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'displacement'");
  std::vector<double> displacements(dof_index_array.size());
  CORE::FE::ExtractMyValues(*disp, displacements, dof_index_array);

  // get values and switches from the condition
  const auto* onoff = &condition.parameters().Get<std::vector<int>>("onoff");
  const auto* val = &condition.parameters().Get<std::vector<double>>("val");
  const auto* spa_func = &condition.parameters().Get<std::vector<int>>("funct");

  // time curve buisiness
  // find out whether we will use a time curve
  double time = -1.0;
  if (ParentElement()->IsParamsInterface())
    time = ParentElement()->StrParamsInterface().GetTotalTime();
  else
    time = params.get<double>("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < node_dof_)
    FOUR_C_THROW(
        "Fewer functions or curves defined than the element's nodal degree of freedoms (6).");

  for (int checkdof = num_dim_; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
    {
      FOUR_C_THROW(
          "Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
    }
  }

  // element geometry update - currently only material configuration
  const int numnode = NumNode();
  CORE::LINALG::SerialDenseMatrix x(numnode, num_dim_);
  MaterialConfiguration(x);

  // integration parameters
  const CORE::FE::IntegrationPoints1D intpoints(CORE::FE::GaussRule1D::line_2point);
  CORE::LINALG::SerialDenseVector shape_functions(numnode);
  CORE::LINALG::SerialDenseMatrix derivatives(1, numnode);
  const CORE::FE::CellType shape = Shape();

  // integration
  for (int gp = 0; gp < intpoints.NumPoints(); ++gp)
  {
    // get shape functions and derivatives of element surface
    const double e = intpoints.qxg[gp][0];
    CORE::FE::shape_function_1D(shape_functions, e, shape);
    CORE::FE::shape_function_1D_deriv1(derivatives, e, shape);

    // covariant basis vectors and metric of shell body
    // g1,g2,g3 stored in Jacobian matrix  = (g1,g2,g3)
    double dL;
    LineIntegration(dL, x, derivatives);
    std::vector<double> a;
    // loop through the dofs of a node
    for (int i = 0; i < num_dim_; ++i)
    {
      // check if this dof is activated
      if ((*onoff)[i])
      {
        // factor given by spatial function
        const int functnum = (spa_func) ? (*spa_func)[i] : -1;
        double functfac = 1.0;

        if (functnum > 0)
        {
          // calculate reference position of gaussian point
          CORE::LINALG::SerialDenseMatrix gp_coord(1, num_dim_);
          gp_coord.multiply(Teuchos::TRANS, Teuchos::NO_TRANS, 1.0, shape_functions, x, 0.0);

          // write coordinates in another datatype
          double gp_coord2[num_dim_];
          for (int k = 0; k < num_dim_; k++) gp_coord2[k] = gp_coord(0, k);
          const double* coordgpref = gp_coord2;  // needed for function evaluation

          // evaluate function at current gauss point
          functfac = GLOBAL::Problem::Instance()
                         ->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functnum - 1)
                         .Evaluate(coordgpref, time, i);
        }

        const double fac = (*val)[i] * intpoints.qwgt[gp] * dL * functfac;

        for (int node = 0; node < numnode; ++node)
        {
          elevec1[node * node_dof_ + i] += shape_functions[node] * fac;
        }
      }
    }
  }

  return 0;
}


void DRT::ELEMENTS::Shell7pLine::LineIntegration(double& dL,
    const CORE::LINALG::SerialDenseMatrix& x, const CORE::LINALG::SerialDenseMatrix& derivatives)
{
  // compute dXYZ / drs
  CORE::LINALG::SerialDenseMatrix dxyzdrs(1, num_dim_);
  dxyzdrs.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, derivatives, x, 0.0);
  dL = 0.0;

  for (int i = 0; i < 3; ++i) dL += dxyzdrs(0, i) * dxyzdrs(0, i);

  dL = sqrt(dL);
}
FOUR_C_NAMESPACE_CLOSE
