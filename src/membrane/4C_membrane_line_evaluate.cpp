/*----------------------------------------------------------------------*/
/*! \file
\brief

\level 3


\brief Nonlinear Membrane Finite Element line evaluation

*----------------------------------------------------------------------*/
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_global_data.hpp"
#include "4C_membrane.hpp"
#include "4C_structure_new_elements_paramsinterface.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)   fbraeu 06/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::MembraneLine<distype>::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  // set params interface pointer in the parent element
  parent_element()->set_params_interface_ptr(params);

  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live
  };
  LoadType ltype;

  const std::string& type = condition.parameters().get<std::string>("type");
  if (type == "neum_live")
  {
    ltype = neum_live;
  }
  else
    FOUR_C_THROW("Unknown type of LineNeumann condition");

  // get values and switches from the condition
  const auto* onoff = &condition.parameters().get<std::vector<int>>("onoff");
  const auto* val = &condition.parameters().get<std::vector<double>>("val");
  const auto* spa_func = &condition.parameters().get<std::vector<int>>("funct");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  double time = -1.0;
  if (parent_element()->is_params_interface())
    time = parent_element()->params_interface_ptr()->get_total_time();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < noddof_)
    FOUR_C_THROW("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = noddof_; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      FOUR_C_THROW(
          "Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // element geometry update - currently only material configuration
  Core::LinAlg::Matrix<numnod_line_, noddof_> x(true);
  for (int i = 0; i < numnod_line_; ++i)
  {
    x(i, 0) = nodes()[i]->x()[0];
    x(i, 1) = nodes()[i]->x()[1];
    x(i, 2) = nodes()[i]->x()[2];
  }

  // allocate vector for shape functions and matrix for derivatives at gp
  Core::LinAlg::Matrix<numnod_line_, 1> shapefcts(true);
  Core::LinAlg::Matrix<1, numnod_line_> derivs(true);

  // integration
  for (int gp = 0; gp < intpointsline_.nquad; ++gp)
  {
    // get gausspoints from integration rule
    double xi_gp = intpointsline_.qxg[gp][0];

    // get gauss weight at current gp
    double gpweight = intpointsline_.qwgt[gp];

    // get shape functions and derivatives in the plane of the element
    Core::FE::shape_function_1D(shapefcts, xi_gp, shape());
    Core::FE::shape_function_1D_deriv1(derivs, xi_gp, shape());

    switch (ltype)
    {
      case neum_live:
      {
        // uniform load on reference configuration

        // compute dXYZ / dr
        Core::LinAlg::Matrix<noddof_, 1> dxyzdr(true);
        dxyzdr.multiply_tt(1.0, x, derivs, 0.0);
        // compute line increment dL
        double dL;
        dL = 0.0;
        for (int i = 0; i < 3; ++i)
        {
          dL += dxyzdr(i) * dxyzdr(i);
        }
        dL = sqrt(dL);

        // loop the dofs of a node
        for (int i = 0; i < noddof_; ++i)
        {
          if ((*onoff)[i])  // is this dof activated?
          {
            // factor given by spatial function
            const int functnum = (spa_func) ? (*spa_func)[i] : -1;
            double functfac = 1.0;

            if (functnum > 0)
            {
              // calculate reference position of GP
              Core::LinAlg::Matrix<noddof_, 1> gp_coord(true);
              gp_coord.multiply_tn(1.0, x, shapefcts, 0.0);

              // write coordinates in another datatype
              double gp_coord2[noddof_];
              for (int k = 0; k < noddof_; k++) gp_coord2[k] = gp_coord(k, 0);
              const double* coordgpref = gp_coord2;  // needed for function evaluation

              // evaluate function at current gauss point
              functfac = Global::Problem::instance()
                             ->function_by_id<Core::UTILS::FunctionOfSpaceTime>(functnum - 1)
                             .evaluate(coordgpref, time, i);
            }

            const double fac = (*val)[i] * gpweight * dL * functfac;
            for (int node = 0; node < numnod_line_; ++node)
            {
              elevec1[noddof_ * node + i] += shapefcts(node) * fac;
            }
          }
        }
        break;
      }

      default:
        FOUR_C_THROW("Unknown type of LineNeumann load");
        break;
    }
  }

  return 0;
}

template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::MembraneLine<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
