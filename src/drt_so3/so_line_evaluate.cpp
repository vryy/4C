/*----------------------------------------------------------------------*/
/*! \file

\brief line element

\level 1

\maintainer Christoph Meier
*----------------------------------------------------------------------*/

#include "so_line.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_elements_paramsinterface.H"


/*-----------------------------------------------------------------------*
 * Integrate a Line Neumann boundary condition (public)         gee 04/08|
 * ----------------------------------------------------------------------*/
int DRT::ELEMENTS::StructuralLine::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // set the interface ptr in the parent element
  ParentElement()->SetParamsInterfacePtr(params);
  // get type of condition
  enum LoadType
  {
    neum_none,
    neum_live
  };
  LoadType ltype;

  // spatial or material configuration depends on the type of load
  // currently only material frame used
  // enum Configuration
  //{
  //  config_none,
  //  config_material,
  //  config_spatial,
  //  config_both
  //};
  // Configuration config = config_none;

  const std::string* type = condition.Get<std::string>("type");
  if (*type == "neum_live")
  {
    ltype = neum_live;
    // config = config_material;
  }
  else
    dserror("Unknown type of LineNeumann condition");

  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");
  const std::vector<int>* spa_func = condition.Get<std::vector<int>>("funct");

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  double time = -1.0;
  if (ParentElement()->IsParamsInterface())
    time = ParentElement()->ParamsInterfacePtr()->GetTotalTime();
  else
    time = params.get("total time", -1.0);

  const int numdim = 3;

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < numdim)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = numdim; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // element geometry update - currently only material configuration
  const int numnode = NumNode();
  LINALG::SerialDenseMatrix x(numnode, numdim);
  MaterialConfiguration(x);

  // integration parameters
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule_);
  LINALG::SerialDenseVector shapefcts(numnode);
  LINALG::SerialDenseMatrix deriv(1, numnode);
  const DRT::Element::DiscretizationType shape = Shape();

  // integration
  for (int gp = 0; gp < intpoints.nquad; ++gp)
  {
    // get shape functions and derivatives of element surface
    const double e = intpoints.qxg[gp][0];
    DRT::UTILS::shape_function_1D(shapefcts, e, shape);
    DRT::UTILS::shape_function_1D_deriv1(deriv, e, shape);
    switch (ltype)
    {
      case neum_live:
      {  // uniform load on reference configuration

        double dL;
        LineIntegration(dL, x, deriv);

        // loop the dofs of a node
        for (int i = 0; i < numdim; ++i)
        {
          if ((*onoff)[i])  // is this dof activated?
          {
            // factor given by spatial function
            const int functnum = (spa_func) ? (*spa_func)[i] : -1;
            double functfac = 1.0;

            if (functnum > 0)
            {
              // calculate reference position of GP
              LINALG::SerialDenseMatrix gp_coord(1, numdim);
              gp_coord.Multiply('T', 'N', 1.0, shapefcts, x, 0.0);

              // write coordinates in another datatype
              double gp_coord2[numdim];
              for (int k = 0; k < numdim; k++) gp_coord2[k] = gp_coord(0, k);
              const double* coordgpref = &gp_coord2[0];  // needed for function evaluation

              // evaluate function at current gauss point
              functfac =
                  DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(i, coordgpref, time);
            }

            const double fac = (*val)[i] * intpoints.qwgt[gp] * dL * functfac;
            for (int node = 0; node < numnode; ++node)
            {
              elevec1[node * numdim + i] += shapefcts[node] * fac;
            }
          }
        }

        break;
      }

      default:
        dserror("Unknown type of LineNeumann load");
        break;
    }
  }

  return 0;
}

/*----------------------------------------------------------------------*
 *  (private)                                                  gee 04/08|
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::StructuralLine::LineIntegration(
    double& dL, const LINALG::SerialDenseMatrix& x, const LINALG::SerialDenseMatrix& deriv)
{
  // compute dXYZ / drs
  LINALG::SerialDenseMatrix dxyzdrs(1, 3);
  dxyzdrs.Multiply('N', 'N', 1.0, deriv, x, 0.0);
  dL = 0.0;
  for (int i = 0; i < 3; ++i) dL += dxyzdrs(0, i) * dxyzdrs(0, i);
  dL = sqrt(dL);
  return;
}
