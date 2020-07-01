/*-----------------------------------------------------------*/
/*! \file
\brief input parameter for Brownian dynamics simulation


\level 2

*/
/*-----------------------------------------------------------*/


#include "inpar_browniandyn.H"

#include "drt_validparameters.H"


void INPAR::BROWNIANDYN::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& browniandyn_list = list->sublist("BROWNIAN DYNAMICS", false, "");

  setStringToIntegralParameter<int>("BROWNDYNPROB", "No", "switch Brownian dynamics on/off",
      yesnotuple, yesnovalue, &browniandyn_list);



  // Reading double parameter for viscosity of background fluid
  DoubleParameter("VISCOSITY", 0.0, "viscosity", &browniandyn_list);

  // Reading double parameter for thermal energy in background fluid (temperature * Boltzmann
  // constant)
  DoubleParameter("KT", 0.0, "thermal energy", &browniandyn_list);

  // cutoff for random forces, which determines the maximal value
  DoubleParameter("MAXRANDFORCE", -1.0,
      "Any random force beyond MAXRANDFORCE*(standard dev.) will be omitted and redrawn. "
      "-1.0 means no bounds.'",
      &browniandyn_list);

  // time interval in which random numbers are constant
  DoubleParameter("TIMESTEP", -1.0,
      "Within this time interval the random numbers remain constant. -1.0 ", &browniandyn_list);

  // the way how damping coefficient values for beams are specified
  setStringToIntegralParameter<int>("BEAMS_DAMPING_COEFF_SPECIFIED_VIA", "cylinder_geometry_approx",
      "In which way are damping coefficient values for beams specified?",
      tuple<std::string>(
          "cylinder_geometry_approx", "Cylinder_geometry_approx", "input_file", "Input_file"),
      tuple<int>(INPAR::BROWNIANDYN::cylinder_geometry_approx,
          INPAR::BROWNIANDYN::cylinder_geometry_approx, INPAR::BROWNIANDYN::input_file,
          INPAR::BROWNIANDYN::input_file),
      &browniandyn_list);

  // values for damping coefficients of beams if they are specified via input file
  // (per unit length, NOT yet multiplied by fluid viscosity)
  StringParameter("BEAMS_DAMPING_COEFF_PER_UNITLENGTH", "0.0 0.0 0.0",
      "values for beam damping coefficients (per unit length and NOT yet multiplied by fluid "
      "viscosity): "
      "translational perpendicular/parallel to beam axis, rotational around axis",
      &browniandyn_list);
}
