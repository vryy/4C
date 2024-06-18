/*---------------------------------------------------------------------*/
/*! \file

\brief Input for RedAirway, RedAcinus, RedInterAcinarDep, RedAirBloodScatra and
RedAirBloodScatraLine3 elements


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_maxwell_0d_acinus.hpp"
#include "4C_red_airways_elementbase.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
| Read in the RED_AIRWAY elements                                       |
*-----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedAirway::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  const int ndim = Global::Problem::Instance()->NDim();
  if (ndim != 3)
    FOUR_C_THROW("Problem defined as %dd, but found Reduced dimensional AIRWAY element.", ndim);

  // Read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // Read the element type, the element specific variables and store them to airwayParams_
  linedef->extract_string("TYPE", elem_type_);
  if (elem_type_ == "Resistive" || elem_type_ == "InductoResistive" ||
      elem_type_ == "ComplientResistive" || elem_type_ == "RLC" ||
      elem_type_ == "ViscoElasticRLC" || elem_type_ == "ConvectiveViscoElasticRLC")
  {
    linedef->extract_string("Resistance", resistance_);
    linedef->extract_string("ElemSolvingType", elemsolving_type_);

    double Ew, tw, A, Ts, Phis, nu, velPow;
    int generation;
    linedef->extract_double("PowerOfVelocityProfile", velPow);
    linedef->extract_double("WallElasticity", Ew);
    linedef->extract_double("PoissonsRatio", nu);
    linedef->extract_double("ViscousTs", Ts);
    linedef->extract_double("ViscousPhaseShift", Phis);
    linedef->extract_double("WallThickness", tw);
    linedef->extract_double("Area", A);
    linedef->extract_int("Generation", generation);

    if (linedef->has_named("AirwayColl"))
    {
      double airwayColl;
      linedef->extract_double("AirwayColl", airwayColl);
      double sc, so, Pcrit_c, Pcrit_o, open_init;
      linedef->extract_double("S_Close", sc);
      linedef->extract_double("S_Open", so);
      linedef->extract_double("Pcrit_Open", Pcrit_o);
      linedef->extract_double("Pcrit_Close", Pcrit_c);
      linedef->extract_double("Open_Init", open_init);
      airway_params_.airway_coll = airwayColl;
      airway_params_.s_close = sc;
      airway_params_.s_open = so;
      airway_params_.p_crit_open = Pcrit_o;
      airway_params_.p_crit_close = Pcrit_c;
      airway_params_.open_init = open_init;
    }

    // Correct the velocity profile power
    // this is because the 2.0 is the minimum energy consumtive laminar profile
    if (velPow < 2.0) velPow = 2.0;
    airway_params_.power_velocity_profile = velPow;
    airway_params_.wall_elasticity = Ew;
    airway_params_.poisson_ratio = nu;
    airway_params_.wall_thickness = tw;
    airway_params_.area = A;
    airway_params_.viscous_Ts = Ts;
    airway_params_.viscous_phase_shift = Phis;
    airway_params_.generation = generation;
    if (linedef->has_named("BranchLength"))
    {
      double l_branch = 0.0;
      linedef->extract_double("BranchLength", l_branch);
      airway_params_.branch_length = l_branch;
    }
  }
  else
  {
    FOUR_C_THROW(
        "Reading type of RED_AIRWAY element failed. Possible types: ComplientResistive/"
        "PoiseuilleResistive/TurbulentPoiseuilleResistive/InductoResistive/RLC/ViscoElasticRLC/"
        "ConvectiveViscoElasticRLC");
    exit(1);
  }

  return true;
}


/*----------------------------------------------------------------------*
| Read in the RED_ACINUS elements                                       |
*-----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedAcinus::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  const int ndim = Global::Problem::Instance()->NDim();
  if (ndim != 3)
    FOUR_C_THROW("Problem defined as %dd, but found Reduced dimensional ACINUS element.", ndim);

  // Read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));

  // Read the element type, the element specific variables and store them to acinusParams_
  linedef->extract_string("TYPE", elem_type_);
  if (elem_type_ == "NeoHookean" || elem_type_ == "Exponential" ||
      elem_type_ == "DoubleExponential" || elem_type_ == "VolumetricOgden")
  {
    double acinusVol, alveolarDuctVol, A;
    const int generation = -1;
    linedef->extract_double("AcinusVolume", acinusVol);
    linedef->extract_double("AlveolarDuctVolume", alveolarDuctVol);
    linedef->extract_double("Area", A);

    acinus_params_.volume_relaxed = acinusVol;
    acinus_params_.alveolar_duct_volume = alveolarDuctVol;
    acinus_params_.area = A;
    acinus_params_.volume_init = acinusVol;
    acinus_params_.generation = generation;

    // Setup material, calls overloaded function setup(linedef) for each Maxwell_0d_acinus material
    Teuchos::RCP<Core::Mat::Material> mat = Material();
    Teuchos::RCP<Mat::Maxwell0dAcinus> acinus_mat =
        Teuchos::rcp_dynamic_cast<Mat::Maxwell0dAcinus>(Material());
    acinus_mat->setup(linedef);
  }
  else
  {
    FOUR_C_THROW(
        "Reading type of RED_ACINUS element failed. Possible types: NeoHookean/ Exponential"
        "/ DoubleExponential/ VolumetricOgden");
    exit(1);
  }

  return true;
}


/*----------------------------------------------------------------------*
| Read in the RED_ACINAR_INTER_DEP elements                             |
*-----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedInterAcinarDep::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  const int ndim = Global::Problem::Instance()->NDim();
  if (ndim != 3)
    FOUR_C_THROW(
        "Problem defined as %dd, but found Reduced dimensional INTER ACINAR DEPENDENCE element.",
        ndim);

  // set generation
  const int generation = -2;
  generation_ = generation;

  // Read number of material model
  int material = 0;
  linedef->extract_int("MAT", material);
  SetMaterial(0, Mat::Factory(material));


  return true;
}


/*----------------------------------------------------------------------*
| Read in the Scatra elements                                           |
*-----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedAirBloodScatra::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  const int ndim = Global::Problem::Instance()->NDim();
  if (ndim != 3)
    FOUR_C_THROW("Problem defined as %dd, but found Reduced dimensional Scatra element.", ndim);

  // read number of material model
  const int generation = -2;
  generation_ = generation;

  double diff = 0.0;
  linedef->extract_double("DiffusionCoefficient", diff);
  elem_params_["DiffusionCoefficient"] = diff;

  double th = 0.0;
  linedef->extract_double("WallThickness", th);
  elem_params_["WallThickness"] = th;

  double percDiffArea = 0.0;
  linedef->extract_double("PercentageOfDiffusionArea", percDiffArea);
  elem_params_["PercentageOfDiffusionArea"] = percDiffArea;

  return true;
}


/*----------------------------------------------------------------------*
| Read in the Scatra elements                                           |
*-----------------------------------------------------------------------*/
bool Discret::ELEMENTS::RedAirBloodScatraLine3::ReadElement(
    const std::string& eletype, const std::string& distype, Input::LineDefinition* linedef)
{
  const int ndim = Global::Problem::Instance()->NDim();
  if (ndim != 3)
    FOUR_C_THROW("Problem defined as %dd, but found Reduced dimensional Scatra element.", ndim);

  double diff = 0.0;
  linedef->extract_double("DiffusionCoefficient", diff);
  elem_params_["DiffusionCoefficient"] = diff;

  double th = 0.0;
  linedef->extract_double("WallThickness", th);
  elem_params_["WallThickness"] = th;

  double percDiffArea = 0.0;
  linedef->extract_double("PercentageOfDiffusionArea", percDiffArea);
  elem_params_["PercentageOfDiffusionArea"] = percDiffArea;

  return true;
}

FOUR_C_NAMESPACE_CLOSE
