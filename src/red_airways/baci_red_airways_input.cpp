/*---------------------------------------------------------------------*/
/*! \file

\brief Input for RedAirway, RedAcinus, RedInterAcinarDep, RedAirBloodScatra and
RedAirBloodScatraLine3 elements


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_global_data.H"
#include "baci_io_linedefinition.H"
#include "baci_mat_maxwell_0d_acinus.H"
#include "baci_red_airways_elementbase.H"

BACI_NAMESPACE_OPEN

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
| Read in the RED_AIRWAY elements                                       |
*-----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirway::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim != 3)
    dserror("Problem defined as %dd, but found Reduced dimensional AIRWAY element.", ndim);

  // Read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // Read the element type, the element specific variables and store them to airwayParams_
  linedef->ExtractString("TYPE", elemType_);
  if (elemType_ == "Resistive" || elemType_ == "InductoResistive" ||
      elemType_ == "ComplientResistive" || elemType_ == "RLC" || elemType_ == "ViscoElasticRLC" ||
      elemType_ == "ConvectiveViscoElasticRLC")
  {
    linedef->ExtractString("Resistance", resistance_);
    linedef->ExtractString("ElemSolvingType", elemsolvingType_);

    double Ew, tw, A, Ts, Phis, nu, velPow;
    int generation;
    linedef->ExtractDouble("PowerOfVelocityProfile", velPow);
    linedef->ExtractDouble("WallElasticity", Ew);
    linedef->ExtractDouble("PoissonsRatio", nu);
    linedef->ExtractDouble("ViscousTs", Ts);
    linedef->ExtractDouble("ViscousPhaseShift", Phis);
    linedef->ExtractDouble("WallThickness", tw);
    linedef->ExtractDouble("Area", A);
    linedef->ExtractInt("Generation", generation);

    if (linedef->HaveNamed("AirwayColl"))
    {
      double airwayColl;
      linedef->ExtractDouble("AirwayColl", airwayColl);
      double sc, so, Pcrit_c, Pcrit_o, open_init;
      linedef->ExtractDouble("S_Close", sc);
      linedef->ExtractDouble("S_Open", so);
      linedef->ExtractDouble("Pcrit_Open", Pcrit_o);
      linedef->ExtractDouble("Pcrit_Close", Pcrit_c);
      linedef->ExtractDouble("Open_Init", open_init);
      airwayParams_.airway_coll = airwayColl;
      airwayParams_.s_close = sc;
      airwayParams_.s_open = so;
      airwayParams_.p_crit_open = Pcrit_o;
      airwayParams_.p_crit_close = Pcrit_c;
      airwayParams_.open_init = open_init;
    }

    // Correct the velocity profile power
    // this is because the 2.0 is the minimum energy consumtive laminar profile
    if (velPow < 2.0) velPow = 2.0;
    airwayParams_.power_velocity_profile = velPow;
    airwayParams_.wall_elasticity = Ew;
    airwayParams_.poisson_ratio = nu;
    airwayParams_.wall_thickness = tw;
    airwayParams_.area = A;
    airwayParams_.viscous_Ts = Ts;
    airwayParams_.viscous_phase_shift = Phis;
    airwayParams_.generation = generation;
    if (linedef->HaveNamed("BranchLength"))
    {
      double l_branch = 0.0;
      linedef->ExtractDouble("BranchLength", l_branch);
      airwayParams_.branch_length = l_branch;
    }
  }
  else
  {
    dserror(
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
bool DRT::ELEMENTS::RedAcinus::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim != 3)
    dserror("Problem defined as %dd, but found Reduced dimensional ACINUS element.", ndim);

  // Read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // Read the element type, the element specific variables and store them to acinusParams_
  linedef->ExtractString("TYPE", elemType_);
  if (elemType_ == "NeoHookean" || elemType_ == "Exponential" || elemType_ == "DoubleExponential" ||
      elemType_ == "VolumetricOgden")
  {
    double acinusVol, alveolarDuctVol, A;
    const int generation = -1;
    linedef->ExtractDouble("AcinusVolume", acinusVol);
    linedef->ExtractDouble("AlveolarDuctVolume", alveolarDuctVol);
    linedef->ExtractDouble("Area", A);

    acinusParams_.volume_relaxed = acinusVol;
    acinusParams_.alveolar_duct_volume = alveolarDuctVol;
    acinusParams_.area = A;
    acinusParams_.volume_init = acinusVol;
    acinusParams_.generation = generation;

    // Setup material, calls overloaded function Setup(linedef) for each Maxwell_0d_acinus material
    Teuchos::RCP<MAT::Material> mat = Material();
    Teuchos::RCP<MAT::Maxwell_0d_acinus> acinus_mat =
        Teuchos::rcp_dynamic_cast<MAT::Maxwell_0d_acinus>(Material());
    acinus_mat->Setup(linedef);
  }
  else
  {
    dserror(
        "Reading type of RED_ACINUS element failed. Possible types: NeoHookean/ Exponential"
        "/ DoubleExponential/ VolumetricOgden");
    exit(1);
  }

  return true;
}


/*----------------------------------------------------------------------*
| Read in the RED_ACINAR_INTER_DEP elements                             |
*-----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedInterAcinarDep::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim != 3)
    dserror(
        "Problem defined as %dd, but found Reduced dimensional INTER ACINAR DEPENDENCE element.",
        ndim);

  // set generation
  const int generation = -2;
  generation_ = generation;

  // Read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);


  return true;
}


/*----------------------------------------------------------------------*
| Read in the Scatra elements                                           |
*-----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirBloodScatra::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim != 3)
    dserror("Problem defined as %dd, but found Reduced dimensional Scatra element.", ndim);

  // read number of material model
  const int generation = -2;
  generation_ = generation;

  double diff = 0.0;
  linedef->ExtractDouble("DiffusionCoefficient", diff);
  elemParams_["DiffusionCoefficient"] = diff;

  double th = 0.0;
  linedef->ExtractDouble("WallThickness", th);
  elemParams_["WallThickness"] = th;

  double percDiffArea = 0.0;
  linedef->ExtractDouble("PercentageOfDiffusionArea", percDiffArea);
  elemParams_["PercentageOfDiffusionArea"] = percDiffArea;

  return true;
}


/*----------------------------------------------------------------------*
| Read in the Scatra elements                                           |
*-----------------------------------------------------------------------*/
bool DRT::ELEMENTS::RedAirBloodScatraLine3::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  const int ndim = DRT::Problem::Instance()->NDim();
  if (ndim != 3)
    dserror("Problem defined as %dd, but found Reduced dimensional Scatra element.", ndim);

  double diff = 0.0;
  linedef->ExtractDouble("DiffusionCoefficient", diff);
  elemParams_["DiffusionCoefficient"] = diff;

  double th = 0.0;
  linedef->ExtractDouble("WallThickness", th);
  elemParams_["WallThickness"] = th;

  double percDiffArea = 0.0;
  linedef->ExtractDouble("PercentageOfDiffusionArea", percDiffArea);
  elemParams_["PercentageOfDiffusionArea"] = percDiffArea;

  return true;
}

BACI_NAMESPACE_CLOSE
