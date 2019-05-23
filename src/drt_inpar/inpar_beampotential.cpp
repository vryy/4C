/*----------------------------------------------------------------------*/
/*!

\brief input parameter definitions for beam potential-based interactions

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_beampotential.H"
#include "inpar_beamcontact.H"
#include "inpar_structure.H"
#include "inpar_tsi.H"
#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::BEAMPOTENTIAL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  /* parameters for potential-based beam interaction */
  Teuchos::ParameterList& beampotential = list->sublist("BEAM POTENTIAL", false, "");

  setNumericStringParameter("POT_LAW_EXPONENT", "1.0",
      "negative(!) exponent(s) m_i of potential law Phi(r) = sum_i (k_i * r^(-m_i)).",
      &beampotential);
  setNumericStringParameter("POT_LAW_PREFACTOR", "0.0",
      "prefactor(s) k_i of potential law Phi(r) = sum_i (k_i * r^(-m_i)).", &beampotential);

  setStringToIntegralParameter<int>("BEAMPOTENTIAL_TYPE", "Surface",
      "Type of potential interaction: surface (default) or volume potential",
      tuple<std::string>("Surface", "surface", "Volume", "volume"),
      tuple<int>(beampot_surf, beampot_surf, beampot_vol, beampot_vol), &beampotential);

  setStringToIntegralParameter<int>("STRATEGY", "DoubleLengthSpecific_LargeSepApprox",
      "strategy to evaluate interaction potential: double/single length specific, "
      "small/large separation approximation, ...",
      tuple<std::string>("DoubleLengthSpecific_LargeSepApprox",
          "DoubleLengthSpecific_SmallSepApprox", "SingleLengthSpecific_SmallSepApprox",
          "SingleLengthSpecific_SmallSepApprox_Simple"),
      tuple<int>(strategy_doublelengthspec_largesepapprox, strategy_doublelengthspec_smallsepapprox,
          strategy_singlelengthspec_smallsepapprox,
          strategy_singlelengthspec_smallsepapprox_simple),
      &beampotential);

  DoubleParameter("CUTOFF_RADIUS", -1.0,
      "Neglect all potential contributions at separation larger"
      "than this cutoff radius",
      &beampotential);

  setStringToIntegralParameter<int>("REGULARIZATION_TYPE", "none",
      "Type of regularization applied to the force law",
      tuple<std::string>("linear_extrapolation", "constant_extrapolation", "None", "none"),
      tuple<int>(
          regularization_linear, regularization_constant, regularization_none, regularization_none),
      &beampotential);

  DoubleParameter("REGULARIZATION_SEPARATION", -1.0,
      "Use regularization of force law at separations "
      "smaller than this separation",
      &beampotential);

  IntParameter("NUM_INTEGRATION_SEGMENTS", 1,
      "Number of integration segments used per beam element", &beampotential);

  IntParameter(
      "NUM_GAUSSPOINTS", 10, "Number of Gauss points used per integration segment", &beampotential);

  setStringToIntegralParameter<int>("AUTOMATIC_DIFFERENTIATION", "No",
      "apply automatic differentiation via FAD?", yesnotuple, yesnovalue, &beampotential);

  setStringToIntegralParameter<int>("BEAMPOT_BTSOL", "No",
      "decide, whether potential-based interaction between beams and solids is considered",
      yesnotuple, yesnovalue, &beampotential);

  setStringToIntegralParameter<int>("BEAMPOT_BTSPH", "No",
      "decide, whether potential-based interaction between beams and spheres is considered",
      yesnotuple, yesnovalue, &beampotential);

  // enable octree search and determine type of bounding box (aabb = axis aligned, spbb = spherical)
  setStringToIntegralParameter<int>("BEAMPOT_OCTREE", "None",
      "octree and bounding box type for octree search routine",
      tuple<std::string>(
          "None", "none", "octree_axisaligned", "octree_cylorient", "octree_spherical"),
      tuple<int>(INPAR::BEAMCONTACT::boct_none, INPAR::BEAMCONTACT::boct_none,
          INPAR::BEAMCONTACT::boct_aabb, INPAR::BEAMCONTACT::boct_cobb,
          INPAR::BEAMCONTACT::boct_spbb),
      &beampotential);

  IntParameter("BEAMPOT_TREEDEPTH", 6, "max, tree depth of the octree", &beampotential);
  IntParameter(
      "BEAMPOT_BOXESINOCT", 8, "max number of bounding boxes in any leaf octant", &beampotential);

  /*------------------------------------------------------------------------*/
  /* parameters for visualization of potential-based beam interactions via vtk output at runtime */

  Teuchos::ParameterList& beampotential_vtk_sublist =
      beampotential.sublist("RUNTIME VTK OUTPUT", false, "");


  // whether to write vtk output for beam contact
  setStringToIntegralParameter<int>("VTK_OUTPUT_BEAM_POTENTIAL", "No",
      "write vtk output for potential-based beam interactions", yesnotuple, yesnovalue,
      &beampotential_vtk_sublist);

  // output interval regarding steps: write output every INTERVAL_STEPS steps
  IntParameter("INTERVAL_STEPS", -1, "write VTK output at runtime every INTERVAL_STEPS steps",
      &beampotential_vtk_sublist);

  // data format for written numeric data
  setStringToIntegralParameter<int>("OUTPUT_DATA_FORMAT", "binary",
      "data format for written numeric data",
      tuple<std::string>("binary", "Binary", "ascii", "ASCII"),
      tuple<int>(INPAR::BEAMPOTENTIAL::binary, INPAR::BEAMPOTENTIAL::binary,
          INPAR::BEAMPOTENTIAL::ascii, INPAR::BEAMPOTENTIAL::ascii),
      &beampotential_vtk_sublist);

  // whether to write output in every iteration of the nonlinear solver
  setStringToIntegralParameter<int>("EVERY_ITERATION", "No",
      "write output in every iteration of the nonlinear solver", yesnotuple, yesnovalue,
      &beampotential_vtk_sublist);

  // whether to write vtp output for forces
  setStringToIntegralParameter<int>("FORCES", "No", "write vtp output for forces", yesnotuple,
      yesnovalue, &beampotential_vtk_sublist);

  // whether to write vtp output for moments
  setStringToIntegralParameter<int>("MOMENTS", "No", "write vtp output for moments", yesnotuple,
      yesnovalue, &beampotential_vtk_sublist);

  // whether to write vtp output for forces/moments separately for each element pair
  setStringToIntegralParameter<int>("WRITE_FORCE_MOMENT_PER_ELEMENTPAIR", "No",
      "write vtp output for forces/moments separately for each element pair", yesnotuple,
      yesnovalue, &beampotential_vtk_sublist);
}

void INPAR::BEAMPOTENTIAL::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*-------------------------------------------------------------------*/
  // beam potential interaction: atom/charge density per unit length on LINE
  Teuchos::RCP<ConditionDefinition> rigidsphere_potential_charge =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT RIGIDSPHERE POTENTIAL CHARGE CONDITIONS",
          "RigidspherePotentialPointCharge", "Rigidsphere_Potential_Point_Charge",
          DRT::Condition::RigidspherePotential_PointCharge, false, DRT::Condition::Point));

  Teuchos::RCP<ConditionDefinition> beam_potential_line_charge =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE BEAM POTENTIAL CHARGE CONDITIONS",
          "BeamPotentialLineCharge", "Beam_Potential_Line_Charge_Density",
          DRT::Condition::BeamPotential_LineChargeDensity, false, DRT::Condition::Line));

  rigidsphere_potential_charge->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("POTLAW")));
  rigidsphere_potential_charge->AddComponent(
      Teuchos::rcp(new IntConditionComponent("potlaw", false, false)));
  rigidsphere_potential_charge->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  rigidsphere_potential_charge->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  rigidsphere_potential_charge->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  rigidsphere_potential_charge->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, true, true)));

  beam_potential_line_charge->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("POTLAW")));
  beam_potential_line_charge->AddComponent(
      Teuchos::rcp(new IntConditionComponent("potlaw", false, false)));
  beam_potential_line_charge->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  beam_potential_line_charge->AddComponent(
      Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  beam_potential_line_charge->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  beam_potential_line_charge->AddComponent(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, true, true)));

  condlist.push_back(rigidsphere_potential_charge);
  condlist.push_back(beam_potential_line_charge);
}
