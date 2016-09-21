/*----------------------------------------------------------------------*/
/*!
\file inpar_potential.cpp

\brief Input parameters for potential

\level 3

\maintainer Thomas Klöppel

*/
/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_potential.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::POTENTIAL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& interaction_potential = list->sublist("INTERACTION POTENTIAL",false,"");

  // read if surfaces , volumes or both including fluid should be considered
  setStringToIntegralParameter<int>("POTENTIAL_TYPE","Surface","Type of interaction potential",
                                tuple<std::string>("Surface",
                                                   "Volume",
                                                   "Surfacevolume",
                                                   "Surface_fsi",
                                                   "Volume_fsi",
                                                   "Surfacevolume_fsi"),
                                tuple<int>(
                                   potential_surface,
                                   potential_volume,
                                   potential_surfacevolume,
                                   potential_surface_fsi,
                                   potential_volume_fsi,
                                   potential_surfacevolume_fsi),
                                &interaction_potential);

  // approximation method
  setStringToIntegralParameter<int>("APPROXIMATION_TYPE","None","Type of approximation",
                                tuple<std::string>("None",
                                                   "Surface_approx",
                                                   "Point_approx"),
                                tuple<int>(
                                           approximation_none,
                                           approximation_surface,
                                           approximation_point),
                                &interaction_potential);

  // switches on the analytical solution computation for two van der waals spheres or membranes
  setStringToIntegralParameter<int>("ANALYTICALSOLUTION","None", "Type of analytical solution"
                                 "computes analytical solutions for two Van der Waals spheres or membranes",
                                 tuple<std::string>("None",
                                                    "Sphere",
                                                    "Membrane"),
                                 tuple<int>(
                                            solution_none,
                                            solution_sphere,
                                            solution_membrane),
                                            &interaction_potential);
  // use 2D integration for pseudo 3D
  setStringToIntegralParameter<int>("PSEUDO3D","no",
                                     "use 2D integration for pseudo 3D",
                                     yesnotuple,yesnovalue,&interaction_potential);

  // radius of can der Waals spheres for analytical testing
  DoubleParameter(  "VDW_RADIUS",0.0,
                    "radius of van der Waals spheres",
                    &interaction_potential);

  // thickness of sphericael Waals membranes for analytical testing
  DoubleParameter(  "THICKNESS",0.0,
                    "membrane thickness",
                    &interaction_potential);

  // number of atoms or molecules offset
  DoubleParameter(  "N_OFFSET",0.0,
                    "number of atoms or molecules offset",
                    &interaction_potential);
}



void INPAR::POTENTIAL::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // Lennard Jones potential volume
  Teuchos::RCP<ConditionDefinition> lj_potential_volume =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential_Volume",
                                         DRT::Condition::LJ_Potential_Volume,
                                         true,
                                         DRT::Condition::Volume));

  lj_potential_volume->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential_volume,"label");
  AddNamedReal(lj_potential_volume,"depth");
  AddNamedReal(lj_potential_volume,"rootDist");
  AddNamedReal(lj_potential_volume,"cutOff");
  AddNamedReal(lj_potential_volume,"exvollength");
  AddNamedReal(lj_potential_volume,"beta");

  condlist.push_back(lj_potential_volume);

  /*--------------------------------------------------------------------*/
  // Lennard Jones potential surface

  Teuchos::RCP<ConditionDefinition> lj_potential_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential_Surface",
                                         DRT::Condition::LJ_Potential_Surface,
                                         true,
                                         DRT::Condition::Surface));

  lj_potential_surface->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential_surface,"label");
  AddNamedReal(lj_potential_surface,"depth");
  AddNamedReal(lj_potential_surface,"rootDist");
  AddNamedReal(lj_potential_surface,"cutOff");
  AddNamedReal(lj_potential_surface,"exvollength");
  AddNamedReal(lj_potential_surface,"beta");

  condlist.push_back(lj_potential_surface);


  /*--------------------------------------------------------------------*/
  // Lennard Jones potential line

  Teuchos::RCP<ConditionDefinition> lj_potential_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE LJ_POTENTIAL CONDITIONS",
                                         "Potential",
                                         "LJ_Potential_Line",
                                         DRT::Condition::LJ_Potential_Line,
                                         true,
                                         DRT::Condition::Line));

  lj_potential_line->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(lj_potential_line,"label");
  AddNamedReal(lj_potential_line,"depth");
  AddNamedReal(lj_potential_line,"rootDist");
  AddNamedReal(lj_potential_line,"cutOff");
  AddNamedReal(lj_potential_line,"exvollength");
  AddNamedReal(lj_potential_line,"beta");

  condlist.push_back(lj_potential_line);


  /*--------------------------------------------------------------------*/
  // Van der Waals potential volume
  Teuchos::RCP<ConditionDefinition> vanderwaals_potential_volume =
    Teuchos::rcp(new ConditionDefinition("DESIGN VOL VAN DER WAALS POTENTIAL CONDITIONS",
                                         "Potential",
                                         "VanDerWaals_Potential_Volume",
                                         DRT::Condition::VanDerWaals_Potential_Volume,
                                         true,
                                         DRT::Condition::Volume));

  vanderwaals_potential_volume->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(vanderwaals_potential_volume,"label");
  AddNamedReal(vanderwaals_potential_volume,"lambda");
  AddNamedReal(vanderwaals_potential_volume,"cutOff");
  AddNamedReal(vanderwaals_potential_volume,"beta");
  AddNamedReal(vanderwaals_potential_volume,"exvollength");

  condlist.push_back(vanderwaals_potential_volume);

  /*--------------------------------------------------------------------*/
  // Van der Waals potential surface

  Teuchos::RCP<ConditionDefinition> vanderwaals_potential_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF VAN DER WAALS POTENTIAL CONDITIONS",
                                         "Potential",
                                         "VanDerWaals_Potential_Surface",
                                         DRT::Condition::VanDerWaals_Potential_Surface,
                                         true,
                                         DRT::Condition::Surface));

  vanderwaals_potential_surface->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(vanderwaals_potential_surface,"label");
  AddNamedReal(vanderwaals_potential_surface,"lambda");
  AddNamedReal(vanderwaals_potential_surface,"cutOff");
  AddNamedReal(vanderwaals_potential_surface,"beta");
  AddNamedReal(vanderwaals_potential_surface,"exvollength");

  condlist.push_back(vanderwaals_potential_surface);


  /*--------------------------------------------------------------------*/
  // Van der Waals line

  Teuchos::RCP<ConditionDefinition> vanderwaals_potential_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE VAN DER WAALS POTENTIAL CONDITIONS",
                                         "Potential",
                                         "VanDerWaals_Potential_Line",
                                         DRT::Condition::VanDerWaals_Potential_Line,
                                         true,
                                         DRT::Condition::Line));

  vanderwaals_potential_line->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(vanderwaals_potential_line,"label");
  AddNamedReal(vanderwaals_potential_line,"lambda");
  AddNamedReal(vanderwaals_potential_line,"cutOff");
  AddNamedReal(vanderwaals_potential_line,"beta");
  AddNamedReal(vanderwaals_potential_line,"exvollength");

  condlist.push_back(vanderwaals_potential_line);


  /*-------------------------------------------------------------------*/
  // Electrostatic Repulsion Surface
  Teuchos::RCP<ConditionDefinition> electro_repulsion_potential_surface =
    Teuchos::rcp(new ConditionDefinition("DESIGN SURF ELECTRO REPULSION CONDITIONS",
                                       "Potential",
                                       "ElectroRepulsion_Potential_Surface",
                                       DRT::Condition::ElectroRepulsion_Potential_Surface,
                                       true,
                                       DRT::Condition::Surface));

  electro_repulsion_potential_surface->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(electro_repulsion_potential_surface,"label");
  AddNamedReal(electro_repulsion_potential_surface,"zeta_param_1");
  AddNamedReal(electro_repulsion_potential_surface,"zeta_param_2");
  AddNamedReal(electro_repulsion_potential_surface,"cutOff");
  AddNamedReal(electro_repulsion_potential_surface,"beta");

  condlist.push_back(electro_repulsion_potential_surface);


  /*-------------------------------------------------------------------*/
  // Electrostatic Repulsion Line
  Teuchos::RCP<ConditionDefinition> electro_repulsion_potential_line =
    Teuchos::rcp(new ConditionDefinition("DESIGN LINE ELECTRO REPULSION CONDITIONS",
                                       "Potential",
                                       "ElectroRepulsion_Potential_Line",
                                       DRT::Condition::ElectroRepulsion_Potential_Line,
                                       true,
                                       DRT::Condition::Line));

  electro_repulsion_potential_line->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
  AddNamedInt(electro_repulsion_potential_line,"label");
  AddNamedReal(electro_repulsion_potential_line,"zeta_param_1");
  AddNamedReal(electro_repulsion_potential_line,"zeta_param_2");
  AddNamedReal(electro_repulsion_potential_line,"cutOff");
  AddNamedReal(electro_repulsion_potential_line,"beta");

  condlist.push_back(electro_repulsion_potential_line);


}
