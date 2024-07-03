/*----------------------------------------------------------------------*/
/*! \file
\file inpar_io.cpp

\brief Input parameters for global IO section


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_io.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_inpar_thermo.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::IO::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& io = list->sublist("IO", false, "");

  Core::UTILS::BoolParameter("OUTPUT_GMSH", "No", "", &io);
  Core::UTILS::BoolParameter("OUTPUT_ROT", "No", "", &io);
  Core::UTILS::BoolParameter("OUTPUT_SPRING", "No", "", &io);
  Core::UTILS::BoolParameter("OUTPUT_BIN", "yes", "Do you want to have binary output?", &io);

  // Output every iteration (for debugging purposes)
  Core::UTILS::BoolParameter("OUTPUT_EVERY_ITER", "no",
      "Do you desire structural displ. output every Newton iteration", &io);
  Core::UTILS::IntParameter(
      "OEI_FILE_COUNTER", 0, "Add an output name affix by introducing a additional number", &io);

  // Structural output
  Core::UTILS::BoolParameter("STRUCT_ELE", "Yes", "Output of element properties", &io);
  Core::UTILS::BoolParameter("STRUCT_DISP", "Yes", "Output of displacements", &io);
  Core::UTILS::BoolParameter("STRUCT_VEL_ACC", "No", "Output of velocity and acceleration", &io);
  Core::UTILS::BoolParameter("STRUCT_CURRENT_VOLUME", "No",
      "Output of current element volume as scalar value for each structural element", &io);
  setStringToIntegralParameter<int>("STRUCT_STRESS", "No", "Output of stress",
      tuple<std::string>("No", "no", "NO", "Yes", "yes", "YES", "Cauchy", "cauchy", "2PK", "2pk"),
      tuple<int>(Inpar::Solid::stress_none, Inpar::Solid::stress_none, Inpar::Solid::stress_none,
          Inpar::Solid::stress_2pk, Inpar::Solid::stress_2pk, Inpar::Solid::stress_2pk,
          Inpar::Solid::stress_cauchy, Inpar::Solid::stress_cauchy, Inpar::Solid::stress_2pk,
          Inpar::Solid::stress_2pk),
      &io);
  // in case of a coupled problem (e.g. TSI) the additional stresses are
  // (TSI: thermal stresses) are printed here
  setStringToIntegralParameter<int>("STRUCT_COUPLING_STRESS", "No", "",
      tuple<std::string>("No", "no", "NO", "Yes", "yes", "YES", "Cauchy", "cauchy", "2PK", "2pk"),
      tuple<int>(Inpar::Solid::stress_none, Inpar::Solid::stress_none, Inpar::Solid::stress_none,
          Inpar::Solid::stress_2pk, Inpar::Solid::stress_2pk, Inpar::Solid::stress_2pk,
          Inpar::Solid::stress_cauchy, Inpar::Solid::stress_cauchy, Inpar::Solid::stress_2pk,
          Inpar::Solid::stress_2pk),
      &io);
  setStringToIntegralParameter<int>("STRUCT_STRAIN", "No", "Output of strains",
      tuple<std::string>(
          "No", "no", "NO", "Yes", "yes", "YES", "EA", "ea", "GL", "gl", "LOG", "log"),
      tuple<int>(Inpar::Solid::strain_none, Inpar::Solid::strain_none, Inpar::Solid::strain_none,
          Inpar::Solid::strain_gl, Inpar::Solid::strain_gl, Inpar::Solid::strain_gl,
          Inpar::Solid::strain_ea, Inpar::Solid::strain_ea, Inpar::Solid::strain_gl,
          Inpar::Solid::strain_gl, Inpar::Solid::strain_log, Inpar::Solid::strain_log),
      &io);
  setStringToIntegralParameter<int>("STRUCT_PLASTIC_STRAIN", "No", "",
      tuple<std::string>("No", "no", "NO", "Yes", "yes", "YES", "EA", "ea", "GL", "gl"),
      tuple<int>(Inpar::Solid::strain_none, Inpar::Solid::strain_none, Inpar::Solid::strain_none,
          Inpar::Solid::strain_gl, Inpar::Solid::strain_gl, Inpar::Solid::strain_gl,
          Inpar::Solid::strain_ea, Inpar::Solid::strain_ea, Inpar::Solid::strain_gl,
          Inpar::Solid::strain_gl),
      &io);
  setStringToIntegralParameter<int>("STRUCT_OPTIONAL_QUANTITY", "No",
      "Output of an optional quantity", tuple<std::string>("No", "no", "NO", "membranethickness"),
      tuple<int>(Inpar::Solid::optquantity_none, Inpar::Solid::optquantity_none,
          Inpar::Solid::optquantity_none, Inpar::Solid::optquantity_membranethickness),
      &io);
  Core::UTILS::BoolParameter("STRUCT_SURFACTANT", "No", "", &io);
  Core::UTILS::BoolParameter("STRUCT_JACOBIAN_MATLAB", "No", "", &io);
  setStringToIntegralParameter<Inpar::Solid::ConditionNumber>("STRUCT_CONDITION_NUMBER", "none",
      "Compute the condition number of the structural system matrix and write it to a text file.",
      tuple<std::string>("gmres_estimate", "max_min_ev_ratio", "one-norm", "inf-norm", "none"),
      tuple<Inpar::Solid::ConditionNumber>(Inpar::Solid::ConditionNumber::gmres_estimate,
          Inpar::Solid::ConditionNumber::max_min_ev_ratio, Inpar::Solid::ConditionNumber::one_norm,
          Inpar::Solid::ConditionNumber::inf_norm, Inpar::Solid::ConditionNumber::none),
      &io);
  Core::UTILS::BoolParameter("FLUID_STRESS", "No", "", &io);
  Core::UTILS::BoolParameter("FLUID_WALL_SHEAR_STRESS", "No", "", &io);
  Core::UTILS::BoolParameter("FLUID_ELEDATA_EVRY_STEP", "No", "", &io);
  Core::UTILS::BoolParameter("FLUID_NODEDATA_FIRST_STEP", "No", "", &io);
  Core::UTILS::BoolParameter("THERM_TEMPERATURE", "No", "", &io);
  setStringToIntegralParameter<int>("THERM_HEATFLUX", "None", "",
      tuple<std::string>("None", "No", "NO", "no", "Current", "Initial"),
      tuple<int>(Inpar::THR::heatflux_none, Inpar::THR::heatflux_none, Inpar::THR::heatflux_none,
          Inpar::THR::heatflux_none, Inpar::THR::heatflux_current, Inpar::THR::heatflux_initial),
      &io);
  setStringToIntegralParameter<int>("THERM_TEMPGRAD", "None", "",
      tuple<std::string>("None", "No", "NO", "no", "Current", "Initial"),
      tuple<int>(Inpar::THR::tempgrad_none, Inpar::THR::tempgrad_none, Inpar::THR::tempgrad_none,
          Inpar::THR::tempgrad_none, Inpar::THR::tempgrad_current, Inpar::THR::tempgrad_initial),
      &io);

  Core::UTILS::IntParameter(
      "FILESTEPS", 1000, "Amount of timesteps written to a single result file", &io);
  Core::UTILS::IntParameter("STDOUTEVRY", 1, "Print to screen every n step", &io);

  Core::UTILS::BoolParameter("WRITE_TO_SCREEN", "Yes", "Write screen output", &io);
  Core::UTILS::BoolParameter("WRITE_TO_FILE", "No", "Write the output into a file", &io);

  Core::UTILS::BoolParameter(
      "WRITE_INITIAL_STATE", "yes", "Do you want to write output for initial state ?", &io);
  Core::UTILS::BoolParameter("WRITE_FINAL_STATE", "no",
      "Enforce to write output/restart data at the final state regardless of the other "
      "output/restart intervals",
      &io);

  Core::UTILS::BoolParameter(
      "PREFIX_GROUP_ID", "No", "Put a <GroupID>: in front of every line", &io);
  Core::UTILS::IntParameter(
      "LIMIT_OUTP_TO_PROC", -1, "Only the specified procs will write output", &io);
  setStringToIntegralParameter<int>("VERBOSITY", "verbose", "",
      tuple<std::string>(
          "minimal", "Minimal", "standard", "Standard", "verbose", "Verbose", "debug", "Debug"),
      tuple<int>(FourC::Core::IO::minimal, FourC::Core::IO::minimal, FourC::Core::IO::standard,
          FourC::Core::IO::standard, FourC::Core::IO::verbose, FourC::Core::IO::verbose,
          FourC::Core::IO::debug, FourC::Core::IO::debug),
      &io);

  Core::UTILS::DoubleParameter("RESTARTWALLTIMEINTERVAL", -1.0,
      "Enforce restart after this walltime interval (in seconds), smaller zero to disable", &io);
  Core::UTILS::IntParameter("RESTARTEVRY", -1, "write restart every RESTARTEVRY steps", &io);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& io_every_iter = io.sublist("EVERY ITERATION", false, "");

  // Output every iteration (for debugging purposes)
  Core::UTILS::BoolParameter(
      "OUTPUT_EVERY_ITER", "No", "Do you wish output every Newton iteration?", &io_every_iter);

  Core::UTILS::IntParameter("RUN_NUMBER", -1,
      "Create a new folder for different runs of the same simulation. "
      "If equal -1, no folder is created.",
      &io_every_iter);

  Core::UTILS::IntParameter("STEP_NP_NUMBER", -1,
      "Give the number of the step (i.e. step_{n+1}) for which you want to write the "
      "debug output. If a negative step number is provided, all steps will"
      "be written.",
      &io_every_iter);

  Core::UTILS::BoolParameter("WRITE_OWNER_EACH_NEWTON_ITER", "No",
      "If yes, the ownership "
      "of elements and nodes are written each Newton step, instead of only once"
      "per time/load step.",
      &io_every_iter);
}

FOUR_C_NAMESPACE_CLOSE
