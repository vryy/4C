/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for cut library

\level 2


*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_cut.hpp"

#include "4C_cut_enum.hpp"
#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Cut::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& cut_general = list->sublist("CUT GENERAL", false, "");

  // intersection precision (double or cln)
  setStringToIntegralParameter<int>("KERNEL_INTERSECTION_FLOATTYPE", "double",
      "The floattype of the cut surface-edge intersection", tuple<std::string>("cln", "double"),
      tuple<int>(Core::Geo::Cut::floattype_cln, Core::Geo::Cut::floattype_double), &cut_general);

  // Computing disctance surface to point precision (double or cln)
  setStringToIntegralParameter<int>("KERNEL_DISTANCE_FLOATTYPE", "double",
      "The floattype of the cut distance computation", tuple<std::string>("cln", "double"),
      tuple<int>(Core::Geo::Cut::floattype_cln, Core::Geo::Cut::floattype_double), &cut_general);

  // A general floattype for Core::Geo::Cut::Position for Embedded Elements (compute_distance)
  // If specified this floattype is used for all computations of Core::Geo::Cut::Position with
  // embedded elements
  setStringToIntegralParameter<int>("GENERAL_POSITON_DISTANCE_FLOATTYPE", "none",
      "A general floattype for Core::Geo::Cut::Position for Embedded Elements (compute_distance)",
      tuple<std::string>("none", "cln", "double"),
      tuple<int>(Core::Geo::Cut::floattype_none, Core::Geo::Cut::floattype_cln,
          Core::Geo::Cut::floattype_double),
      &cut_general);

  // A general floattype for Core::Geo::Cut::Position for Elements (ComputePosition)
  // If specified this floattype is used for all computations of Core::Geo::Cut::Position
  setStringToIntegralParameter<int>("GENERAL_POSITON_POSITION_FLOATTYPE", "none",
      "A general floattype for Core::Geo::Cut::Position Elements (ComputePosition)",
      tuple<std::string>("none", "cln", "double"),
      tuple<int>(Core::Geo::Cut::floattype_none, Core::Geo::Cut::floattype_cln,
          Core::Geo::Cut::floattype_double),
      &cut_general);

  // Specifiy which Referenceplanes are used in DirectDivergence
  setStringToIntegralParameter<int>("DIRECT_DIVERGENCE_REFPLANE", "all",
      "Specifiy which Referenceplanes are used in DirectDivergence",
      tuple<std::string>("all", "diagonal_side", "facet", "diagonal", "side", "none"),
      tuple<int>(Core::Geo::Cut::DirDiv_refplane_all, Core::Geo::Cut::DirDiv_refplane_diagonal_side,
          Core::Geo::Cut::DirDiv_refplane_facet, Core::Geo::Cut::DirDiv_refplane_diagonal,
          Core::Geo::Cut::DirDiv_refplane_side, Core::Geo::Cut::DirDiv_refplane_none),
      &cut_general);

  // Specifiy is Cutsides are triangulated
  Core::UTILS::BoolParameter(
      "SPLIT_CUTSIDES", "Yes", "Split Quad4 CutSides into Tri3-Subtriangles?", &cut_general);

  // Do the Selfcut before standard CUT
  Core::UTILS::BoolParameter("DO_SELFCUT", "Yes", "Do the SelfCut?", &cut_general);

  // Do meshcorrection in Selfcut
  Core::UTILS::BoolParameter(
      "SELFCUT_DO_MESHCORRECTION", "Yes", "Do meshcorrection in the SelfCut?", &cut_general);

  // Selfcut meshcorrection multiplicator
  Core::UTILS::IntParameter("SELFCUT_MESHCORRECTION_MULTIPLICATOR", 30,
      "ISLANDS with maximal size of the bounding box of h*multiplacator will be removed in the "
      "meshcorrection",
      &cut_general);

  // Cubaturedegree utilized for the numerical integration on the CUT BoundaryCells.
  Core::UTILS::IntParameter("BOUNDARYCELL_CUBATURDEGREE", 20,
      "Cubaturedegree utilized for the numerical integration on the CUT BoundaryCells.",
      &cut_general);

  // Integrate inside volume cells
  Core::UTILS::BoolParameter("INTEGRATE_INSIDE_CELLS", "Yes",
      "Should the integration be done on inside cells", &cut_general);
}

FOUR_C_NAMESPACE_CLOSE
