/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for volmortar

\level 1

*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_volmortar.hpp"

#include "4C_coupling_volmortar.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::VolMortar::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /* parameters for volmortar */
  Teuchos::ParameterList& volmortar = list->sublist("VOLMORTAR COUPLING", false, "");

  setStringToIntegralParameter<int>("INTTYPE", "Elements", "Type of numerical integration scheme",
      tuple<std::string>("Elements", "elements", "Segments", "segments"),
      tuple<int>(Coupling::VolMortar::inttype_elements, Coupling::VolMortar::inttype_elements,
          Coupling::VolMortar::inttype_segments, Coupling::VolMortar::inttype_segments),
      &volmortar);

  setStringToIntegralParameter<int>("COUPLINGTYPE", "Volmortar", "Type of coupling",
      tuple<std::string>("Volmortar", "volmortar", "consistentinterpolation", "consint"),
      tuple<int>(Coupling::VolMortar::couplingtype_volmortar,
          Coupling::VolMortar::couplingtype_volmortar, Coupling::VolMortar::couplingtype_coninter,
          Coupling::VolMortar::couplingtype_coninter),
      &volmortar);

  setStringToIntegralParameter<int>("SHAPEFCN", "Dual", "Type of employed set of shape functions",
      tuple<std::string>("Dual", "dual", "Standard", "standard", "std"),
      tuple<int>(Coupling::VolMortar::shape_dual, Coupling::VolMortar::shape_dual,
          Coupling::VolMortar::shape_std, Coupling::VolMortar::shape_std,
          Coupling::VolMortar::shape_std),
      &volmortar);

  setStringToIntegralParameter<int>("CUTTYPE", "dd",
      "Type of cut procedure/ integration point calculation",
      tuple<std::string>(
          "dd", "directdivergence", "DirectDivergence", "tessellation", "t", "Tessellation"),
      tuple<int>(Coupling::VolMortar::cuttype_directdivergence,
          Coupling::VolMortar::cuttype_directdivergence,
          Coupling::VolMortar::cuttype_directdivergence, Coupling::VolMortar::cuttype_tessellation,
          Coupling::VolMortar::cuttype_tessellation, Coupling::VolMortar::cuttype_tessellation),
      &volmortar);

  setStringToIntegralParameter<int>("DUALQUAD", "nomod",
      "Type of dual shape function for weighting function for quadr. problems",
      tuple<std::string>("nm", "nomod", "lm", "lin_mod", "qm", "quad_mod"),
      tuple<int>(Coupling::VolMortar::dualquad_no_mod, Coupling::VolMortar::dualquad_no_mod,
          Coupling::VolMortar::dualquad_lin_mod, Coupling::VolMortar::dualquad_lin_mod,
          Coupling::VolMortar::dualquad_quad_mod, Coupling::VolMortar::dualquad_quad_mod),
      &volmortar);

  Core::UTILS::bool_parameter(
      "MESH_INIT", "No", "If chosen, mesh initialization procedure is performed", &volmortar);

  Core::UTILS::bool_parameter("KEEP_EXTENDEDGHOSTING", "Yes",
      "If chosen, extended ghosting is kept for simulation", &volmortar);
}

FOUR_C_NAMESPACE_CLOSE
