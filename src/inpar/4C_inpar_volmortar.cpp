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



void INPAR::VOLMORTAR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /* parameters for volmortar */
  Teuchos::ParameterList& volmortar = list->sublist("VOLMORTAR COUPLING", false, "");

  setStringToIntegralParameter<int>("INTTYPE", "Elements", "Type of numerical integration scheme",
      tuple<std::string>("Elements", "elements", "Segments", "segments"),
      tuple<int>(CORE::VOLMORTAR::inttype_elements, CORE::VOLMORTAR::inttype_elements,
          CORE::VOLMORTAR::inttype_segments, CORE::VOLMORTAR::inttype_segments),
      &volmortar);

  setStringToIntegralParameter<int>("COUPLINGTYPE", "Volmortar", "Type of coupling",
      tuple<std::string>("Volmortar", "volmortar", "consistentinterpolation", "consint"),
      tuple<int>(CORE::VOLMORTAR::couplingtype_volmortar, CORE::VOLMORTAR::couplingtype_volmortar,
          CORE::VOLMORTAR::couplingtype_coninter, CORE::VOLMORTAR::couplingtype_coninter),
      &volmortar);

  setStringToIntegralParameter<int>("SHAPEFCN", "Dual", "Type of employed set of shape functions",
      tuple<std::string>("Dual", "dual", "Standard", "standard", "std"),
      tuple<int>(CORE::VOLMORTAR::shape_dual, CORE::VOLMORTAR::shape_dual,
          CORE::VOLMORTAR::shape_std, CORE::VOLMORTAR::shape_std, CORE::VOLMORTAR::shape_std),
      &volmortar);

  setStringToIntegralParameter<int>("CUTTYPE", "dd",
      "Type of cut procedure/ integration point calculation",
      tuple<std::string>(
          "dd", "directdivergence", "DirectDivergence", "tessellation", "t", "Tessellation"),
      tuple<int>(CORE::VOLMORTAR::cuttype_directdivergence,
          CORE::VOLMORTAR::cuttype_directdivergence, CORE::VOLMORTAR::cuttype_directdivergence,
          CORE::VOLMORTAR::cuttype_tessellation, CORE::VOLMORTAR::cuttype_tessellation,
          CORE::VOLMORTAR::cuttype_tessellation),
      &volmortar);

  setStringToIntegralParameter<int>("DUALQUAD", "nomod",
      "Type of dual shape function for weighting function for quadr. problems",
      tuple<std::string>("nm", "nomod", "lm", "lin_mod", "qm", "quad_mod"),
      tuple<int>(CORE::VOLMORTAR::dualquad_no_mod, CORE::VOLMORTAR::dualquad_no_mod,
          CORE::VOLMORTAR::dualquad_lin_mod, CORE::VOLMORTAR::dualquad_lin_mod,
          CORE::VOLMORTAR::dualquad_quad_mod, CORE::VOLMORTAR::dualquad_quad_mod),
      &volmortar);

  CORE::UTILS::BoolParameter(
      "MESH_INIT", "No", "If chosen, mesh initialization procedure is performed", &volmortar);

  CORE::UTILS::BoolParameter("KEEP_EXTENDEDGHOSTING", "Yes",
      "If chosen, extended ghosting is kept for simulation", &volmortar);
}

FOUR_C_NAMESPACE_CLOSE
