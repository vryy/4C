/*----------------------------------------------------------------------*/
/*!

\brief Input parameters for volmortar

\level 1

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_volmortar.H"



void INPAR::VOLMORTAR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  /* parameters for volmortar */
  Teuchos::ParameterList& volmortar = list->sublist("VOLMORTAR COUPLING", false, "");

  setStringToIntegralParameter<int>("INTTYPE", "Elements", "Type of numerical integration scheme",
      tuple<std::string>("Elements", "elements", "Segments", "segments"),
      tuple<int>(inttype_elements, inttype_elements, inttype_segments, inttype_segments),
      &volmortar);

  setStringToIntegralParameter<int>("COUPLINGTYPE", "Volmortar", "Type of coupling",
      tuple<std::string>("Volmortar", "volmortar", "consistentinterpolation", "consint"),
      tuple<int>(couplingtype_volmortar, couplingtype_volmortar, couplingtype_coninter,
          couplingtype_coninter),
      &volmortar);

  setStringToIntegralParameter<int>("SHAPEFCN", "Dual", "Type of employed set of shape functions",
      tuple<std::string>("Dual", "dual", "Standard", "standard", "std"),
      tuple<int>(shape_dual, shape_dual, shape_std, shape_std, shape_std), &volmortar);

  setStringToIntegralParameter<int>("CUTTYPE", "dd",
      "Type of cut procedure/ integration point calculation",
      tuple<std::string>(
          "dd", "directdivergence", "DirectDivergence", "tessellation", "t", "Tessellation"),
      tuple<int>(cuttype_directdivergence, cuttype_directdivergence, cuttype_directdivergence,
          cuttype_tessellation, cuttype_tessellation, cuttype_tessellation),
      &volmortar);

  setStringToIntegralParameter<int>("DUALQUAD", "nomod",
      "Type of dual shape function for weighting function for quadr. problems",
      tuple<std::string>("nm", "nomod", "lm", "lin_mod", "qm", "quad_mod"),
      tuple<int>(dualquad_no_mod, dualquad_no_mod, dualquad_lin_mod, dualquad_lin_mod,
          dualquad_quad_mod, dualquad_quad_mod),
      &volmortar);

  setStringToIntegralParameter<int>("MESH_INIT", "No",
      "If chosen, mesh initialization procedure is performed", yesnotuple, yesnovalue, &volmortar);

  setStringToIntegralParameter<int>("KEEP_EXTENDEDGHOSTING", "Yes",
      "If chosen, extended ghosting is kept for simulation", yesnotuple, yesnovalue, &volmortar);
}
