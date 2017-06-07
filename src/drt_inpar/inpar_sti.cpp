/*----------------------------------------------------------------------*/
/*!
\file inpar_sti.cpp

\brief input quantities and globally accessible enumerations for scatra-thermo interaction

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "inpar_sti.H"

#include "drt_validparameters.H"
#include "inpar_scatra.H"

#include "../drt_lib/drt_conditiondefinition.H"

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void INPAR::STI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& stidyn = list->sublist(
      "STI DYNAMIC",
      false,
      "general control parameters for scatra-thermo interaction problems"
      );

  // type of scalar transport
  setStringToIntegralParameter<int>(
      "SCATRATYPE",
      "Undefined",
      "type of scalar transport",
       tuple<std::string>(
           "Undefined",
           "ConvectionDiffusion",
           "Elch"
           ),
       tuple<int>(
           INPAR::SCATRA::impltype_undefined,
           INPAR::SCATRA::impltype_std,
           INPAR::SCATRA::impltype_elch_diffcond   // we abuse this enumeration entry here to indicate electrochemistry in general
           ),
       &stidyn
       );

  // specification of initial temperature field
  setStringToIntegralParameter<int>(
      "THERMO_INITIALFIELD",
      "zero_field",
      "initial temperature field for scatra-thermo interaction problems",
      tuple<std::string>(
          "zero_field",
          "field_by_function",
          "field_by_condition"
          ),
      tuple<int>(
          INPAR::SCATRA::initfield_zero_field,
          INPAR::SCATRA::initfield_field_by_function,
          INPAR::SCATRA::initfield_field_by_condition
          ),
      &stidyn
      );

  // function number for initial temperature field
  IntParameter("THERMO_INITFUNCNO",-1,"function number for initial temperature field for scatra-thermo interaction problems",&stidyn);

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<int>(
      "MATRIXTYPE",
      "block",
      "type of global system matrix in global system of equations",
      tuple<std::string>(
          "block",
          "sparse"
          ),
      tuple<int>(
          INPAR::STI::matrix_block,
          INPAR::STI::matrix_sparse
          ),
      &stidyn
      );

  // ID of linear solver for global system of equations
  IntParameter("LINEAR_SOLVER",-1,"ID of linear solver for global system of equations",&stidyn);

  // ID of linear solver for temperature field
  IntParameter("THERMO_LINEAR_SOLVER",-1,"ID of linear solver for temperature field",&stidyn);

  // flag for double condensation of linear equations associated with temperature field
  BoolParameter("THERMO_CONDENSATION","No","flag for double condensation of linear equations associated with temperature field",&stidyn);
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void INPAR::STI::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  return;
}
