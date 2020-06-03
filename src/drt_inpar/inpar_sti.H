/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-thermo interaction

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_STI_H
#define INPAR_STI_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

namespace INPAR
{
  namespace STI
  {
    //! type of coupling between scatra and thermo fields
    enum class CouplingType
    {
      coupling_undefined,
      coupling_monolithic,
      coupling_oneway_scatratothermo,
      coupling_oneway_thermotoscatra,
      coupling_twoway_scatratothermo,
      coupling_twoway_scatratothermo_aitken,
      coupling_twoway_scatratothermo_aitken_dofsplit,
      coupling_twoway_thermotoscatra,
      coupling_twoway_thermotoscatra_aitken
    };

    //! type of monolithic system matrix in global system of equations
    enum class MatrixType
    {
      matrix_block,
      matrix_sparse
    };

    //! type of scalar transport time integration
    enum class ScaTraTimIntType
    {
      scatratiminttype_standard,
      scatratiminttype_elch
    };

    //! set valid parameters for scatra-thermo interaction
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    //! set valid conditions for scatra-thermo interaction
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace STI
}  // namespace INPAR
#endif  // #ifndef INPAR_STI_H
