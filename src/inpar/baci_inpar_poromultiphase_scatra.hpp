/*----------------------------------------------------------------------*/
/*! \file
 \brief input parameters for porous multiphase problem with scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_POROMULTIPHASE_SCATRA_HPP
#define FOUR_C_INPAR_POROMULTIPHASE_SCATRA_HPP

#include "baci_config.hpp"

#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
{
  namespace POROMULTIPHASESCATRA
  {
    /// Type of coupling strategy for poro scatra problems
    enum SolutionSchemeOverFields
    {
      solscheme_twoway_partitioned_nested,
      solscheme_twoway_partitioned_sequential,
      solscheme_twoway_monolithic
    };

    /// type of finite difference check
    enum FDCheck
    {
      fdcheck_none,
      fdcheck_global
    };

    /// type of norm to be calculated
    enum VectorNorm
    {
      norm_undefined,
      norm_l1,         //!< L1/linear norm
      norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
      norm_l2,         //!< L2/Euclidean norm
      norm_rms,        //!< root mean square (RMS) norm
      norm_inf         //!< Maximum/infinity norm
    };

    //! Handling of non-converged nonlinear solver
    enum DivContAct
    {
      divcont_stop,     ///< abort simulation
      divcont_continue  ///< continue nevertheless
    };

    //! map enum term to std::string
    static inline std::string VectorNormString(const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case INPAR::POROMULTIPHASESCATRA::norm_l1:
          return "L1";
          break;
        case INPAR::POROMULTIPHASESCATRA::norm_l1_scaled:
          return "L1_scaled";
          break;
        case INPAR::POROMULTIPHASESCATRA::norm_l2:
          return "L2";
          break;
        case INPAR::POROMULTIPHASESCATRA::norm_rms:
          return "Rms";
          break;
        case INPAR::POROMULTIPHASESCATRA::norm_inf:
          return "Inf";
          break;
        default:
          dserror("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }

    /// set the poromultiphasescatra parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the poromultiphasescatra conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace POROMULTIPHASESCATRA

}  // namespace INPAR



BACI_NAMESPACE_CLOSE

#endif  // INPAR_POROMULTIPHASE_SCATRA_H
