/*----------------------------------------------------------------------*/
/*! \file
 \brief input parameters for porous multiphase problem with scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_POROMULTIPHASE_SCATRA_HPP
#define FOUR_C_INPAR_POROMULTIPHASE_SCATRA_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}
/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace PoroMultiPhaseScaTra
  {
    /// Type of coupling strategy for poro scatra problems
    enum SolutionSchemeOverFields
    {
      solscheme_twoway_partitioned_nested,
      solscheme_twoway_partitioned_sequential,
      solscheme_twoway_monolithic
    };

    /// type of finite difference check
    enum FdCheck
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
    static inline std::string vector_norm_string(const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case Inpar::PoroMultiPhaseScaTra::norm_l1:
          return "L1";
          break;
        case Inpar::PoroMultiPhaseScaTra::norm_l1_scaled:
          return "L1_scaled";
          break;
        case Inpar::PoroMultiPhaseScaTra::norm_l2:
          return "L2";
          break;
        case Inpar::PoroMultiPhaseScaTra::norm_rms:
          return "Rms";
          break;
        case Inpar::PoroMultiPhaseScaTra::norm_inf:
          return "Inf";
          break;
        default:
          FOUR_C_THROW("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }

    /// set the poromultiphasescatra parameters
    void set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the poromultiphasescatra conditions
    void set_valid_conditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

  }  // namespace PoroMultiPhaseScaTra

}  // namespace Inpar



FOUR_C_NAMESPACE_CLOSE

#endif
