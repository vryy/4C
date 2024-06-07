/*----------------------------------------------------------------------*/
/*! \file
 \brief input parameters for porous multiphase problem

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_POROMULTIPHASE_HPP
#define FOUR_C_INPAR_POROMULTIPHASE_HPP


#include "4C_config.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace Inpar
{
  namespace POROMULTIPHASE
  {
    /// Type of coupling strategy for POROMULTIPHASE problems
    enum SolutionSchemeOverFields
    {
      solscheme_undefined,
      solscheme_twoway_partitioned,
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

    //! relaxation methods for partitioned coupling
    enum RelaxationMethods
    {
      relaxation_none,
      relaxation_constant,
      relaxation_aitken
    };

    //! map enum term to std::string
    static inline std::string VectorNormString(const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case Inpar::POROMULTIPHASE::norm_l1:
          return "L1";
          break;
        case Inpar::POROMULTIPHASE::norm_l1_scaled:
          return "L1_scaled";
          break;
        case Inpar::POROMULTIPHASE::norm_l2:
          return "L2";
          break;
        case Inpar::POROMULTIPHASE::norm_rms:
          return "Rms";
          break;
        case Inpar::POROMULTIPHASE::norm_inf:
          return "Inf";
          break;
        default:
          FOUR_C_THROW("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }


    /// set the POROMULTIPHASE parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace POROMULTIPHASE

}  // namespace Inpar



FOUR_C_NAMESPACE_CLOSE

#endif
