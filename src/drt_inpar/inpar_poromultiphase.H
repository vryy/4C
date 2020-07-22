/*----------------------------------------------------------------------*/
/*! \file
 \brief input parameters for porous multiphase problem

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_INPAR_INPAR_POROMULTIPHASE_H_
#define SRC_DRT_INPAR_INPAR_POROMULTIPHASE_H_


#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_dserror.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace INPAR
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

    //! relaxation methods for partitioned coupling
    enum RelaxationMethods
    {
      relaxation_none,
      relaxation_constant,
      relaxation_aitken
    };

    //! map enum term to std::string
    static inline const std::string VectorNormString(
        const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case INPAR::POROMULTIPHASE::norm_l1:
          return "L1";
          break;
        case INPAR::POROMULTIPHASE::norm_l1_scaled:
          return "L1_scaled";
          break;
        case INPAR::POROMULTIPHASE::norm_l2:
          return "L2";
          break;
        case INPAR::POROMULTIPHASE::norm_rms:
          return "Rms";
          break;
        case INPAR::POROMULTIPHASE::norm_inf:
          return "Inf";
          break;
        default:
          dserror("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }


    /// set the POROMULTIPHASE parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace POROMULTIPHASE

}  // namespace INPAR



#endif /* SRC_DRT_INPAR_INPAR_POROMULTIPHASE_H_ */
