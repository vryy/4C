/*----------------------------------------------------------------------*/
/*!
 \file inpar_porofluidmultiphase.H

 \brief input parameters for porous multiphase fluid problem

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_INPAR_INPAR_POROFLUIDMULTIPHASE_H_
#define SRC_DRT_INPAR_INPAR_POROFLUIDMULTIPHASE_H_

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
  namespace POROFLUIDMULTIPHASE
  {
    /// time integration schemes
    enum TimeIntegrationScheme
    {
      timeint_one_step_theta
    };

    /// compute error compared to analytical solution
    enum CalcError
    {
      calcerror_no,
      calcerror_byfunction
    };

    //! type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_abs,  //!< absolute norm
      convnorm_rel,  //!< relative norm
      convnorm_mix   //!< mixed absolute-relative norm
    };

    //! type of vector norm used for error/residual vectors
    enum VectorNorm
    {
      norm_undefined,
      norm_l1,         //!< L1/linear norm
      norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
      norm_l2,         //!< L2/Euclidean norm
      norm_rms,        //!< root mean square (RMS) norm
      norm_inf         //!< Maximum/infinity norm
    };

    /// type of finite difference check
    enum FDCheck
    {
      fdcheck_none,
      fdcheck_global
    };

    /// initial field for scalar transport problem
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function,
      initfield_field_by_condition
    };

    /// Handling of non-converged nonlinear solver
    enum DivContAct
    {
      divcont_stop,     ///< abort simulation
      divcont_continue  ///< continue nevertheless
    };

    //! reconstruction type of gradients (e.g. velocity gradient)
    enum FluxReconstructionMethod
    {
      gradreco_none,
      // gradreco_spr, super-convergent patch recovery not activated yet
      gradreco_l2
    };

    //! map enum term to std::string
    static inline const std::string VectorNormString(
        const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case INPAR::POROFLUIDMULTIPHASE::norm_l1:
          return "L1";
          break;
        case INPAR::POROFLUIDMULTIPHASE::norm_l1_scaled:
          return "L1_scaled";
          break;
        case INPAR::POROFLUIDMULTIPHASE::norm_l2:
          return "L2";
          break;
        case INPAR::POROFLUIDMULTIPHASE::norm_rms:
          return "Rms";
          break;
        case INPAR::POROFLUIDMULTIPHASE::norm_inf:
          return "Inf";
          break;
        default:
          dserror("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }

    /// set the lubrication parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace POROFLUIDMULTIPHASE

}  // namespace INPAR



#endif /* SRC_DRT_INPAR_INPAR_POROFLUIDMULTIPHASE_H_ */
