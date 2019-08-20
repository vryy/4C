/*--------------------------------------------------------------------------*/
/*! \file
\brief Lubrication dynamic parameters

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_INPAR_INPAR_LUBRICATION_H_
#define SRC_DRT_INPAR_INPAR_LUBRICATION_H_

#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_dserror.H"

namespace INPAR
{
  namespace LUBRICATION
  {
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
      norm_vague = 0,  //!< undetermined norm
      norm_l1,         //!< L1/linear norm
      norm_l2,         //!< L2/Euclidean norm
      norm_rms,        //!< root mean square (RMS) norm
      norm_inf         //!< Maximum/infinity norm
    };

    //! map enum term to std::string
    static inline const std::string VectorNormString(
        const enum VectorNorm norm  //!< input enum term
    )
    {
      switch (norm)
      {
        case INPAR::LUBRICATION::norm_vague:
          return "Vague";
          break;
        case INPAR::LUBRICATION::norm_l1:
          return "L1";
          break;
        case INPAR::LUBRICATION::norm_l2:
          return "L2";
          break;
        case INPAR::LUBRICATION::norm_rms:
          return "Rms";
          break;
        case INPAR::LUBRICATION::norm_inf:
          return "Inf";
          break;
        default:
          dserror("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }

    /// set the lubrication parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace LUBRICATION
}  // namespace INPAR

#endif /* SRC_DRT_INPAR_INPAR_LUBRICATION_H_ */
