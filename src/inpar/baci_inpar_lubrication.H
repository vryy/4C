/*--------------------------------------------------------------------------*/
/*! \file
\brief Lubrication dynamic parameters

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_INPAR_LUBRICATION_H
#define BACI_INPAR_LUBRICATION_H

#include "baci_config.H"

#include "baci_inpar_parameterlist_utils.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

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

    /// compute velocity by function
    enum VelocityField
    {
      velocity_zero,
      velocity_function,
      velocity_EHL
    };

    /// compute height by function
    enum HeightField
    {
      height_zero,
      height_function,
      height_EHL
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
    static inline std::string VectorNormString(const enum VectorNorm norm  //!< input enum term
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

BACI_NAMESPACE_CLOSE

#endif  // INPAR_LUBRICATION_H
