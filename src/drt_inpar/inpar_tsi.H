/*----------------------------------------------------------------------*/
/*! \file

\brief TSI input parameters
\level 1

*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 11/09 |
 *----------------------------------------------------------------------*/
#ifndef INPAR_TSI_H
#define INPAR_TSI_H

#include "inpar_parameterlist_utils.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |                                                           dano 11/09 |
 *----------------------------------------------------------------------*/
namespace INPAR
{
  namespace TSI
  {
    //! Type of coupling strategy for TSI problems
    enum SolutionSchemeOverFields
    {
      OneWay,
      SequStagg,
      IterStagg,
      IterStaggAitken,
      IterStaggAitkenIrons,
      IterStaggFixedRel,
      Monolithic
    };

    //! Type of coupling strategy for BACI-INCA coupling (TFSI)
    enum BaciIncaCoupling
    {
      TFSI,
      FSI,
      ConjHeatTransfer,
      NoIncaFSI
    };

    //! @name Solution technique and related stuff

    //! type of norm to check for convergence
    enum ConvNorm
    {
      convnorm_abs,  //!< absolute norm
      convnorm_rel,  //!< relative norm of TSI problem with inital TSI rhs
      convnorm_mix   //!< mixed absolute-relative norm
    };

    //! type of norm to check for convergence
    enum BinaryOp
    {
      bop_and,              //!< and
      bop_or,               //!< or
      bop_coupl_or_singl,   //!< either TSI problem or single field problems converged
      bop_coupl_and_singl,  //!< either TSI problem or single field problems converged
      bop_and_singl,        //!< and in single field problems
      bop_or_singl          //!< or in single field problems
    };

    //! type of solution techniques
    enum NlnSolTech
    {
      soltech_newtonfull,  //!< full Newton-Raphson iteration
      soltech_ptc,         //!< pseudo transient continuation nonlinear iteration
    };

    //! type of line-search strategy
    enum LineSearch
    {
      LS_none = 0,   //!< no line search
      LS_structure,  //!< line-search based on structural residual
      LS_thermo,     //!< line-search based on thermal residual
      LS_or,         //!< line-search based on structural or thermal residual
      LS_and,        //!< line-search based on structural and thermal residual
    };

    //! Map solution technique enum to std::string
    static inline const std::string NlnSolTechString(
        const enum NlnSolTech name  //!< enum to convert
    )
    {
      switch (name)
      {
        case soltech_newtonfull:
          return "fullnewton";
          break;
        case soltech_ptc:
          return "ptc";
          break;
        default:
          dserror("Cannot make std::string for solution technique %d", name);
          return "";
      }
    }

    //@}

    //! @name General
    //@{

    //! type of vector norm used for error/residual vectors
    enum VectorNorm
    {
      norm_vague = 0,  //!< undetermined norm
      norm_l1,         //!< L1/linear norm
      norm_l1_scaled,  //!< L1/linear norm scaled by length of vector
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
        case INPAR::TSI::norm_vague:
          return "Vague";
          break;
        case INPAR::TSI::norm_l1:
          return "L1";
          break;
        case INPAR::TSI::norm_l1_scaled:
          return "L1_scaled";
          break;
        case INPAR::TSI::norm_l2:
          return "L2";
          break;
        case INPAR::TSI::norm_rms:
          return "Rms";
          break;
        case INPAR::TSI::norm_inf:
          return "Inf";
          break;
        default:
          dserror("Cannot make std::string to vector norm %d", norm);
          return "";
      }
    }

    //! Method used to calculate plastic dissipation
    enum DissipationMode
    {
      pl_multiplier,  //!< Dissipation = yield stress times plastic multipler
      pl_flow,        //!< Dissipation = Mandel stress : sym(L^p)
      Taylor_Quinney  //!< Dissipation based on Taylor Quinney factor
    };

    //@}

    /// set the tsi parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);


  }  // namespace TSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_TSI_H
